# COMPREHENSIVE MIRNA ANALYSIS MASTER SCRIPT (aSAH / DCI)
# INCLUDES 4 CONTRASTS, MEAN/SD CALCULATIONS, FORMATTED TABLES, VOLCANO PLOTS

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(svglite)
library(dplyr)

# Step 1: Environment and data preparation
output_dir <- "Wyniki_Analizy_Master_Dynamika"
dir.create(output_dir, showWarnings = FALSE)

dane <- readRDS("00_dane_zintegrowane.rds")
counts_filtered <- dane$counts_filtered
logTMM <- dane$logTMM 
targets <- dane$targets

# Step 2: Descriptive statistics calculation
calculate_stats <- function(matrix_data, sample_ids) {
  sub_data <- matrix_data[, sample_ids, drop=FALSE]
  list(Mean = rowMeans(sub_data), SD = apply(sub_data, 1, sd))
}

idx_DCI_d2    <- targets$SampleID[targets$Group == "DCI" & targets$Timepoint == "Day2"]
idx_DCI_d6    <- targets$SampleID[targets$Group == "DCI" & targets$Timepoint == "Day6"]
idx_nonDCI_d2 <- targets$SampleID[targets$Group == "nonDCI" & targets$Timepoint == "Day2"]
idx_nonDCI_d6 <- targets$SampleID[targets$Group == "nonDCI" & targets$Timepoint == "Day6"]

stat_DCI_d2    <- calculate_stats(logTMM, idx_DCI_d2)
stat_DCI_d6    <- calculate_stats(logTMM, idx_DCI_d6)
stat_nonDCI_d2 <- calculate_stats(logTMM, idx_nonDCI_d2)
stat_nonDCI_d6 <- calculate_stats(logTMM, idx_nonDCI_d6)

stats_table <- data.frame(
  miRNA = rownames(logTMM),
  Mean_DCI_D2 = stat_DCI_d2$Mean,       SD_DCI_D2 = stat_DCI_d2$SD,
  Mean_nonDCI_D2 = stat_nonDCI_d2$Mean, SD_nonDCI_D2 = stat_nonDCI_d2$SD,
  Mean_DCI_D6 = stat_DCI_d6$Mean,       SD_DCI_D6 = stat_DCI_d6$SD,
  Mean_nonDCI_D6 = stat_nonDCI_d6$Mean, SD_nonDCI_D6 = stat_nonDCI_d6$SD
)

# Step 3: Paired Limma model execution
# STATISTICAL JUSTIFICATION: DGEList initialization and TMM normalization apply composition adjustments directly to the raw counts, essential before linear modeling.
y <- DGEList(counts = counts_filtered, samples = targets)
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ 0 + Group_Time, data = targets)
colnames(design) <- levels(targets$Group_Time)

# BIOINFORMATIC JUSTIFICATION (VOOM): The voom transformation computes precision weights based on the mean-variance trend of the count data. 
# This renders the counts suitable for linear modeling (which assumes constant variance) by down-weighting highly variable observations.
v_pre <- voom(y, design, plot = FALSE)

# STATISTICAL JUSTIFICATION (DUPLICATE CORRELATION): Because Day 2 and Day 6 measurements originate from the same patients, they are not statistically independent. 
# Estimating the consensus correlation models the intra-patient pairing, preventing false positives caused by baseline individual biological differences.
corfit <- duplicateCorrelation(v_pre, design, block = targets$PatientID)
v <- voom(y, design, plot = FALSE, block = targets$PatientID, correlation = corfit$consensus)
fit <- lmFit(v, design, block = targets$PatientID, correlation = corfit$consensus)

# Export filtered counts and Voom transformation matrices
write.csv(counts_filtered, file.path(output_dir, "Macierz_Counts_Filtered.csv"), row.names = TRUE)
write.csv(v$E, file.path(output_dir, "Macierz_Voom_Expression.csv"), row.names = TRUE)
write.csv(v$weights, file.path(output_dir, "Macierz_Voom_Weights.csv"), row.names = TRUE)

# STATISTICAL JUSTIFICATION: Explicit contrasts isolate the specific biological questions.
contrast_matrix <- makeContrasts(
  Doba2_DCI_vs_nonDCI   = DCI_Day2 - nonDCI_Day2,
  Doba6_DCI_vs_nonDCI   = DCI_Day6 - nonDCI_Day6,
  Dynamika_Interakcja   = (DCI_Day6 - DCI_Day2) - (nonDCI_Day6 - nonDCI_Day2),
  Czas_Ogolna_Zmiana    = ((DCI_Day6 + nonDCI_Day6)/2) - ((DCI_Day2 + nonDCI_Day2)/2),
  levels = design
)

# BIOINFORMATIC JUSTIFICATION (EBAYES): Empirical Bayes moderation squeezes the gene-wise residual variances towards a common global variance. 
# This protects against false positives arising from artificially small sample variances, which is crucial in low-N transcriptomics.
fit_ebayes <- eBayes(contrasts.fit(fit, contrast_matrix))

# Step 4: Export helper functions
format_comma_decimal <- function(df) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) format(x, scientific = FALSE, dec = ","))
  df
}

prepare_slide_table <- function(df_results, file_name, contrast_type) {
  tabela <- df_results %>%
    filter(P.Value < 0.05 & abs(logFC) > 1) %>%
    mutate(
      Expression_Trend = case_when(
        contrast_type == "Czas_Ogolna_Zmiana" & logFC > 0 ~ "Upregulated on Day 6",
        contrast_type == "Czas_Ogolna_Zmiana" & logFC < 0 ~ "Downregulated on Day 6",
        contrast_type == "Dynamika_Interakcja" & logFC > 0 ~ "Relative increase in DCI",
        contrast_type == "Dynamika_Interakcja" & logFC < 0 ~ "Relative decrease in DCI",
        TRUE & logFC > 0 ~ "Higher expression in DCI",
        TRUE & logFC < 0 ~ "Lower expression in DCI"
      ),
      logFC = sprintf("%.4f", logFC),
      P.Value = sprintf("%.4f", P.Value),
      adj.P.Val = sprintf("%.4f", adj.P.Val) # FDR correction via Benjamini-Hochberg applied automatically by limma
    ) %>%
    select(miRNA, Expression_Trend, logFC, P.Value, adj.P.Val)
  
  write.csv(tabela, file = file.path(output_dir, file_name), row.names = FALSE)
}

# Step 5: Main analysis and visualization loop
for (contrast in colnames(contrast_matrix)) {
  res_limma <- topTable(fit_ebayes, coef = contrast, number = Inf, sort.by = "P")
  res_limma$miRNA <- rownames(res_limma)
  res_full <- merge(stats_table, res_limma, by = "miRNA")
  res_full <- res_full[order(res_full$P.Value), ]
  rownames(res_full) <- res_full$miRNA
  
  write.csv(format_comma_decimal(res_full), file.path(output_dir, paste0("SZEF_Pelne_Wyniki_i_Statystyki_", contrast, ".csv")), row.names = FALSE)
  prepare_slide_table(res_full, paste0("SLAJD_Tylko_Istotne_", contrast, ".csv"), contrast)
  
  if (contrast == "Czas_Ogolna_Zmiana") {
    legend_title <- "Expression (Day 6 vs Day 2)"
    x_axis_label <- "Log2 Fold Change (Day 6 vs Day 2)"
    leg_levels <- c("Not significant", "Downregulated Day 6 (Raw P < 0.05)", "Upregulated Day 6 (Raw P < 0.05)", "Downregulated Day 6 (FDR < 0.05)", "Upregulated Day 6 (FDR < 0.05)")
  } else if (contrast == "Dynamika_Interakcja") {
    legend_title <- "Dynamics (DCI vs nonDCI)"
    x_axis_label <- "Log2 Fold Change (DCI trajectory shift)"
    leg_levels <- c("Not significant", "Relative decrease in DCI (Raw P < 0.05)", "Relative increase in DCI (Raw P < 0.05)", "Relative decrease in DCI (FDR < 0.05)", "Relative increase in DCI (FDR < 0.05)")
  } else {
    legend_title <- "Expression (DCI vs nonDCI)"
    x_axis_label <- "Log2 Fold Change (DCI vs nonDCI)"
    leg_levels <- c("Not significant", "Lower expression in DCI (Raw P < 0.05)", "Higher expression in DCI (Raw P < 0.05)", "Lower expression in DCI (FDR < 0.05)", "Higher expression in DCI (FDR < 0.05)")
  }
  
  legend_colors <- setNames(c("grey80", "#1F78B4", "#E31A1C", "#08306B", "#67000D"), leg_levels)
  
  res_full$Significance <- leg_levels[1]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC < -1] <- leg_levels[2]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC >  1] <- leg_levels[3]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC < -1] <- leg_levels[4]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC >  1] <- leg_levels[5]
  res_full$Significance <- factor(res_full$Significance, levels = leg_levels)
  
  res_full$Label <- ""
  sig_idx <- which(res_full$P.Value < 0.05 & abs(res_full$logFC) > 1)
  res_full$Label[sig_idx] <- rownames(res_full)[sig_idx]
  
  dummy_df <- data.frame(logFC = rep(0, 5), P.Value = rep(1, 5), Significance = factor(leg_levels, levels = leg_levels))
  
  p_volc <- ggplot(res_full, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.85, size = 2.5, show.legend = FALSE) +
    geom_point(data = dummy_df, aes(x = logFC, y = -log10(P.Value), color = Significance), alpha = 0, size = 2.5, show.legend = TRUE) +
    scale_color_manual(values = legend_colors, drop = FALSE) +
    guides(color = guide_legend(title = legend_title, override.aes = list(shape = 16, color = unname(legend_colors), size = 4, alpha = 1))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_text_repel(data = subset(res_full, Label != ""), aes(label = Label), size = 3.2, max.overlaps = Inf, box.padding = 0.6, point.padding = 0.4, force = 4, segment.color = "grey50", segment.size = 0.3, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    labs(x = x_axis_label, y = "-Log10(Raw P-value)") +
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("Volcano_", contrast, ".png")), plot = p_volc, width = 9.5, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0("Volcano_", contrast, ".svg")), plot = p_volc, width = 9.5, height = 6)
}

cat("\nAnalysis complete. All tables and plots are available in directory:", output_dir, "\n")