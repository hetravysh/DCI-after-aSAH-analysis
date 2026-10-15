# TWO-FACTOR LIMMA MODEL (DCI * TIMEPOINT) WITH DUPLICATE CORRELATION
# COMPREHENSIVE OUTPUT: 4 CONTRASTS, SIGNIFICANCE TABLES, VOLCANO PLOTS

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(svglite)
library(dplyr)

# Step 1: Initialize output directory
output_dir <- "Wyniki_Analizy_Kontrasty_Osobne"
dir.create(output_dir, showWarnings = FALSE)

# Step 2: Load integrated datasets
dane <- readRDS("00_dane_zintegrowane.rds")
counts_filtered <- dane$counts_filtered
targets <- dane$targets

# Step 3: Construct experimental design matrix
design <- model.matrix(~ 0 + Group_Time, data = targets)
colnames(design) <- levels(targets$Group_Time)

# Step 4: TMM normalization and intra-patient correlation estimation
# STATISTICAL JUSTIFICATION: DGEList initialization and TMM normalization apply composition adjustments directly to the raw counts, essential before linear modeling.
y <- DGEList(counts = counts_filtered)
y <- calcNormFactors(y, method = "TMM")

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

# Step 5: Define 4 primary clinical contrasts
# STATISTICAL JUSTIFICATION: Explicit contrasts isolate the specific biological questions.
contrast_matrix <- makeContrasts(
  DCI_vs_nonDCI_Day2  = DCI_Day2 - nonDCI_Day2,
  DCI_vs_nonDCI_Day6  = DCI_Day6 - nonDCI_Day6,
  Interakcja_Dynamika = (DCI_Day6 - DCI_Day2) - (nonDCI_Day6 - nonDCI_Day2),
  Ogolna_Zmiana_Czasowa = ((DCI_Day6 + nonDCI_Day6)/2) - ((DCI_Day2 + nonDCI_Day2)/2),
  levels = design
)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# BIOINFORMATIC JUSTIFICATION (EBAYES): Empirical Bayes moderation squeezes the gene-wise residual variances towards a common global variance. 
# This protects against false positives arising from artificially small sample variances, which is crucial in low-N transcriptomics.
fit_ebayes <- eBayes(fit_contrasts)

# Step 6: Presentation-ready table generator
prepare_slide_table <- function(fit_obj, coef_name, file_name) {
  res <- topTable(fit_obj, coef = coef_name, number = Inf, sort.by = "P")
  
  tabela <- res %>%
    filter(P.Value < 0.05 & abs(logFC) > 1) %>%
    mutate(
      miRNA = rownames(.),
      Expression_Trend = case_when(
        coef_name == "Ogolna_Zmiana_Czasowa" & logFC > 0 ~ "Upregulated on Day 6",
        coef_name == "Ogolna_Zmiana_Czasowa" & logFC < 0 ~ "Downregulated on Day 6",
        coef_name == "Interakcja_Dynamika"   & logFC > 0 ~ "Relative increase in DCI",
        coef_name == "Interakcja_Dynamika"   & logFC < 0 ~ "Relative decrease in DCI",
        TRUE & logFC > 0 ~ "Higher expression in DCI",
        TRUE & logFC < 0 ~ "Lower expression in DCI"
      ),
      logFC = sprintf("%.4f", logFC),
      P.Value = sprintf("%.4f", P.Value),
      adj.P.Val = sprintf("%.4f", adj.P.Val) # FDR correction via Benjamini-Hochberg applied automatically by limma
    ) %>%
    select(miRNA, Expression_Trend, logFC, P.Value, adj.P.Val)
  
  write.csv(tabela, file = file.path(output_dir, file_name), row.names = FALSE)
  return(tabela)
}

# Automatically generate significant results tables
for (contrast in colnames(contrast_matrix)) {
  prepare_slide_table(fit_ebayes, contrast, paste0("Tabela_Slajd_TYLKO_ISTOTNE_", contrast, ".csv"))
}

# Step 7: Loop to generate complete results CSVs and Volcano plots
for (contrast in colnames(contrast_matrix)) {
  res_full <- topTable(fit_ebayes, coef = contrast, number = Inf, sort.by = "P")
  write.csv(res_full, file = file.path(output_dir, paste0("Cale_Wyniki_Model_", contrast, ".csv")), row.names = TRUE)
  
  # Dynamic legend label assignment
  if (contrast == "Ogolna_Zmiana_Czasowa") {
    legend_title <- "Expression (Day 6 vs Day 2)"
    x_axis_label <- "Log2 Fold Change (Day 6 vs Day 2)"
    legend_levels <- c("Not significant", "Downregulated Day 6 (Raw P < 0.05)", "Upregulated Day 6 (Raw P < 0.05)", "Downregulated Day 6 (FDR < 0.05)", "Upregulated Day 6 (FDR < 0.05)")
  } else if (contrast == "Interakcja_Dynamika") {
    legend_title <- "Relative Dynamics (DCI vs nonDCI)"
    x_axis_label <- "Log2 Fold Change (Trajectory shift)"
    legend_levels <- c("Not significant", "Relative decrease in DCI (Raw P < 0.05)", "Relative increase in DCI (Raw P < 0.05)", "Relative decrease in DCI (FDR < 0.05)", "Relative increase in DCI (FDR < 0.05)")
  } else {
    legend_title <- "Expression (DCI vs nonDCI)"
    x_axis_label <- "Log2 Fold Change (DCI vs nonDCI)"
    legend_levels <- c("Not significant", "Lower expression in DCI (Raw P < 0.05)", "Higher expression in DCI (Raw P < 0.05)", "Lower expression in DCI (FDR < 0.05)", "Higher expression in DCI (FDR < 0.05)")
  }
  
  legend_colors <- setNames(c("grey80", "#1F78B4", "#E31A1C", "#08306B", "#67000D"), legend_levels)
  
  res_full$Significance <- legend_levels[1]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC < -1] <- legend_levels[2]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC >  1] <- legend_levels[3]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC < -1] <- legend_levels[4]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC >  1] <- legend_levels[5]
  res_full$Significance <- factor(res_full$Significance, levels = legend_levels)
  
  res_full$Label <- ""
  sig_idx <- which(res_full$P.Value < 0.05 & abs(res_full$logFC) > 1)
  res_full$Label[sig_idx] <- rownames(res_full)[sig_idx]
  
  dummy_df <- data.frame(logFC = rep(0, 5), P.Value = rep(1, 5), Significance = factor(legend_levels, levels = legend_levels))
  
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
  
  ggsave(file.path(output_dir, paste0("Volcano_", contrast, ".png")), plot = p_volc, width = 8.5, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0("Volcano_", contrast, ".svg")), plot = p_volc, width = 8.5, height = 6)
}

cat("\nAnalysis complete. Filtered tables and plots generated.\n")