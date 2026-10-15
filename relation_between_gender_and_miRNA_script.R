# DIFFERENTIAL GENDER ANALYSIS (DAY 2 AND DAY 6)

library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(svglite)
library(dplyr)

output_dir <- "Wyniki_Analizy_Plec"
dir.create(output_dir, showWarnings = FALSE)

dane <- readRDS("00_dane_zintegrowane.rds")
clinical_data <- read.csv("roboczy_kliniczny.csv", sep = ";")

# Define standard legend for gender analysis
legend_levels <- c(
  "Not significant", 
  "Higher expression in Females (Raw P < 0.05)", 
  "Higher expression in Males (Raw P < 0.05)", 
  "Higher expression in Females (FDR < 0.05)", 
  "Higher expression in Males (FDR < 0.05)"
)

legend_colors <- setNames(c("grey80", "#1F78B4", "#E31A1C", "#08306B", "#67000D"), legend_levels)

dummy_df <- data.frame(
  logFC = rep(0, 5),
  P.Value = rep(1, 5), 
  Significance = factor(legend_levels, levels = legend_levels)
)

for(timepoint in c("Day2", "Day6")) {
  
  # Step 1: Filter and merge data
  targets_tp <- subset(dane$targets, Timepoint == timepoint)
  mirna_tp <- dane$counts_filtered[, targets_tp$SampleID]
  
  valid_patient_ids <- as.numeric(targets_tp$PatientID)
  clinical_tp <- clinical_data[clinical_data$Lp %in% valid_patient_ids, ]
  clinical_tp <- clinical_tp[match(valid_patient_ids, clinical_tp$Lp), ]
  
  clinical_tp$Sex <- ifelse(clinical_tp[, 2] == 1, "Male", "Female")
  clinical_tp$Sex <- factor(clinical_tp$Sex, levels = c("Female", "Male"))
  
  # Step 2: Mathematical modeling
  # STATISTICAL JUSTIFICATION: Evaluating gender differences via linear modeling (limma/voom) handles heteroscedasticity across the transcriptome better than standalone non-parametric tests.
  design <- model.matrix(~ 0 + Sex, data = clinical_tp)
  colnames(design) <- levels(clinical_tp$Sex)
  
  y <- DGEList(counts = mirna_tp)
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design, plot = FALSE)
  
  fit <- lmFit(v, design)
  contrast_matrix <- makeContrasts(Male_vs_Female = Male - Female, levels = design)
  
  # BIOINFORMATIC JUSTIFICATION (EBAYES): Empirical Bayes moderation is crucial here due to the small sample sizes within gender subgroups. It leverages global variance to stabilize gene-wise estimates.
  fit_ebayes <- eBayes(contrasts.fit(fit, contrast_matrix))
  
  # Export Voom metrics and filtered counts
  write.csv(v$E, file.path(output_dir, paste0("Macierz_Voom_Expression_", timepoint, ".csv")), row.names = TRUE)
  write.csv(v$weights, file.path(output_dir, paste0("Macierz_Voom_Weights_", timepoint, ".csv")), row.names = TRUE)
  write.csv(mirna_tp, file.path(output_dir, paste0("Macierz_Counts_Filtered_", timepoint, ".csv")), row.names = TRUE)
  
  res_full <- topTable(fit_ebayes, coef = "Male_vs_Female", number = Inf, sort.by = "P")
  write.csv(res_full, file.path(output_dir, paste0("Wyniki_Plec_Pelne_", timepoint, ".csv")), row.names = TRUE)
  
  # Step 3: Significant results table
  significant_table <- res_full %>%
    filter(P.Value < 0.05 & abs(logFC) > 1) %>%
    mutate(
      miRNA = rownames(.),
      Expression_Trend = ifelse(logFC > 0, "Higher in Males", "Higher in Females"),
      logFC = sprintf("%.4f", logFC),
      P.Value = sprintf("%.4f", P.Value),
      adj.P.Val = sprintf("%.4f", adj.P.Val)
    ) %>%
    select(miRNA, Expression_Trend, logFC, P.Value, adj.P.Val)
  
  write.csv(significant_table, file.path(output_dir, paste0("Tabela_Slajd_Plec_TYLKO_ISTOTNE_", timepoint, ".csv")), row.names = FALSE)
  
  # Step 4: Classify points for Volcano plot
  res_full$Significance <- legend_levels[1]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC < -1] <- legend_levels[2]
  res_full$Significance[res_full$P.Value < 0.05 & res_full$logFC >  1] <- legend_levels[3]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC < -1] <- legend_levels[4]
  res_full$Significance[res_full$adj.P.Val < 0.05 & res_full$logFC >  1] <- legend_levels[5]
  res_full$Significance <- factor(res_full$Significance, levels = legend_levels)
  
  res_full$Label <- ""
  sig_idx <- which(res_full$P.Value < 0.05 & abs(res_full$logFC) > 1)
  res_full$Label[sig_idx] <- rownames(res_full)[sig_idx]
  
  # Step 5: Volcano plot execution
  p_volc <- ggplot(res_full, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.85, size = 2.5, show.legend = FALSE) +
    geom_point(data = dummy_df, aes(x = logFC, y = -log10(P.Value), color = Significance), alpha = 0, size = 2.5, show.legend = TRUE) +
    scale_color_manual(values = legend_colors, drop = FALSE) +
    guides(color = guide_legend(title = paste("Expression (Gender) -", timepoint), override.aes = list(shape = 16, color = unname(legend_colors), size = 4, alpha = 1))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_text_repel(data = subset(res_full, Label != ""), aes(label = Label), size = 3.2, max.overlaps = Inf, box.padding = 0.6, point.padding = 0.4, force = 4, segment.color = "grey50", segment.size = 0.3, show.legend = FALSE) +
    theme_minimal(base_size = 12) +
    labs(x = "Log2 Fold Change (Males vs Females)", y = "-Log10(Raw P-value)") +
    theme(plot.title = element_blank(), plot.subtitle = element_blank(), legend.position = "right")
  
  ggsave(file.path(output_dir, paste0("Volcano_Plec_", timepoint, ".png")), plot = p_volc, width = 8.5, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0("Volcano_Plec_", timepoint, ".svg")), plot = p_volc, width = 8.5, height = 6)
}

cat("\nGender analysis complete.\n")