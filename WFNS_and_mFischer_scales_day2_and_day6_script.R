# SPEARMAN RANK CORRELATION AGAINST CLINICAL SCALES

library(ggplot2)
library(svglite)
library(dplyr)
library(tidyr)

dane <- readRDS("00_dane_zintegrowane.rds")
cpm_TMM <- dane$cpm_TMM
targets <- dane$targets

# Export log-transformed counts log2(CPM + 1) used for plotting
# BIOINFORMATIC JUSTIFICATION: For visualizations, log2(CPM + 1) stabilizes the variance across the expression spectrum, avoiding heavily skewed, unintelligible scatterplots.
output_dir_matrices <- "Wyniki_Korelacje_Macierze"
dir.create(output_dir_matrices, showWarnings = FALSE)
log2_cpm_plus_1 <- log2(cpm_TMM + 1)
write.csv(cpm_TMM, file.path(output_dir_matrices, "Macierz_CPM_TMM.csv"), row.names = TRUE)
write.csv(log2_cpm_plus_1, file.path(output_dir_matrices, "Macierz_log2_CPM_plus_1.csv"), row.names = TRUE)

clinical_data <- read.csv("roboczy_kliniczny.csv", sep = ";")
clinical_data$WFNS <- as.numeric(gsub(",", ".", as.character(clinical_data[,6])))
clinical_data$mFisher <- as.numeric(gsub(",", ".", as.character(clinical_data$mFisher)))

# Step 1: Calculate temporal expression deltas
cpm_df <- as.data.frame(cpm_TMM)
cpm_df$miRNA <- rownames(cpm_df)

targets_long <- targets %>% select(SampleID, PatientID, Group, Timepoint)
cpm_long <- cpm_df %>% 
  pivot_longer(cols = -miRNA, names_to = "SampleID", values_to = "CPM") %>%
  left_join(targets_long, by = "SampleID")

cpm_wide <- cpm_long %>%
  select(miRNA, PatientID, Group, Timepoint, CPM) %>%
  pivot_wider(names_from = Timepoint, values_from = CPM) %>%
  na.omit() %>%
  mutate(Delta = Day6 - Day2)

deltas_summary <- cpm_wide %>%
  group_by(miRNA) %>%
  summarise(
    Mean_Delta_DCI = mean(Delta[Group == "DCI"], na.rm = TRUE),
    Mean_Delta_nonDCI = mean(Delta[Group == "nonDCI"], na.rm = TRUE),
    Delta_All = mean(Delta, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Iterate through clinical scales and compute Spearman correlations
format_table <- function(df) {
  df$miRNA <- gsub("\\.", "-", df$miRNA)
  df$Rho <- sprintf("%.4f", df$Rho)
  df$P_Value <- sprintf("%.5f", df$P_Value)
  df$FDR <- sprintf("%.4f", df$FDR)
  if("Mean_Delta_DCI" %in% names(df)) {
    df$Mean_Delta_DCI <- sprintf("%.3f", df$Mean_Delta_DCI)
    df$Mean_Delta_nonDCI <- sprintf("%.3f", df$Mean_Delta_nonDCI)
    df$Delta_All <- sprintf("%.3f", df$Delta_All)
  }
  return(df)
}

for (scale_name in c("WFNS", "mFisher")) {
  for (timepoint in c("Day2", "Day6")) {
    
    targets_timepoint <- subset(targets, Timepoint == timepoint)
    cpm_timepoint <- cpm_TMM[, targets_timepoint$SampleID]
    
    valid_patient_ids <- as.numeric(targets_timepoint$PatientID)
    clinical_tp <- clinical_data[clinical_data$Lp %in% valid_patient_ids, ]
    clinical_tp <- clinical_tp[match(valid_patient_ids, clinical_tp$Lp), ]
    
    cor_results <- data.frame(miRNA = rownames(cpm_timepoint), Rho = NA, P_Value = NA)
    
    for (i in seq_len(nrow(cpm_timepoint))) {
      if (sum(cpm_timepoint[i, ]) > 0) {
        # STATISTICAL JUSTIFICATION: Spearman's rank correlation is chosen over Pearson's because WFNS and mFisher are ordinal, non-parametric severity scales, violating the normality and continuous interval assumptions required for Pearson.
        test <- cor.test(clinical_tp[[scale_name]], cpm_timepoint[i, ], method = "spearman", exact = FALSE)
        cor_results$Rho[i] <- test$estimate
        cor_results$P_Value[i] <- test$p.value
      }
    }
    
    cor_results <- na.omit(cor_results)
    
    # STATISTICAL JUSTIFICATION: Benjamini-Hochberg (BH) False Discovery Rate (FDR) adjustment corrects for multiple hypothesis testing across the entire transcriptome, limiting type I errors.
    cor_results$FDR <- p.adjust(cor_results$P_Value, method = "BH")
    
    cor_results <- merge(cor_results, deltas_summary, by = "miRNA")
    
    res_pos <- head(cor_results[cor_results$Rho > 0, ][order(cor_results[cor_results$Rho > 0, ]$P_Value), ], 15)
    res_neg <- head(cor_results[cor_results$Rho < 0, ][order(cor_results[cor_results$Rho < 0, ]$P_Value), ], 15)
    
    cor_dir <- paste0("Wyniki_Korelacje_", scale_name, "_", timepoint)
    dir.create(cor_dir, showWarnings = FALSE)
    
    write.csv(format_table(res_pos), file.path(cor_dir, "Top15_Positive_with_Deltas.csv"), row.names = FALSE)
    write.csv(format_table(res_neg), file.path(cor_dir, "Top15_Negative_with_Deltas.csv"), row.names = FALSE)
    
    generate_plot <- function(row_data, direction, rank_idx) {
      mirna_id <- row_data$miRNA
      plot_df <- data.frame(
        Scale = clinical_tp[[scale_name]], 
        LogCPM = log2_cpm_plus_1[mirna_id, targets_timepoint$SampleID]
      )
      
      mirna_plot_name <- gsub("\\.", "-", mirna_id)
      
      # STATISTICAL JUSTIFICATION: Trendlines (geom_smooth) were intentionally omitted. Drawing linear models across strictly ordinal x-axis categories (1, 2, 3, 4, 5) misrepresents the non-parametric nature of the test.
      p <- ggplot(plot_df, aes(x = Scale, y = LogCPM)) +
        geom_point(color = ifelse(direction == "Positive", "firebrick", "royalblue"), size = 3, alpha = 0.85) +
        theme_minimal(base_size = 13) +
        labs(title = mirna_plot_name, 
             subtitle = paste0("Rho = ", sprintf("%.3f", row_data$Rho), " | P = ", sprintf("%.3g", row_data$P_Value)),
             x = paste("Score", scale_name, "(", timepoint, ")"), 
             y = paste(mirna_plot_name, "log2(CPM + 1)"))
      
      safe_name <- gsub("[^A-Za-z0-9_-]", "_", mirna_plot_name)
      ggsave(file.path(cor_dir, paste0("Top", rank_idx, "_", direction, "_", safe_name, ".png")), plot = p, width = 7, height = 5, dpi = 300)
      ggsave(file.path(cor_dir, paste0("Top", rank_idx, "_", direction, "_", safe_name, ".svg")), plot = p, width = 7, height = 5)
    }
    
    for (i in 1:min(5, nrow(res_pos))) generate_plot(res_pos[i, ], "Positive", i)
    for (i in 1:min(5, nrow(res_neg))) generate_plot(res_neg[i, ], "Negative", i)
  }
}

cat("Correlations generated successfully.\n")