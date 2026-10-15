# EXPLORATORY HEATMAP AND Z-SCORE STANDARDIZATION
# Top 50 most variable miRNAs in SVG Format

library(pheatmap)

output_dir <- "Wyniki_Analizy_Heatmap"
dir.create(output_dir, showWarnings = FALSE)

# Load integrated dataset
dane <- readRDS("00_dane_zintegrowane.rds")
logTMM <- dane$logTMM
targets <- dane$targets

# Remove "hsa" prefix for visualization clarity
rownames(logTMM) <- gsub("^hsa[.-]", "", rownames(logTMM), ignore.case = TRUE)

# Export cleaned logTMM matrix
write.csv(logTMM, file.path(output_dir, "Macierz_logTMM_Czyste_Nazwy.csv"), row.names = TRUE)

# Export Z-score standardized matrix (per row/miRNA scaling)
zscore_matrix <- t(scale(t(logTMM)))
write.csv(zscore_matrix, file.path(output_dir, "Macierz_Zscore_Standaryzowana.csv"), row.names = TRUE)

# Calculate variance and extract top 50
# BIOINFORMATIC JUSTIFICATION: Subsetting the top 50 most variable miRNAs isolates the most biologically dynamic signals across the cohort, filtering out constitutive 'housekeeping' genes that provide no discriminatory power.
miRNA_variance <- apply(logTMM, 1, var)
top_50_miRNAs <- names(sort(miRNA_variance, decreasing = TRUE))[1:50]
heatmap_data_horizontal <- t(logTMM[top_50_miRNAs, ])

# Construct annotation data frame
annotation_row <- data.frame(
  DCI = ifelse(targets$Group == "DCI", "DCI", "Non-DCI"),
  Timepoint = ifelse(targets$Timepoint == "Day2", "Day 2", "Day 6"),
  row.names = targets$SampleID
)

# Define legend color mapping
annotation_colors <- list(
  DCI = c("DCI" = "#b2182b", "Non-DCI" = "#2166ac"),
  Timepoint = c("Day 2" = "#e0e0e0", "Day 6" = "#4d4d4d")
)

# Initialize SVG graphics device
svg(file.path(output_dir, "Exploratory_Heatmap_Top50_Pozioma_Poprawiona.svg"), width = 12, height = 8)

# STATISTICAL JUSTIFICATION: Scale="column" automatically computes Z-scores across each miRNA profile. 
# This normalizes absolute expression differences, centering the mean at 0 and standard deviation at 1, allowing clear visual comparison of relative up/down-regulation trends irrespective of basal abundance.
pheatmap(
  mat = heatmap_data_horizontal,
  scale = "column",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8
)

dev.off()
cat("\nHeatmap generated successfully. LogTMM and Z-score matrices exported to CSV.\n")