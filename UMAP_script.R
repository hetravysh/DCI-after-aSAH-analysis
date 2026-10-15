# DIMENSIONALITY REDUCTION: UMAP QUALITY CONTROL ON NORMALIZED logTMM DATA

library(ggplot2)
library(umap)
library(ggrepel)

output_dir <- "Wyniki_Analizy_UMAP"
dir.create(output_dir, showWarnings = FALSE)

dane <- readRDS("00_dane_zintegrowane.rds")
logTMM <- dane$logTMM
targets <- dane$targets

# BIOINFORMATIC JUSTIFICATION (UMAP): Uniform Manifold Approximation and Projection (UMAP) preserves non-linear topological relationships and local dataset structures better than PCA, exposing subtle batch effects or complex patient groupings.
set.seed(42)
umap_result <- umap(t(logTMM))
umap_data <- data.frame(
  SampleID = rownames(umap_result$layout),
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2]
)
umap_data <- merge(umap_data, targets, by = "SampleID")

umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Group, shape = Timepoint)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = PatientID), show.legend = FALSE, size = 3.5) + 
  theme_minimal() +
  labs(x = "UMAP 1 (Non-linear representation)", y = "UMAP 2 (Non-linear representation)") +
  scale_color_manual(values = c("DCI" = "#E31A1C", "nonDCI" = "#1F78B4"))

print(umap_plot)

# Export to various formats
ggsave(filename = file.path(output_dir, "UMAP_Kontrola_Jakosci.pdf"), plot = umap_plot, width = 8, height = 6)
ggsave(filename = file.path(output_dir, "UMAP_Kontrola_Jakosci.png"), plot = umap_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = file.path(output_dir, "UMAP_Kontrola_Jakosci.svg"), plot = umap_plot, width = 8, height = 6)