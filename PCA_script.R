# DIMENSIONALITY REDUCTION: PCA AND LONGITUDINAL SHIFT ANALYSIS

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(svglite) 

# Step 0: Initialize output directory
output_dir <- "Wyniki_Analizy_PCA"
dir.create(output_dir, showWarnings = FALSE)

# Step 1: Load integrated datasets
dane <- readRDS("00_dane_zintegrowane.rds")
logTMM <- dane$logTMM
targets <- dane$targets

# Step 2: Principal Component Analysis
# STATISTICAL JUSTIFICATION (PCA): Principal Component Analysis provides linear orthogonal transformations. 
# Using scale.=TRUE standardizes features to zero mean and unit variance, preventing miRNAs with exceptionally high baseline expression from dominating the principal components.
pca <- prcomp(t(logTMM), scale. = TRUE)

# Prepare dataframe for plotting
pca_data <- data.frame(SampleID = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
pca_data <- merge(pca_data, targets, by = "SampleID")

# Generate main PCA scatter plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = SampleID), size = 3, show.legend = FALSE, max.overlaps = Inf) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
  ) +
  scale_color_manual(values = c("DCI" = "#E31A1C", "nonDCI" = "#1F78B4")) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank()
  )

print(pca_plot)
ggsave(filename = file.path(output_dir, "PCA_labels_no_title.png"), plot = pca_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = file.path(output_dir, "PCA_labels_no_title.svg"), plot = pca_plot, width = 8, height = 6)

# Step 3: Longitudinal shift analysis
pca_wide <- pca_data %>%
  select(PatientID, Group, Timepoint, PC1, PC2) %>%
  pivot_wider(names_from = Timepoint, values_from = c(PC1, PC2))

# Calculate PC1 differential and Euclidean distance
pca_distances <- pca_wide %>%
  mutate(
    PC1_diff = PC1_Day6 - PC1_Day2,
    PCA_distance = sqrt((PC1_Day6 - PC1_Day2)^2 + (PC2_Day6 - PC2_Day2)^2)
  )

pos <- position_jitter(width = 0.2, seed = 42)

# Boxplot A: Differential shift in PC1 axis only
pc1_diff_plot <- ggplot(pca_distances, aes(x = Group, y = PC1_diff, fill = Group, label = PatientID)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position = pos, size = 3, alpha = 0.8, shape = 21, color = "black") +
  geom_text_repel(position = pos, size = 3, show.legend = FALSE, max.overlaps = Inf) +
  theme_minimal() +
  scale_fill_manual(values = c("DCI" = "#E31A1C", "nonDCI" = "#1F78B4")) +
  labs(x = "Group", y = "Change in PC1 (Day 6 - Day 2)") +
  theme(legend.position = "none")

print(pc1_diff_plot)
ggsave(filename = file.path(output_dir, "PC1_diff_boxplot.png"), plot = pc1_diff_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = file.path(output_dir, "PC1_diff_boxplot.svg"), plot = pc1_diff_plot, width = 8, height = 6)

# Boxplot B: Overall Euclidean distance (PC1 + PC2)
pca_dist_plot <- ggplot(pca_distances, aes(x = Group, y = PCA_distance, fill = Group, label = PatientID)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(position = pos, size = 3, alpha = 0.8, shape = 21, color = "black") +
  geom_text_repel(position = pos, size = 3, show.legend = FALSE, max.overlaps = Inf) +
  theme_minimal() +
  scale_fill_manual(values = c("DCI" = "#E31A1C", "nonDCI" = "#1F78B4")) +
  labs(x = "Group", y = "PCA Euclidean Distance (Day 2 vs Day 6)") +
  theme(legend.position = "none")

print(pca_dist_plot)
ggsave(filename = file.path(output_dir, "PCA_distance_boxplot.png"), plot = pca_dist_plot, width = 8, height = 6, dpi = 300)
ggsave(filename = file.path(output_dir, "PCA_distance_boxplot.svg"), plot = pca_dist_plot, width = 8, height = 6)

# Step 4: Statistical Verification
# STATISTICAL JUSTIFICATION (WILCOXON): The Euclidean distance delta between Day 2 and Day 6 in PCA space is compared using the Wilcoxon rank-sum test. 
# This non-parametric test is robust to extreme outliers and does not assume a normal distribution of trajectory distances.
cat("\nWilcoxon test for shift in PC1 axis (DCI vs nonDCI)\n")
wilcox_pc1 <- wilcox.test(PC1_diff ~ Group, data = pca_distances)
print(wilcox_pc1)

cat("\nWilcoxon test for total PCA distance (DCI vs nonDCI)\n")
wilcox_dist <- wilcox.test(PCA_distance ~ Group, data = pca_distances)
print(wilcox_dist)

# Save distance results
write.table(pca_distances, file = file.path(output_dir, "PCA_Odleglosci_Pacjentow.csv"), sep = ";", dec = ",", row.names = FALSE)