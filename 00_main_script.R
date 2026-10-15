# SCRIPT 00: CENTRAL MASTER
# LIQUID BIOPSY FILTERING AND TMM NORMALIZATION

library(edgeR)

output_dir <- "Wyeksportowane_Macierze_Master"
dir.create(output_dir, showWarnings = FALSE)

# Load and transpose raw count matrix
# BIOINFORMATIC JUSTIFICATION: edgeR and limma-voom strictly require unnormalized, raw integer counts to accurately model the mean-variance relationship of sequencing data.
counts <- read.csv("mature_counts.csv", row.names=1)
counts <- as.data.frame(t(counts))
counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))

# Correct naming error for sample 1
sample_names <- colnames(counts)
sample_names[sample_names == "1_4"] <- "1_6"
colnames(counts) <- sample_names

# Construct metadata dataframe
patient_ids <- sapply(strsplit(sample_names, "_"), function(x) x[1])
timepoints <- sapply(strsplit(sample_names, "_"), function(x) x[2])
targets <- data.frame(SampleID = sample_names, PatientID = patient_ids, Timepoint = paste0("Day", timepoints))

dci_list <- c("1", "2", "7", "10", "14", "15", "18", "38", "43", "48")
targets$Group <- ifelse(targets$PatientID %in% dci_list, "DCI", "nonDCI")
targets$Group <- factor(targets$Group, levels = c("nonDCI", "DCI"))
targets$Timepoint <- factor(targets$Timepoint, levels = c("Day2", "Day6"))
rownames(targets) <- targets$SampleID
targets$Group_Time <- factor(paste(targets$Group, targets$Timepoint, sep="_"))

# Export raw counts and metadata before any filtering
write.csv(counts, file.path(output_dir, "01_Raw_Counts_Corrected.csv"), row.names = TRUE)
write.csv(targets, file.path(output_dir, "02_Metadata_Targets.csv"), row.names = FALSE)

# Initialize edgeR object
y <- DGEList(counts = counts, samples = targets)

# STRICT FILTERING FOR LIQUID BIOPSY
# STATISTICAL JUSTIFICATION: Liquid biopsies inherently suffer from low RNA yield and high background noise. 
# Filtering low-expressed miRNAs (min.count = 50) prevents the mean-variance relationship from being skewed by noise and reduces the multiple testing burden (FDR), increasing statistical power.
keep_global <- filterByExpr(y, group = targets$Group_Time, min.count = 50)
y_filtered <- y[keep_global, , keep.lib.sizes=FALSE]

cat("\nGLOBAL STATISTICS (LIQUID BIOPSY FILTER)\n")
cat("Rejected (FALSE) and retained (TRUE) miRNAs in the dataset:\n")
print(table(keep_global))
cat("\n")

# Export filtered raw counts
write.csv(y_filtered$counts, file.path(output_dir, "03_Counts_Filtered.csv"), row.names = TRUE)

# TMM NORMALIZATION AND TRANSFORMATIONS
# BIOINFORMATIC JUSTIFICATION: Trimmed Mean of M-values (TMM) accounts for RNA compositional bias across samples, which is critical when a few highly expressed miRNAs could artificially suppress the apparent expression of others.
y_filtered <- calcNormFactors(y_filtered, method = "TMM")

# Generate log2-CPM for downstream visualization (PCA, UMAP, Heatmaps)
logTMM <- cpm(y_filtered, log = TRUE) 
# Generate non-log CPM for non-parametric rank correlation tests (Spearman)
cpm_TMM <- cpm(y_filtered, log = FALSE) 

# Export normalized matrices
write.csv(logTMM, file.path(output_dir, "04_logTMM_Normalized.csv"), row.names = TRUE)
write.csv(cpm_TMM, file.path(output_dir, "05_CPM_TMM_Normalized.csv"), row.names = TRUE)

# Save to a single serialized R object for data integrity across downstream scripts
dane_do_analizy <- list(
  counts_filtered = y_filtered$counts,
  logTMM = logTMM,
  cpm_TMM = cpm_TMM,
  targets = targets
)
saveRDS(dane_do_analizy, "00_dane_zintegrowane.rds")
cat("\nSuccessfully generated 00_dane_zintegrowane.rds file. Ready for downstream analysis.\n")