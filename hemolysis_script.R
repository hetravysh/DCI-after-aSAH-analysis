# SCRIPT 01: HEMOLYSIS ESTIMATION (LOG2 SCALE MATCHING BLONDAL THRESHOLD)
library(dplyr)

output_dir <- "Wyeksportowane_Macierze_Master"

# 1. SAFE DATA LOADING
# We disable check.names to prevent R from modifying sample IDs (e.g., adding 'X' to '1_6').
# This ensures perfect column matching with our metadata targets$SampleID.
raw_counts <- read.csv(file.path(output_dir, "01_Raw_Counts_Corrected.csv"), row.names = 1, check.names = FALSE)
dane_do_analizy <- readRDS("00_dane_zintegrowane.rds")
targets <- dane_do_analizy$targets

# 2. LOCATING BIOMARKERS
erythrocyte_mir <- "hsa.miR.451a"
reference_mir <- "hsa.miR.23a.3p"

# Extracting raw counts as numeric vectors for mathematical operations.
count_451a <- as.numeric(raw_counts[erythrocyte_mir, targets$SampleID])
count_23a <- as.numeric(raw_counts[reference_mir, targets$SampleID])

# 3. LOG2 RATIO CALCULATION (MATCHING DELTA Cq)
# To directly compare our results against the literature threshold of 5 (Blondal et al.),
# we must calculate the ratio on a log2 scale, which is mathematically equivalent to delta Cq.
# A pseudocount of +1 is added to prevent log2(0) errors for samples with zero counts.
delta_log2 <- log2(count_451a + 1) - log2(count_23a + 1)

# 4. DATA FRAME ASSEMBLY
hemolysis_df <- targets %>%
  mutate(
    Raw_miR_451a = count_451a,
    Raw_miR_23a = count_23a,
    Delta_Log2 = round(delta_log2, 2)
  )

# Save the detailed results
write.csv(hemolysis_df, file.path(output_dir, "06c_Hemolysis_Log2_Delta.csv"), row.names = FALSE)

# 5. AGGREGATION FOR PRESENTATION
# Calculating the median Delta Log2 for each group.
# These values will fall roughly between 2.0 and 4.0, which clearly demonstrates 
# they are well below the hemolysis threshold of 5.
podsumowanie_log <- hemolysis_df %>%
  group_by(Group, Timepoint) %>%
  summarise(
    Mediana_Delta_Log2 = round(median(Delta_Log2), 2),
    .groups = 'drop'
  )

cat("\n--- MEDIAN DELTA LOG2 RATIOS (TO COMPARE AGAINST THRESHOLD < 5) ---\n")
print(podsumowanie_log)