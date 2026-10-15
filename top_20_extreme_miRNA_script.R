# TOP 20 HIGHEST AND LOWEST EXPRESSION SUMMARY

library(dplyr)

output_dir <- "Wyniki_Analizy_Top20"
dir.create(output_dir, showWarnings = FALSE)

# Step 1: Load integrated datasets
dane <- readRDS("00_dane_zintegrowane.rds")
logTMM <- dane$logTMM
targets <- dane$targets

samples_day2 <- targets$SampleID[targets$Timepoint == "Day2"]
samples_day6 <- targets$SampleID[targets$Timepoint == "Day6"]

# Step 2: Calculate means
expr_summary <- data.frame(
  miRNA = rownames(logTMM),
  Mean_logCPM_Global = rowMeans(logTMM),
  Mean_logCPM_Day2   = rowMeans(logTMM[, samples_day2]),
  Mean_logCPM_Day6   = rowMeans(logTMM[, samples_day6])
)

# Step 3: Extract extremes
top20_highest <- expr_summary %>%
  arrange(desc(Mean_logCPM_Global)) %>%
  slice_head(n = 20)

top20_lowest <- expr_summary %>%
  arrange(Mean_logCPM_Global) %>%
  slice_head(n = 20)

# Step 4: Format decimal places
format_us <- function(df) {
  df$Mean_logCPM_Global <- sprintf("%.2f", df$Mean_logCPM_Global)
  df$Mean_logCPM_Day2   <- sprintf("%.2f", df$Mean_logCPM_Day2)
  df$Mean_logCPM_Day6   <- sprintf("%.2f", df$Mean_logCPM_Day6)
  return(df)
}

high_fmt <- format_us(top20_highest)
low_fmt  <- format_us(top20_lowest)

colnames(high_fmt) <- c("miRNA_Highest", "Global_Highest", "Day2_Highest", "Day6_Highest")
colnames(low_fmt)  <- c("miRNA_Lowest", "Global_Lowest", "Day2_Lowest", "Day6_Lowest")

# Combine into 20 rows
summary_table <- cbind(
  Rank = 1:20,
  high_fmt,
  low_fmt
)

# Step 5: Save in CSV format
write.csv(summary_table, 
          file = file.path(output_dir, "Tabela_Top20_Highest_vs_Lowest_logCPM.csv"), 
          row.names = FALSE)