# CONSENSUS RANKING FOR REFERENCE GENE TRIPLETS (NORMIRAZOR OUTPUT)

library(readr)
library(dplyr)
library(writexl)
library(here)

# Step 1: Load triplet files (CL=3) and rename columns for table output
genorm <- read_tsv(here("genorm.txt"), show_col_types = FALSE) %>% 
  select(Kombinacja = Name, Stability_GeNorm = Stability, Rank_G = Ranking)

normfinder <- read_tsv(here("normfinder.txt"), show_col_types = FALSE) %>% 
  select(Kombinacja = Name, Stability_NormFinder = Stability, Rank_N = Ranking)

bestkeeper <- read_tsv(here("bestkeeper.txt"), show_col_types = FALSE) %>% 
  select(Kombinacja = Name, Stability_BestKeeper = Stability, Rank_B = Ranking)

# Step 2: Merge results, calculate global ranking, and extract Top 15
tabela_trojek <- genorm %>%
  inner_join(normfinder, by = "Kombinacja") %>%
  inner_join(bestkeeper, by = "Kombinacja") %>%
  # BIOINFORMATIC JUSTIFICATION: Summing the independent rankings from geNorm, NormFinder, and BestKeeper generates a robust consensus metric, reducing the bias inherent in relying on any single algorithm's stability calculation. Lower sum indicates better global stability.
  mutate(Suma_Rang = Rank_G + Rank_N + Rank_B) %>%
  arrange(Suma_Rang) %>%
  slice_head(n = 15) %>%
  select(Kombinacja, Stability_GeNorm, Stability_NormFinder, Stability_BestKeeper)

# Step 3: Console preview and export to Excel and CSV formats
print(tabela_trojek)

write_xlsx(tabela_trojek, here("Top15_Kombinacje3_Stabilnosc.xlsx"))
write_csv(tabela_trojek, here("Top15_Kombinacje3_Stabilnosc.csv"))

cat("\nConsensus ranking complete. Top 15 reference combinations exported successfully.\n")