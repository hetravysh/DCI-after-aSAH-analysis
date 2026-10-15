# MACHINE LEARNING LOOCV FOLD DETAILS AND ELASTIC NET

library(dplyr)
library(caret) 
library(pROC) 
library(ggplot2)
library(glmnet)

output_dir <- "Wyniki_Analizy_ML_Folds"
dir.create(output_dir, showWarnings = FALSE)

# Step 1: Load master data
dane <- readRDS("00_dane_zintegrowane.rds")
logTMM <- dane$logTMM 
targets <- dane$targets

clinical_data <- read.csv("roboczy_kliniczny.csv", sep = ";")
clinical_data$WFNS <- as.numeric(gsub(",", ".", as.character(clinical_data[,6])))
clinical_data$mFisher <- as.numeric(gsub(",", ".", as.character(clinical_data[,5])))

targets_d2 <- subset(targets, Timepoint == "Day2")
valid_patient_ids <- as.numeric(targets_d2$PatientID)

clinical_d2 <- clinical_data[clinical_data$Lp %in% valid_patient_ids, ]
clinical_d2 <- clinical_d2[match(valid_patient_ids, clinical_d2$Lp), ]

# Step 2: Elastic Net Selection
X_mirna <- t(logTMM[, targets_d2$SampleID]) 
X_clin <- as.matrix(clinical_d2[, c("WFNS", "mFisher")])
colnames(X_clin) <- c("WFNS", "mFisher")

X_matrix <- cbind(X_clin, X_mirna)
Y_vector <- factor(targets_d2$Group, levels = c("nonDCI", "DCI"))

# STATISTICAL JUSTIFICATION: Penalty factors allow the algorithm to shrink redundant miRNA coefficients towards zero while forcing established clinical variables (WFNS, mFisher) into the model unpenalized (penalty = 0).
# Elastic Net (alpha = 0.5) is chosen over pure Lasso to handle high collinearity typical in miRNA co-expression clusters.
p_fac <- c(0, 0, rep(1, ncol(X_mirna)))

set.seed(123)
cv_model <- cv.glmnet(X_matrix, Y_vector, family = "binomial", alpha = 0.5, 
                      nfolds = nrow(X_matrix), type.measure = "deviance", 
                      penalty.factor = p_fac, grouped = FALSE)

coefs <- coef(cv_model, s = "lambda.min")
non_zero <- coefs[which(coefs != 0), ]
selected_mirna_codes <- names(non_zero)[!names(non_zero) %in% c("(Intercept)", "WFNS", "mFisher")]

mirna_coefficients <- non_zero[selected_mirna_codes]
best_mirna <- names(which.max(abs(mirna_coefficients)))

# Step 3: Construct matrix and filter patient 43 for training
expr_selected <- logTMM[best_mirna, targets_d2$SampleID]

df_model_full <- data.frame(
  PatientID = targets_d2$PatientID,
  Outcome = factor(targets_d2$Group, levels = c("nonDCI", "DCI")), 
  WFNS = clinical_d2$WFNS,
  mFisher = clinical_d2$mFisher,
  moje_miRNA = expr_selected
)

# STATISTICAL JUSTIFICATION: Outlier removal. Patient 43 exhibited atypical, overly rapid DCI onset, which skews the predictive clinical timeline and introduces confounding variance to the model.
excluded_patients <- c("43") 
df_model <- df_model_full[!df_model_full$PatientID %in% excluded_patients, ]

# Retain row indices for LOOCV fold identification
patient_indices <- df_model$PatientID
df_model <- df_model %>% select(-PatientID)

# Step 4: Configure LOOCV parameters
# STATISTICAL JUSTIFICATION: Leave-One-Out Cross-Validation (LOOCV) is selected due to the limited sample size (n=19). 
# It provides an nearly unbiased estimate of generalization error by training on n-1 samples and testing on the single held-out sample, maximizing the training data efficiency.
train_control <- trainControl(
  method = "LOOCV",             
  classProbs = TRUE,            
  summaryFunction = twoClassSummary, 
  savePredictions = "all"
)

# Step 5: Train validation models
set.seed(123)
m_mir   <- train(Outcome ~ moje_miRNA, data = df_model, method="glm", family="binomial", trControl=train_control, metric="ROC")
m_wfns  <- train(Outcome ~ WFNS, data = df_model, method="glm", family="binomial", trControl=train_control, metric="ROC")
m_mfish <- train(Outcome ~ mFisher, data = df_model, method="glm", family="binomial", trControl=train_control, metric="ROC")
m_comb1 <- train(Outcome ~ WFNS + moje_miRNA, data = df_model, method="glm", family="binomial", trControl=train_control, metric="ROC")
m_comb2 <- train(Outcome ~ mFisher + moje_miRNA, data = df_model, method="glm", family="binomial", trControl=train_control, metric="ROC")

clean_name <- gsub("\\.", "-", best_mirna)

# Step 6: Aggregate results and fold details
# BIOINFORMATIC JUSTIFICATION: Evaluating both Sensitivity and Specificity alongside AUC prevents misleading interpretations in imbalanced datasets. Absolute counts of TP/TN/FP/FN provide clinical context to the errors.
get_full_statistics <- function(caret_model, model_name) {
  pred_t <- caret_model$pred
  tp <- sum(pred_t$obs == "DCI" & pred_t$pred == "DCI")
  tn <- sum(pred_t$obs == "nonDCI" & pred_t$pred == "nonDCI")
  fp <- sum(pred_t$obs == "nonDCI" & pred_t$pred == "DCI")
  fn <- sum(pred_t$obs == "DCI" & pred_t$pred == "nonDCI")
  
  data.frame(Model = model_name, TP = tp, TN = tn, FP = fp, FN = fn)
}

full_table <- bind_rows(
  get_full_statistics(m_mir, paste("Baseline:", clean_name)),
  get_full_statistics(m_wfns, "Baseline: WFNS Score"),
  get_full_statistics(m_mfish, "Baseline: mFisher Score"),
  get_full_statistics(m_comb1, paste("Combined: WFNS +", clean_name)),
  get_full_statistics(m_comb2, paste("Combined: mFisher +", clean_name))
)

model_results <- data.frame(
  Model = full_table$Model,
  AUC = sprintf("%.3f", c(m_mir$results$ROC, m_wfns$results$ROC, m_mfish$results$ROC, m_comb1$results$ROC, m_comb2$results$ROC)),
  Sensitivity = sprintf("%.3f", c(m_mir$results$Sens, m_wfns$results$Sens, m_mfish$results$Sens, m_comb1$results$Sens, m_comb2$results$Sens)),
  Specificity = sprintf("%.3f", c(m_mir$results$Spec, m_wfns$results$Spec, m_mfish$results$Spec, m_comb1$results$Spec, m_comb2$results$Spec)),
  TP = full_table$TP, TN = full_table$TN, FP = full_table$FP, FN = full_table$FN
)
write.csv(model_results, file.path(output_dir, paste0("Ostateczna_Zsumowana_Macierz_", clean_name, ".csv")), row.names = FALSE, quote = FALSE)

# Detailed results per specific FOLD (Iteration)
fold_details <- bind_rows(
  m_mir$pred %>% mutate(Model = paste("Baseline:", clean_name)),
  m_wfns$pred %>% mutate(Model = "Baseline: WFNS Score"),
  m_mfish$pred %>% mutate(Model = "Baseline: mFisher Score"),
  m_comb1$pred %>% mutate(Model = paste("Combined: WFNS +", clean_name)),
  m_comb2$pred %>% mutate(Model = paste("Combined: mFisher +", clean_name))
) %>%
  mutate(
    Classification_Type = case_when(
      obs == "DCI" & pred == "DCI" ~ "TP",
      obs == "nonDCI" & pred == "nonDCI" ~ "TN",
      obs == "nonDCI" & pred == "DCI" ~ "FP",
      obs == "DCI" & pred == "nonDCI" ~ "FN"
    ),
    Test_Patient_ID = patient_indices[rowIndex]
  ) %>%
  select(Model, Fold_Number = rowIndex, Test_Patient_ID, Actual_Class = obs, Predicted_Class = pred, DCI_Probability = DCI, Classification_Type) %>%
  arrange(Model, Test_Patient_ID)

write.csv(fold_details, file.path(output_dir, paste0("Szczegoly_Wszystkich_Foldow_", clean_name, ".csv")), row.names = FALSE, quote = FALSE)

# Step 7: ROC Curves composed of all LOOCV folds combined
roc_mir   <- roc(m_mir$pred$obs, m_mir$pred$DCI, quiet=TRUE)
roc_wfns  <- roc(m_wfns$pred$obs, m_wfns$pred$DCI, quiet=TRUE)
roc_mfish <- roc(m_mfish$pred$obs, m_mfish$pred$DCI, quiet=TRUE)
roc_comb1 <- roc(m_comb1$pred$obs, m_comb1$pred$DCI, quiet=TRUE)
roc_comb2 <- roc(m_comb2$pred$obs, m_comb2$pred$DCI, quiet=TRUE)

plot_single_roc <- function(roc_obj, plot_title, line_color, file_name) {
  auc_val <- sprintf("%.3f", auc(roc_obj))
  p <- ggroc(roc_obj, color = line_color, linewidth = 1.2) +
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1, color="black", linetype="dashed", linewidth=0.5) +
    theme_minimal(base_size = 14) +
    labs(title = plot_title, subtitle = paste("AUC (Aggregated LOOCV) =", auc_val),
         x = "Specificity", y = "Sensitivity")
  
  ggsave(file.path(output_dir, paste0(file_name, ".png")), plot = p, width = 6, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(file_name, ".svg")), plot = p, width = 6, height = 6, device = "svg")
}

plot_single_roc(roc_mir, paste("Baseline:", clean_name), "#1F78B4", paste0("ROC_1_", clean_name))
plot_single_roc(roc_wfns, "Baseline: WFNS Score", "grey40", "ROC_2_WFNS")
plot_single_roc(roc_mfish, "Baseline: mFisher Score", "grey40", "ROC_3_mFisher")
plot_single_roc(roc_comb1, paste("Combined: WFNS +", clean_name), "#E31A1C", paste0("ROC_4_WFNS_", clean_name))
plot_single_roc(roc_comb2, paste("Combined: mFisher +", clean_name), "#E31A1C", paste0("ROC_5_mFisher_", clean_name))