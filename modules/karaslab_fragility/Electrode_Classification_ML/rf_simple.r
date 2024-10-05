source("read_data.r")
library(glue)
library(h2o)
#h2o.shutdown()
h2o.init(max_mem_size = "10g")

all_column_names <- colnames(dt_final)
time_columns <- grep("^X", all_column_names, value = TRUE)

## Candidate: features_freq, features
my_features <- time_columns
## leave one out
predictions <- list()

# patient <- unique(dt_weight$patient_code)[1]
for (patient in unique(dt_final$patient_code)) {
    message(patient)
    dt_train <- dt_final|>
        dplyr::filter(patient_code!=patient)
    dt_test <- dt_final|>
        dplyr::filter(patient_code==patient)

    dt_train_h2o <- as.h2o(dt_train)
    dt_test_h2o <- as.h2o(dt_test)
    rf_model <- h2o.randomForest(x = my_features, y = "resect",
                            training_frame = dt_train_h2o, nfolds = 0,
                            seed = 1234)

    rf_probabilities <- predict(rf_model, newdata = dt_test_h2o, type = "prob")
    preds <- as.data.frame(rf_probabilities)$p1
    dt_test[["preds"]] <- preds
    predictions[[patient]] <- dt_test|>
        select(patient_code, condition, electrode, resect, preds)
}

predictions <- do.call(rbind, predictions)

## use average to aggregate the predictions
predictions_voted <- predictions|>
group_by(patient_code, electrode)|>
summarise(preds = mean(preds),
          resect = unique(resect))


## Performance metrics
library(pROC)
roc_score=roc(predictions_voted$resect, predictions_voted$preds)
coordinates <- coords(roc_score, x = "all", input = "threshold", ret = "all")
thresholds <- coordinates[,"threshold"]
ppv <- coordinates[,"ppv"]
npv <- coordinates[,"npv"]
tp <- coordinates[,"tp"]
fp <- coordinates[,"fp"]


## Positive Predictive Value vs True Positive
# ggplot(coordinates, aes(x = tp, y = ppv)) +
#   geom_point() +
#   geom_line() +
#   ylim(0,1)+
#   scale_x_continuous(breaks = seq(0,100,by=2), limits  = c(0,100))


# Identify the threshold with the highest TP where PPV is above 95%
valid_thresholds <- coordinates[ppv > 0.95,]
best_threshold <- valid_thresholds[which.max(valid_thresholds$tp),]

# Add a column for labels that only appear at every 0.1 threshold
coordinates <- coordinates %>%
  mutate(across(c(ppv, tp, fp), ~ifelse(!is.finite(.), NA, .))) %>%
  mutate(label_rounded = round(threshold, 1)) %>%
  group_by(label_rounded) %>%
  mutate(label = ifelse(threshold == min(threshold) & !is.na(ppv) & threshold >= 0.05, sprintf("%.1f", label_rounded), NA)) %>%
  ungroup()

coordinates <- coordinates %>%
  dplyr::filter(if_all(c(tp, ppv, fp, threshold), is.finite))

# Plotting code
plot <- ggplot(coordinates, aes(x = tp, y = ppv)) +
  geom_point() +  # Plot all points
  geom_line() +  # Connect the points with a line
  geom_point(data = best_threshold, aes(x = tp, y = ppv), colour = "red", size = 3) +
  geom_text(aes(label = label), vjust = 1.5, hjust = 2, check_overlap = TRUE, color = "blue", size = 5, fontface = "bold") +
  geom_text(data = best_threshold, aes(x = tp + 20, y = ppv, label = sprintf("Threshold: %.2f\n  TP: %d, FP: %d\n  PPV: %.2f%%", threshold, tp, fp, ppv * 100)), vjust = 1.5, hjust = -0.4, color = "red", size = 8, fontface = "bold") +  # Annotate the best threshold with detailed information
  geom_vline(xintercept = best_threshold$tp, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = best_threshold$ppv, linetype = "dashed", color = "red", size = 1) +
  ylim(0, 1) +
  scale_x_continuous(name = "True Positives", breaks = seq(min(coordinates$tp), max(coordinates$tp), by = 50), limits = c(min(coordinates$tp), max(coordinates$tp))) +  # Customize X-axis for raw count of TPs
  scale_y_continuous(name = "PPV") +
  ggtitle("Positive Predictive Value (PPV) vs. True Positives") +
  theme_minimal(base_size = 25) +
  theme(plot.title = element_text(size = 25))

print(plot)


plot <- seq(0.1, 1, by = 0.0001)
tp_fp_data <- data.frame(threshold = numeric(), true_positive = numeric(), false_positive = numeric())

## Voting for each threshold
for (threshold in plot) {

  set.seed(123)

  # predictions: patient_code condition electrode resect  preds

  predictions$pred_threshold <- ifelse(predictions$preds > threshold, 1, 0)

  result <- data.frame(patient_code = character(),
                       electrode = numeric(),
                       prediction = numeric(),
                       outcome = numeric())

  for (patient in unique(predictions$patient_code)) {
    df_patient <- predictions[predictions$patient_code==patient,]

    conditions <- unique(df_patient$condition)
    temp_cond <- df_patient[df_patient$condition==conditions[1],]
    electrodes <- temp_cond$electrode
    resect <- temp_cond$resect
    voting <- rep(0, length(electrodes))

    for (condition in unique(df_patient$condition)) {
      df_condition <- df_patient[df_patient$condition==condition,]
      voting <- voting + df_condition$pred_threshold
    }

    voting <- voting / length(unique(df_patient$condition))
    voting <- ifelse(voting >= 0.5, 1, 0)
    # voting <- ifelse(voting >= 0.5, 1, 0)
    # Repeat single values to match the length of the vectors
    patient_repeated <- rep(patient, length(electrodes))

    # Combine into a data frame
    row_data <- data.frame(patient = patient_repeated,
                           electrodes = electrodes,
                           voting = voting,
                           resect = resect)

    result <- rbind(result, row_data)
  }

  colnames(result) <- c("patient_code", "electrode", "prediction", "outcome")
  result$prediction <- factor(result$prediction)

  # Calculate TP, FP, TN, FN
  conf_matrix <- table(result$outcome, result$prediction)
  print(conf_matrix)

  if (ncol(conf_matrix) > 1) {
    # Extract TP, FP, TN, FN
    TP <- conf_matrix[2, 2]
    FP <- conf_matrix[1, 2]
    TN <- conf_matrix[1, 1]
    FN <- conf_matrix[2, 1]

    # Calculate metrics
    accuracy <- (TP + TN) / sum(conf_matrix)
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    ppv <- TP / (TP + FP)
    npv <- TN / (TN + FN)

    # Calculate AUC
    result$prediction_numeric <- as.numeric(result$prediction) - 1
    auc <- roc(result$outcome, result$prediction_numeric)$auc

    # Print metrics
    cat("Overall AUC:", auc, "\n")
    cat("Accuracy:", accuracy, "\n")
    cat("Sensitivity:", sensitivity, "\n")
    cat("Specificity:", specificity, "\n")
    cat("Positive Predictive Value (PPV):", ppv, "\n")
    cat("Negative Predictive Value (NPV):", npv, "\n")
    cat("Number of True Positives:",TP)
    cat("Number of False Positives:",FP)
    cat("Number of True Negatives:",TN)
    cat("Number of False Negatives:",FN)
  } else if (colnames(conf_matrix) == "1" ) {
    FP <- conf_matrix[1, 1]
    TP <- conf_matrix[2, 1]
  } else if (colnames(conf_matrix) == "0" ) {
    FP <- 0
    TP <- 0
  }

  tp_fp_data <- rbind(tp_fp_data, data.frame(threshold = threshold, true_positive = TP, false_positive = FP))

}

print(tp_fp_data)

tp_fp_data$PPV <- tp_fp_data$true_positive / (tp_fp_data$true_positive + tp_fp_data$false_positive)

# Identify the point with the most true positives and no false positives
optimal_point <- tp_fp_data[tp_fp_data$PPV > 0.95, ]
optimal_point <- optimal_point[which.max(optimal_point$true_positive),]

# Plotting
tp_fp_data$threshold <- round(tp_fp_data$threshold, 2)
# Define desired thresholds
desired_thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

# Filter data for these thresholds and select the point with the highest PPV for each threshold
label_data <- tp_fp_data %>%
  dplyr::filter(threshold %in% desired_thresholds) %>%
  group_by(threshold) %>%
  slice(which.max(PPV)) %>%
  ungroup()

# Plotting code with ensured unique labels
plot <- ggplot(tp_fp_data, aes(x = true_positive, y = PPV)) +
  geom_point() +
  geom_line() +
  geom_point(data = optimal_point, aes(x = true_positive, y = PPV), colour = "red", size = 3) +
  geom_text(data = label_data, aes(label = sprintf("%.1f", threshold)), vjust = 1.5, hjust = 2, check_overlap = TRUE, color = "blue", size = 5, fontface = "bold") +
  geom_text(data = optimal_point, aes(x = true_positive, y = PPV,
                                      label = sprintf("Threshold: %.2f%%\n   TP: %g, FP: %g\n   PPV: %.2f%%", threshold * 100, true_positive, false_positive, PPV * 100)),
            vjust = 1.5, hjust = -0.4, color = "red", size = 8, fontface = "bold") +
  geom_vline(xintercept = optimal_point$true_positive, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = optimal_point$PPV, linetype = "dashed", color = "red", size = 1) +
  ylim(0, 1) +
  scale_x_continuous(name = "True Positives") +
  scale_y_continuous(name = "PPV") +
  ggtitle("Positive Predictive Value (PPV) vs. True Positives") +
  theme_minimal(base_size = 25) +
  theme(plot.title = element_text(size = 25))

print(plot)

h2o.shutdown()





