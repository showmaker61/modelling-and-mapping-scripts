library(ggplot2)
library(dplyr)
library(ranger)

rm(list = ls()[!ls() %in% c("x_tr_final", "y_full", "best_row")]) 
gc()

full_data_loco <- data.frame(x_tr_final, Y_target = y_full)


full_data_loco$Compound_ID <- as.factor(round(full_data_loco$logKow, 2))
unique_compounds <- unique(full_data_loco$Compound_ID)
validation_results <- data.frame(
  Compound = character(),
  Sample_Size = numeric(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

for (comp in unique_compounds) {
  

  
  idx_test  <- which(full_data_loco$Compound_ID == comp)
  idx_train <- which(full_data_loco$Compound_ID != comp)
  
  train_set <- full_data_loco[idx_train, ]
  test_set  <- full_data_loco[idx_test, ]
  

  model <- ranger(
    formula = Y_target ~ ., 
    data = train_set %>% select(-Compound_ID),
    num.trees = if(exists("best_row")) best_row$num.trees else 500,
    mtry = if(exists("best_row")) best_row$mtry else 3,
    min.node.size = if(exists("best_row")) best_row$min.node.size else 5,
    splitrule = "extratrees",
    verbose = FALSE
  )
  

  preds <- predict(model, data = test_set %>% select(-Compound_ID, -Y_target))$predictions
  

  r2 <- cor(preds, test_set$Y_target)^2
  rmse <- sqrt(mean((preds - test_set$Y_target)^2))
  

  validation_results[nrow(validation_results)+1, ] <- c(
    as.character(comp), nrow(test_set), r2, rmse
  )
  
  cat(sprintf("  >>  %s (N=%d, R2=%.3f)\n", comp, nrow(test_set), r2))
}
