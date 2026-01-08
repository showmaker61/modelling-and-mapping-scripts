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
  MAE = numeric(),       
  Bias = numeric(),     
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
  obs <- test_set$Y_target
  

  residuals <- preds - obs
  
  r2   <- cor(preds, obs)^2
  rmse <- sqrt(mean(residuals^2))
  mae  <- mean(abs(residuals))   
  bias <- mean(residuals)        
  

  validation_results[nrow(validation_results)+1, ] <- c(
    as.character(comp), nrow(test_set), r2, rmse, mae, bias
  )
}


cols_to_numeric <- c("Sample_Size", "R2", "RMSE", "MAE", "Bias")
validation_results[cols_to_numeric] <- lapply(validation_results[cols_to_numeric], as.numeric)


print(validation_results)

avg_stats <- colMeans(validation_results[, c("R2", "RMSE", "MAE", "Bias")])
print(round(avg_stats, 4))