rm(list = ls()) 
gc()

packages <- c("dplyr", "caret", "gbm", "RANN", "missForest")

for(p in packages){
  if(!require(p, character.only = T)) install.packages(p)
  library(p, character.only = T)
}

set.seed(123) 
soil_data <- read.csv("C:/Users/电击小子/Desktop/SoilID.csv", header = TRUE)

colnames(soil_data)[1] <- "Soil_ID"
soil_data$Soil_ID <- as.character(soil_data$Soil_ID)

if("T" %in% colnames(soil_data)) {
  soil_data <- soil_data %>% rename(Temperature = T)
}

target_col <- "ln_Qe"
feature_cols <- c("OC", "pH", "Clay", "CEC", "logKow", "lnCe", "Temperature", "Ratio")

clean_df <- soil_data %>% 
  select(Soil_ID, all_of(target_col), all_of(feature_cols)) %>%
  filter(!is.na(.[[target_col]]))

for(col in feature_cols) {
  clean_df[[col]] <- as.numeric(clean_df[[col]])
}

unique_soils <- unique(clean_df$Soil_ID)

split_ratio <- 0.9 
train_ids <- sample(unique_soils, size = floor(split_ratio * length(unique_soils)))

train_data <- clean_df %>% filter(Soil_ID %in% train_ids)
test_data  <- clean_df %>% filter(!Soil_ID %in% train_ids)

hyper_grid <- expand.grid(
  interaction.depth = c(1, 3, 5, 7),
  n.trees = c(500, 1000, 1500, 2000),
  shrinkage = c(0.01, 0.05, 0.1),
  n.minobsinnode = c(5, 10, 15, 20),
  bag.fraction = c(0.5, 0.65, 0.8, 0.95)
)

grid_results <- data.frame(hyper_grid, Val_RMSE = NA, Val_R2 = NA, Train_R2 = NA)
folds_index <- groupKFold(train_data$Soil_ID, k = 9) 
n_folds <- length(folds_index)



for(j in 1:nrow(hyper_grid)) {
  
  d <- hyper_grid$interaction.depth[j]
  nt <- hyper_grid$n.trees[j]
  lr <- hyper_grid$shrinkage[j]
  min_obs <- hyper_grid$n.minobsinnode[j]
  bag <- hyper_grid$bag.fraction[j]
  
  fold_res <- data.frame(Val_RMSE=numeric(n_folds), Val_R2=numeric(n_folds), Train_R2=numeric(n_folds))
  
  for(i in 1:n_folds) {
    
    train_idx <- folds_index[[i]]
    f_train <- train_data[train_idx, ]
    f_val   <- train_data[-train_idx, ]
    
    x_tr <- f_train[, feature_cols]; y_tr <- f_train[[target_col]]
    x_val <- f_val[, feature_cols];  y_val <- f_val[[target_col]]
    

    nzv_pre <- preProcess(x_tr, method = "nzv")
    x_tr_clean <- predict(nzv_pre, x_tr)
    x_val_clean <- predict(nzv_pre, x_val)
    
    combined_x <- rbind(x_tr_clean, x_val_clean)
    n_tr_rows <- nrow(x_tr_clean)
    
    mf_obj <- missForest(combined_x, ntree = 20, verbose = FALSE) 
    combined_imp <- mf_obj$ximp
    
    x_tr_imp <- combined_imp[1:n_tr_rows, ]
    x_val_imp <- combined_imp[(n_tr_rows + 1):nrow(combined_imp), ]
    
    train_df_gbm <- data.frame(x_tr_imp, Y_target = y_tr)
    val_df_gbm   <- data.frame(x_val_imp) 
    
    m <- gbm(
      formula = Y_target ~ .,
      data = train_df_gbm,
      distribution = "gaussian",
      n.trees = nt,
      interaction.depth = d,
      shrinkage = lr,
      n.minobsinnode = min_obs,
      bag.fraction = bag,
      verbose = FALSE,
      n.cores = 1
    )
    
    p_tr <- predict(m, newdata = train_df_gbm, n.trees = nt)
    p_val <- predict(m, newdata = val_df_gbm, n.trees = nt)
    
    fold_res$Val_RMSE[i] <- sqrt(mean((p_val - y_val)^2))
    fold_res$Val_R2[i] <- cor(p_val, y_val)^2
    fold_res$Train_R2[i] <- cor(p_tr, y_tr)^2
  }
  
  grid_results$Val_RMSE[j] <- mean(fold_res$Val_RMSE)
  grid_results$Val_R2[j]   <- mean(fold_res$Val_R2)
  grid_results$Train_R2[j] <- mean(fold_res$Train_R2)
  
  if(j %% 1 == 0) {
    gap <- grid_results$Train_R2[j] - grid_results$Val_R2[j]
    cat(sprintf("组合 %2d (D:%d, Obs:%d) | Val: %.4f | Train: %.4f | Gap: %.3f\n", 
                j, d, min_obs, grid_results$Val_R2[j], grid_results$Train_R2[j], gap))
  }
}
best_idx <- which.min(grid_results$Val_RMSE)
best_row <- grid_results[best_idx, ]

x_full <- train_data[, feature_cols]; y_full <- train_data[[target_col]]
x_test <- test_data[, feature_cols];  y_test <- test_data[[target_col]]

final_nzv <- preProcess(x_full, method = "nzv")
x_full_clean <- predict(final_nzv, x_full)
x_test_clean <- predict(final_nzv, x_test)

n_full_rows <- nrow(x_full_clean)
combined_final <- rbind(x_full_clean, x_test_clean)

mf_final <- missForest(combined_final, ntree = 100, verbose = FALSE)
combined_final_imp <- mf_final$ximp

x_tr_final <- combined_final_imp[1:n_full_rows, ]
x_te_final <- combined_final_imp[(n_full_rows + 1):nrow(combined_final_imp), ]

train_final_df <- data.frame(x_tr_final, Y_target = y_full)
test_final_df  <- data.frame(x_te_final)

final_m_gbm <- gbm(
  formula = Y_target ~ .,
  data = train_final_df,
  distribution = "gaussian",
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction,
  verbose = FALSE
)

f_preds <- predict(final_m_gbm, newdata = test_final_df, n.trees = best_row$n.trees)
f_preds_train <- predict(final_m_gbm, newdata = train_final_df, n.trees = best_row$n.trees)
test_rmse <- sqrt(mean((f_preds - y_test)^2))
test_r2 <- cor(f_preds, y_test)^2
train_rmse <- sqrt(mean((f - y_full)^2))

cat("Training Set R² :", round(best_row$Train_R2, 4), "\n")
cat("Test Set R²     :", round(test_r2, 4), "\n")
cat("Test Set RMSE   :", round(test_rmse, 4), "\n")