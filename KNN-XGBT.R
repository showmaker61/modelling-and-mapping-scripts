rm(list = ls()) 
gc()

packages <- c("dplyr", "caret", "xgboost", "RANN")

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
train_ids <- sample(unique_soils, size = floor(0.9 * length(unique_soils)))

train_data <- clean_df %>% filter(Soil_ID %in% train_ids)
test_data  <- clean_df %>% filter(!Soil_ID %in% train_ids)


hyper_grid <- expand.grid(
       max_depth = c(3, 5, 7, 9),        
       min_child_weight = c(1, 5, 10, 12, 15, 18),     
       eta = c(0.01, 0.05, 0.1, 0.3),              
       gamma = c(0, 0.1, 0.5, 1, 2),
       alpha = c(0, 0.1, 1, 3, 5),
       colsample_bytree = c(0.5, 0.6, 0.7, 0.8, 0.9),      
       subsample = c(0.5, 0.6, 0.7, 0.8, 0.9),
       nrounds = c(300, 500, 1000, 1500, 1800)
)


grid_results <- data.frame(hyper_grid, Val_RMSE = NA, Val_R2 = NA, Train_R2 = NA)

folds_index <- groupKFold(train_data$Soil_ID, k = 9)
n_folds <- length(folds_index)

for(j in 1:nrow(hyper_grid)) {
  current_nrounds <- hyper_grid$nrounds[j]
  params <- list(
    booster = "gbtree", objective = "reg:squarederror",
    eta = hyper_grid$eta[j], 
    max_depth = hyper_grid$max_depth[j],
    subsample = hyper_grid$subsample[j], 
    colsample_bytree = hyper_grid$colsample_bytree[j],
    min_child_weight = hyper_grid$min_child_weight[j], 
    gamma = hyper_grid$gamma[j],
    alpha = hyper_grid$alpha[j]
  )
  
  fold_rmse <- numeric(n_folds)
  fold_r2   <- numeric(n_folds)
  fold_tr_r2 <- numeric(n_folds)
  
  for(i in 1:n_folds) {
    
    train_idx <- folds_index[[i]]
    f_train <- train_data[train_idx, ]
    f_val   <- train_data[-train_idx, ]
    
    x_tr <- f_train[, feature_cols]; y_tr <- f_train[[target_col]]
    x_val <- f_val[, feature_cols];  y_val <- f_val[[target_col]]

    imputer <- preProcess(x_tr, method = c("nzv", "knnImpute"), k = 5) 
    
    x_tr_imp <- predict(imputer, x_tr)
    x_val_imp <- predict(imputer, x_val)
    
    dtr <- xgb.DMatrix(as.matrix(x_tr_imp), label = y_tr)
    dval <- xgb.DMatrix(as.matrix(x_val_imp), label = y_val)

    
    m <- xgb.train(
      params = params, 
      data = dtr, 
      nrounds = current_nrounds,
      verbose = 0
    )
    
    
    p_tr <- predict(m, dtr)
    p_val <- predict(m, dval)
    
    fold_rmse[i] <- sqrt(mean((p_val - y_val)^2))
    fold_r2[i]   <- cor(p_val, y_val)^2
    fold_tr_r2[i] <- cor(p_tr, y_tr)^2
  }
  
  grid_results$Val_RMSE[j] <- mean(fold_rmse)
  grid_results$Val_R2[j]   <- mean(fold_r2)
  grid_results$Train_R2[j] <- mean(fold_tr_r2)
  
  if(j %% 1 == 0) {
    gap <- grid_results$Train_R2[j] - grid_results$Val_R2[j]
    cat(sprintf("组合 %2d (D:%d, CW:%d) | Val: %.4f | Train: %.4f | Gap: %.3f\n", 
                j, hyper_grid$max_depth[j], hyper_grid$min_child_weight[j],
                grid_results$Val_R2[j], grid_results$Train_R2[j], gap))
  }
}


best_idx <- which.min(grid_results$Val_RMSE)
best_row <- grid_results[best_idx, ]


x_full <- train_data[, feature_cols]; y_full <- train_data[[target_col]]
x_test <- test_data[, feature_cols];  y_test <- test_data[[target_col]]

final_imp <- preProcess(x_full, method = c("nzv", "knnImpute"), k = 5)
x_tr_final <- predict(final_imp, x_full)
x_te_final <- predict(final_imp, x_test)

f_params <- list(
  booster="gbtree", objective="reg:squarederror",
  eta=best_row$eta, max_depth=best_row$max_depth,
  subsample=best_row$subsample, colsample_bytree=best_row$colsample_bytree,
  min_child_weight=best_row$min_child_weight, gamma=best_row$gamma,
  alpha=best_row$alpha
)

final_nrounds <- best_row$nrounds


final_m <- xgb.train(
  params = f_params, 
  data = xgb.DMatrix(as.matrix(x_tr_final), label = y_full), 
  nrounds = final_nrounds, 
  verbose = 0
)

f_preds <- predict(final_m, xgb.DMatrix(as.matrix(x_te_final)))

test_rmse <- sqrt(mean((f_preds - y_test)^2))
test_r2 <- cor(f_preds, y_test)^2

cat("Training Set R² :", round(best_row$Train_R2, 4), "\n")
cat("Test Set R²     :", round(test_r2, 4), "\n")
cat("Test Set RMSE   :", round(test_rmse, 4), "\n")