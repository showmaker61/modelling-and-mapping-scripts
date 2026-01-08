library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)


train_df_lm <- as.data.frame(x_tr_final)
train_df_lm$ln_Qe <- y_full

test_df_lm <- as.data.frame(x_te_final)
test_df_lm$ln_Qe <- y_test

model_freundlich <- lm(ln_Qe ~ lnCe, data = train_df_lm)

model_simple_mlr <- lm(ln_Qe ~ lnCe + OC, data = train_df_lm)

model_lfer <- lm(ln_Qe ~ lnCe + OC + logKow, data = train_df_lm)

pred_freundlich <- predict(model_freundlich, newdata = test_df_lm)
pred_simple     <- predict(model_simple_mlr, newdata = test_df_lm)
pred_lfer       <- predict(model_lfer, newdata = test_df_lm)
pred_et         <- f_preds 

results_df <- data.frame(
  Observed = y_test,
  ET_ML    = pred_et,
  LFER     = pred_lfer,
  Simple   = pred_simple,
  Freundlich = pred_freundlich,
  logKow   = test_df_lm$logKow
)

calc_metrics <- function(obs, pred, name) {
  r2 <- cor(obs, pred)^2
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  return(c(Model = name, R2 = r2, RMSE = rmse, MAE = mae))
}

perf_table <- rbind(
  calc_metrics(y_test, pred_et, "Machine Learning (ET)"),
  calc_metrics(y_test, pred_lfer, "LFER (Koc-based)"),
  calc_metrics(y_test, pred_simple, "Simple MLR (OC+Ce)"),
  calc_metrics(y_test, pred_freundlich, "Freundlich Isotherm")
)
perf_table <- as.data.frame(perf_table)
perf_table$R2 <- as.numeric(perf_table$R2)
perf_table$RMSE <- as.numeric(perf_table$RMSE)

print(perf_table)

rmse_base <- perf_table$RMSE[perf_table$Model == "LFER (Koc-based)"]
rmse_ml <- perf_table$RMSE[perf_table$Model == "Machine Learning (ET)"]
improvement <- (rmse_base - rmse_ml) / rmse_base * 100


library(dplyr)

train_df_lm <- as.data.frame(x_tr_final)
train_df_lm$ln_Qe <- y_full
test_df_lm <- as.data.frame(x_te_final)

# (A) Freundlich 
mod_freundlich <- lm(ln_Qe ~ lnCe, data = train_df_lm)
# (B) Simple MLR (OC + lnCe)
mod_simple <- lm(ln_Qe ~ lnCe + OC, data = train_df_lm)
# (C) Koc (OC + lnCe + logKow)
mod_lfer <- lm(ln_Qe ~ lnCe + OC + logKow, data = train_df_lm)

pred_freundlich <- predict(mod_freundlich, newdata = test_df_lm)
pred_simple     <- predict(mod_simple, newdata = test_df_lm)
pred_lfer       <- predict(mod_lfer, newdata = test_df_lm)
pred_et         <- f_preds

get_stats <- function(observed, predicted, model_name) {
  res <- predicted - observed
  r2 <- cor(predicted, observed)^2
  rmse <- sqrt(mean(res^2))
  mae <- mean(abs(res))
  bias <- mean(res)
  
  return(data.frame(
    Model = model_name,
    R2    = round(r2, 4),
    RMSE  = round(rmse, 4),
    MAE   = round(mae, 4),
    Bias  = round(bias, 4)
  ))
}

final_table <- rbind(
  get_stats(y_test, pred_freundlich, "Freundlich Model"),
  get_stats(y_test, pred_simple,     "Simple MLR (OC+Ce)"),
  get_stats(y_test, pred_lfer,       "LFER (Koc Model)"),
  get_stats(y_test, pred_et,         "Machine Learning (ET)")
)

print(final_table)