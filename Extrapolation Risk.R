rm(list = ls()[!ls() %in% c("modeldata", "china_sample", "china_vect")]) 
gc()

packages <- c("ggplot2", "dplyr", "viridis", "sf", "terra") 
for(p in packages) library(p, character.only = T)

train_df <- modeldata[, 1:4]
colnames(train_df) <- c("OC", "pH", "CEC", "Clay")

target_cols <- c("OC", "pH", "CEC", "Clay")

if(all(target_cols %in% colnames(china_sample))) {
  pred_df <- china_sample[, target_cols]
} else {
  pred_df <- china_sample[, 3:6]
  colnames(pred_df) <- target_cols
}

coords <- china_sample[, 1:2]
colnames(coords) <- c("x", "y")
calculate_mess <- function(train, pred) {
  
  vars <- colnames(train)
  n_pred <- nrow(pred)
  n_vars <- length(vars)
  
  sim_scores <- matrix(NA, nrow = n_pred, ncol = n_vars)
  colnames(sim_scores) <- vars
  
  for(v in vars) {
    tr_vals <- train[[v]]
    min_tr <- min(tr_vals, na.rm=T)
    max_tr <- max(tr_vals, na.rm=T)
    range_tr <- max_tr - min_tr
    
    p_vals <- pred[[v]]
    
    score_min <- (p_vals - min_tr) / range_tr * 100
    score_max <- (max_tr - p_vals) / range_tr * 100
    
    sim_scores[, v] <- pmin(score_min, score_max)
  }
  
  mess_final <- apply(sim_scores, 1, min)
  most_limiting <- colnames(sim_scores)[max.col(-sim_scores)]
  
  return(list(MESS = mess_final, Limiting_Factor = most_limiting))
}

mess_res <- calculate_mess(train_df, pred_df)

plot_df <- data.frame(
  x = coords$x,
  y = coords$y,
  MESS = mess_res$MESS,
  Factor = mess_res$Limiting_Factor
)

plot_df$Status <- ifelse(plot_df$MESS < 0, "Extrapolation (Negative)", "Reliable (Positive)")


pct_extrap <- mean(plot_df$MESS < 0) * 100

print(summary(plot_df$MESS)) 
p_mess <- ggplot() +
  geom_point(data = plot_df, aes(x = x, y = y, color = MESS), size = 0.5, alpha = 0.8) +
  geom_sf(data = st_as_sf(china_vect), fill = NA, color = "black", linewidth = 0.4) +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, 
                        name = "MESS Value") +
  labs(title = "MESS Map (Multivariate Environmental Similarity)",
       subtitle = paste0("Negative values (Red) indicate extrapolation risk (", round(pct_extrap, 1), "%)"),
       x = NULL, y = NULL) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        panel.grid = element_blank(), legend.position = "right")

ggsave("Figure_MESS_Risk.png", p_mess, width = 10, height = 8, dpi = 300)

plot_df_risk <- plot_df %>% filter(MESS < 0)

if(nrow(plot_df_risk) > 0) {
  p_factor <- ggplot() +
    geom_sf(data = st_as_sf(china_vect), fill = "gray95", color = "black", linewidth = 0.4) +
    
    geom_point(data = plot_df_risk, aes(x = x, y = y, color = Factor), size = 0.5, alpha = 0.8) +
    
    scale_color_brewer(palette = "Set1", name = "Driving Variable") +
    
    labs(title = "Map of Most Limiting Factors",
         subtitle = "Which variable caused the extrapolation?",
         x = NULL, y = NULL) +
    
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), 
          panel.grid = element_blank(), legend.position = "right")
  
  ggsave("Figure1.png", p_factor, width = 10, height = 8, dpi = 300)
}