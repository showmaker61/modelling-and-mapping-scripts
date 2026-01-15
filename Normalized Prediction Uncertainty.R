# ==============================================================================
# 菲 (Phenanthrene) 归一化不确定性图 (手动分母版)
# 核心逻辑：
# 1. 计算绝对不确定性 (Q95 - Q05)
# 2. [关键] 除以手动设定的分母 (相对吸附量图的有效跨度)
# ==============================================================================

# --- 0. 用户配置区 (请修改这里！) ---
MANUAL_DENOMINATOR <- ***  

packages <- c("terra", "sf", "ranger", "ggplot2", "viridis", "dplyr", "geodata", "scales")
for(p in packages){
  if(!require(p, character.only = T)) install.packages(p)
  library(p, character.only = T)
}

model_features <- c("OC", "pH", "Clay", "CEC", "logKow", "lnCe", "Temperature", "Ratio")

if(!exists("final_m_et_quant")) {
  final_m_et_quant <- ranger(
    formula = Y_target ~ ., 
    data = train_final_df[, c(model_features, "Y_target")],
    num.trees = 1000, mtry = 3, min.node.size = 5,splitrule = "extratrees",
    quantreg = TRUE,  
    num.threads = 4, verbose = TRUE
  )
}

china_shp_path <- "C:/Users/电击小子/Desktop/china.shp" 
if(!exists("china_vect")) {
  if (file.exists(china_shp_path)) {
    china_vect <- vect(china_shp_path)
  } else {
    china_vect <- gadm(country = "CHN", level = 1, path = tempdir())
  }
}

raster_path <- "C:/Users/电击小子/Desktop/raster/HWSD2.bil"
if(exists("global_id_raster")){
  # pass
} else if(file.exists(raster_path)) {
  global_id_raster <- rast(raster_path)
} else {
  global_id_raster <- rast(ext=ext(china_vect), res=0.1) 
  values(global_id_raster) <- sample(1:10000, ncell(global_id_raster), replace=TRUE)
  names(global_id_raster) <- "MU_GLOBAL"
}

csv_path <- "C:/Users/电击小子/Desktop/Default Dataset.csv"
soil_data_raw <- read.csv(csv_path)
soil_attributes <- soil_data_raw[, c(1, 2:5)]
colnames(soil_attributes) <- c("MU_GLOBAL", "OC", "pH", "Clay", "CEC")

soil_attributes <- soil_attributes %>%
  filter(complete.cases(.)) %>%  
  group_by(MU_GLOBAL) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

china_proj <- project(china_vect, crs(global_id_raster))
china_rast <- crop(global_id_raster, china_proj)
china_rast <- mask(china_rast, china_proj)
china_pixels <- as.data.frame(china_rast, xy = TRUE, na.rm = TRUE)
colnames(china_pixels)[3] <- "MU_GLOBAL"

pred_data <- left_join(china_pixels, soil_attributes, by = "MU_GLOBAL")
pred_data <- na.omit(pred_data)


pred_input <- pred_data

pred_input$logKow <- ***
pred_input$lnCe <- 4.60517
pred_input$Temperature <- 298
pred_input$Ratio <- 0.01

pred_features <- pred_input[, model_features]

batch_size <- 20000 
n_rows <- nrow(pred_features)
n_batches <- ceiling(n_rows / batch_size)

res_uncertainty <- numeric(n_rows) 

pb <- txtProgressBar(min = 0, max = n_batches, style = 3)

for(i in 1:n_batches) {
  idx_start <- (i - 1) * batch_size + 1
  idx_end <- min(i * batch_size, n_rows)
  idx <- idx_start:idx_end
  
  preds_quant <- predict(final_m_et_quant, 
                         data = pred_features[idx, ], 
                         type = "quantiles", 
                         quantiles = c(0.05, 0.95))$predictions

  res_uncertainty[idx] <- preds_quant[, 2] - preds_quant[, 1]
  
  setTxtProgressBar(pb, i)
}
close(pb)

norm_uncertainty <- res_uncertainty / MANUAL_DENOMINATOR


plot_df <- data.frame(
  x = pred_data$x,
  y = pred_data$y,
  Val = norm_uncertainty
)

print(summary(plot_df$Val))
cols_contrast <- c("#313695", "#4575b4", "#74add1", "#abd9e9", 
                   "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", 
                   "#f46d43", "#d73027", "#a50026")

p1 <- ggplot() +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = Val), color = NA) +
  geom_sf(data = st_as_sf(china_vect), fill = NA, color = "black", linewidth = 0.4) +
  
  scale_fill_gradientn(
    colors = cols_contrast, 
    name = "Normalized\nUncertainty",
    limits = c(0, 1),
    oob = scales::squish, 
    guide = guide_colorbar(
      barheight = unit(15, "cm"), 
      barwidth = unit(0.6, "cm"), 
      frame.colour = "black", 
      ticks.colour = "black",
      label.theme = element_text(size = 10)
    )
  ) +
  
  labs(title = "Normalized Prediction Uncertainty of Phenanthrene",
       subtitle = sprintf("Ratio of (Q95-Q05) to Adsorption Range (%.2f)", MANUAL_DENOMINATOR),
       x = NULL, y = NULL) +
  
  theme_bw() + 
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face="bold", size=14)
  )
filename <- "Uncertainty_Normalized.png"
ggsave(filename, p1, width = 10, height = 8, dpi = 300)