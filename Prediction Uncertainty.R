packages <- c("terra", "sf", "ranger", "ggplot2", "viridis", "dplyr", "geodata", "scales", "ggspatial")
for(p in packages){
  if(!require(p, character.only = T)) install.packages(p)
  library(p, character.only = T)
}

model_features <- c("OC", "pH", "Clay", "CEC", "logKow", "lnCe", "Temperature", "Ratio")

final_m_et_quant <- ranger(
  formula = Y_target ~ ., 
  data = train_final_df[, c(model_features, "Y_target")],
  num.trees = 1000,        
  mtry = 3,                
  min.node.size = 5,       
  quantreg = TRUE,         
  num.threads = 4,
  verbose = TRUE
)

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
} else if(file.exists(raster_path)) {
  global_id_raster <- rast(raster_path)
} else {
  global_id_raster <- rast(ext=ext(china_vect), res=0.1) 
  values(global_id_raster) <- sample(1:10000, ncell(global_id_raster), replace=TRUE)
  names(global_id_raster) <- "MU_GLOBAL"
}

csv_path <- "C:/Users/电击小子/Desktop/Default Dataset.csv"
soil_data_raw <- read.csv(csv_path)
soil_cols <- c("MU_GLOBAL", "OC", "pH", "Clay", "CEC")
soil_attributes <- soil_data_raw[, c(1, 2:5)]
colnames(soil_attributes) <- soil_cols

soil_attributes <- soil_attributes %>%
  filter(complete.cases(.)) %>%  
  group_by(MU_GLOBAL) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  ungroup()

china_proj <- project(china_vect, crs(global_id_raster))
china_rast <- crop(global_id_raster, china_proj)
china_rast <- mask(china_rast, china_proj)
china_pixels <- as.data.frame(china_rast, xy = TRUE, na.rm = TRUE)
colnames(china_pixels)[3] <- "MU_GLOBAL"
pred_data <- left_join(china_pixels, soil_attributes, by = "MU_GLOBAL")
pred_data <- na.omit(pred_data)

pred_input <- pred_data

pred_input$logKow <- 3.36
pred_input$lnCe <- 4.60517
pred_input$Temperature <- 298
pred_input$Ratio <- 0.01

pred_features <- pred_input[, model_features]

batch_size <- 20000 
n_rows <- nrow(pred_features)
n_batches <- ceiling(n_rows / batch_size)

res_uncertainty_ln <- numeric(n_rows)

pb <- txtProgressBar(min = 0, max = n_batches, style = 3)

for(i in 1:n_batches) {
  idx_start <- (i - 1) * batch_size + 1
  idx_end <- min(i * batch_size, n_rows)
  idx <- idx_start:idx_end
  
  batch_data <- pred_features[idx, ]
  
  preds_quant <- predict(final_m_et_quant, 
                         data = batch_data, 
                         type = "quantiles", 
                         quantiles = c(0.05, 0.95))$predictions
  
  res_uncertainty_ln[idx] <- preds_quant[, 2] - preds_quant[, 1]
  
  setTxtProgressBar(pb, i)
  gc() 
}
close(pb)

plot_idx <- 1:n_rows
plot_df <- data.frame(
  x = pred_data$x[plot_idx],
  y = pred_data$y[plot_idx],
  Uncertainty = res_uncertainty_ln[plot_idx]
)

val_min <- min(plot_df$Uncertainty, na.rm=TRUE)
val_max <- max(plot_df$Uncertainty, na.rm=TRUE)


china_sf <- st_as_sf(china_vect)

p1 <- ggplot() +
  geom_tile(data = plot_df, aes(x = x, y = y, fill = Uncertainty), color = NA) +
  
  geom_sf(data = china_sf, fill = NA, color = "black", linewidth = 0.5) +
  
  scale_fill_viridis(
    option = "plasma",   
    name = "Uncertainty\n(lnQe 90% PI)", 

    guide = guide_colorbar(
      barheight = unit(12, "cm"),   
      barwidth = unit(0.6, "cm"),   
      frame.colour = "black",        
      ticks.colour = "black",        
      frame.linewidth = 0.5,        
      ticks.linewidth = 0.5,        
      title.position = "top",        
      title.hjust = 0.5             
    )
  ) +

  annotation_scale(
    location = "bl",
    style = "bar",
    bar_cols = c("black", "white"),
    width_hint = 0.5, 
    height = unit(0.4, "cm"),
    line_width = 1,
    text_cex = 1.2,
    text_pad = unit(0.2, "cm"),
    pad_x = unit(1.0, "cm"),
    pad_y = unit(0.8, "cm")
  ) +
  
  annotation_north_arrow(
    location = "tl",            
    which_north = "true",
    style = north_arrow_fancy_orienteering, 
    height = unit(1.5, "cm"),
    width = unit(1.5, "cm"),
    pad_x = unit(0.5, "cm"),    
    pad_y = unit(0.5, "cm")     
  ) +

  labs(
    title = "Prediction Uncertainty of Phenanthrene",
    subtitle = paste0("Raw ln(Qe) 90% PI Width (Range: ", round(val_min, 2), " - ", round(val_max, 2), ")"),
    x = NULL, y = NULL
  ) +
  
  theme_bw() + 
  theme(
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    legend.margin = margin(l = 10),
    plot.title = element_text(face="bold", size=16, hjust=0.5),
    plot.subtitle = element_text(size=12, hjust=0.5, color="grey30")
  )

filename <- file.path(output_dir, "China_Map_Uncertainty_Phenanthrene_lnQe_Plasma.png")

ggsave(filename, p1, width = 10, height = 8, dpi = 300)