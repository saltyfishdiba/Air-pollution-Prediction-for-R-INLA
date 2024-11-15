library(sf)
library(raster)
library(sp) 
library(dplyr)
library(data.table)
library(gstat)
library(ggplot2)
library(viridis) 

shapefile_path <- "D:/biostat article/st.louis project/LimaBorder/LimaBorder.shp"
lima_shapefile <- st_read(shapefile_path)

lima_shapefile_singlepart <- st_cast(lima_shapefile, "POLYGON")

target_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
lima_shapefile_transformed <- st_transform(lima_shapefile_singlepart, crs = target_crs)

csv_file_path <- "D:/biostat article/st.louis project/selected_columns_data.csv"
csv_data <- fread(csv_file_path)

csv_data$Year <- 2017

# missing x and y coordinates with the mean values from surrounding points
csv_data <- csv_data %>%
  mutate(x = ifelse(is.na(x), mean(x, na.rm = TRUE), x),
         y = ifelse(is.na(y), mean(y, na.rm = TRUE), y))


csv_data$mean_PM25 <- as.numeric(csv_data$mean_PM25)

# csv_data to SpatialPointsDataFrame
coordinates(csv_data) <- ~x + y
proj4string(csv_data) <- CRS(target_crs)

# Create a blank raster grid 
raster_template <- raster(extent(lima_shapefile_transformed), 
                          crs = target_crs, 
                          resolution = c(1000, 1000))  # Adjust resolution as needed

#Perform IDW interpolation
idw_result <- gstat::idw(mean_PM25 ~ 1, csv_data, newdata = as(raster_template, 'SpatialGrid'))

# IDW output to raster
pm25_raster <- raster(idw_result)

# Mask the interpolated raster to the boundaries of the shapefile
pm25_raster <- mask(pm25_raster, st_as_sf(lima_shapefile_transformed))

pm25_df <- as.data.frame(rasterToPoints(pm25_raster))
colnames(pm25_df) <- c("x", "y", "PM25")

breaks <- c(0, 10, 15, 20, 24.8, 28.2, 31.5, 35, 38.5, 55)
custom_colors <- colorRampPalette(c("#0000ff", "#0033cc", "#0066ff", "#00cc66", 
                                    "#66ff33", "#ffff00", "#ffcc00", "#ff9900", 
                                    "#ff6600"))(length(breaks) - 1)

ggplot() +
  geom_tile(data = pm25_df, aes(x = x, y = y, fill = cut(PM25, breaks = breaks))) +
  scale_fill_manual(values = custom_colors, name = "PM2.5") +
  geom_sf(data = lima_shapefile_transformed, fill = NA, color = "black") +  # Add boundary lines
  labs(title = "Annual PM2.5 - 2017") +
  theme_minimal()

ggsave("D:/biostat_article/PM2.5_Interpolated.png")

