# Load required libraries
library(randomForest)
library(ggplot2)
library(sf)
library(dplyr)

# Step 1: Load Lima city shapefile
shapefile_path <- "D:/biostat article/st.louis project/LimaBorder/LimaBorder.shp"
lima_shapefile <- st_read(shapefile_path)

# Step 2: Transform the coordinate reference system (CRS)
lima_shapefile_transformed <- st_transform(lima_shapefile, crs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")

# Step 3: Load the dataset (PM2.5 data)
data_path <- "D:/biostat article/st.louis project/Prediction_selected_data_2017.csv"
data <- read.csv(data_path)

# Step 4: Preprocess the data
data <- na.omit(data)

# Skip 'x' and 'y' (spatial coordinates) from the variable selection
features <- c("tree", "shrub", "grass", "crop", "builtup", "bare", "snowice", 
              "water", "wetland", "Pop", "roadl", "elevation", "AOD", "PBLH", 
              "HFX", "LH", "CLDFRA", "SWDOWN", "VEGFRA", "T2", "PSFC", "U10", "V10")

# Use the selected features (excluding x and y) as predictors
X <- data[features]
y <- data$PM2_5_DRY

# Step 5: Train the Random Forest model with memory optimizations
set.seed(42)
rf_model <- randomForest(X, y, ntree = 100, mtry = floor(sqrt(ncol(X))), sampsize = floor(0.8 * nrow(X)))

# Step 6: Predict PM2.5 values
data$PM2_5_Pred <- predict(rf_model, newdata = X)

# Step 7: Create spatial points from the predicted data (use x and y for coordinates)
coordinates <- data.frame(lon = data$x, lat = data$y)  # Assuming x and y represent longitude and latitude
points_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = st_crs(lima_shapefile_transformed))

# Step 8: Plot the predictions on the map of Lima city
ggplot() +
  geom_sf(data = lima_shapefile_transformed, fill = "white", color = "black") +
  geom_sf(data = points_sf, aes(color = data$PM2_5_Pred), size = 1) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "PM2.5 Prediction Map for Lima City", color = "PM2.5 Levels") +
  theme_minimal()


