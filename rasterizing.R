# Libraries
library(terra)

# load geopackage
rivers <- vect("data/vectors/gme_rivers.gpkg")

# Reproject the rivers to UTM Zone 36S (EPSG 32736
rivers_utm <- project(rivers, "EPSG:32736")

# Create a template raster in meters based on the new UTM extent
# Setting resolution = 30 creates 30m x 30m grid cells
template_raster <- rast(ext(rivers_utm), resolution = 30)
crs(template_raster) <- crs(rivers_utm)

# Calculate Euclidean distance from every cell to the nearest river
distance_raster <- distance(template_raster, rivers_utm)

#  Review the result visually
plot(distance_raster, main = "Distance to Rivers (Meters)")
plot(rivers_utm, add = TRUE, col = "blue")

# Save the final GeoTIFF to your working directory
writeRaster(distance_raster, "data/rasters/distance_rivers.tif", overwrite = TRUE)

# clip
library(terra)

library(terra)

# 1. Load the rasters
prev_river <- rast("data/rasters/rivers.tif")
current_river <- rast("data/rasters/distance_rivers.tif")

# 2. Match CRS if they differ
if (crs(current_river) != crs(prev_river)) {
  current_river <- project(current_river, crs(prev_river))
}

# 3. Apply the 200m buffer filter first (0.2 km)
current_river[current_river > 0.2] <- NA

# 4. CROP: Trim the outer bounding box to match prev_river
cropped_river <- crop(current_river, ext(prev_river))

# 5. MASK: Force it to take the exact irregular pixel shape of prev_river
# Any pixel that is NA in prev_river will now become NA in your final output
final_shaped_river <- mask(cropped_river, prev_river)

# 6. Save the perfectly shaped raster
writeRaster(final_shaped_river, "data/rasters/current_river_final_shape.tif", overwrite = TRUE)
