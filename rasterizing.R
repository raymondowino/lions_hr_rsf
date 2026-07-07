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

