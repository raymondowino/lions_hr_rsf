# Libraries
library(elevatr)
library(terra)
library(sf)

# load geopackage

## river
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

## road
roads <- vect("data/vectors/gme_roads.gpkg")

# Reproject the rivers to UTM Zone 36S (EPSG 32736
roads_utm <- project(roads, "EPSG:32736")

# Create a template raster in meters based on the new UTM extent
# Setting resolution = 30 creates 30m x 30m grid cells
template_raster <- rast(ext(roads_utm), resolution = 30)
crs(template_raster) <- crs(roads_utm)

# Calculate Euclidean distance from every cell to the nearest river
distance_raster <- distance(template_raster, roads_utm)

#  Review the result visually
plot(distance_raster, main = "Distance to Roads (Meters)")
plot(roads_utm, add = TRUE, col = "blue")

# Save the final GeoTIFF to your working directory
writeRaster(distance_raster, "data/rasters/distance_roads.tif", overwrite = TRUE)


## settlement 
settlements <- vect("data/vectors/gme_settlments.gpkg")

# Reproject the settlements to UTM Zone 36S (EPSG 32736)
settlements_utm <- project(settlements, "EPSG:32736")

# Filter out non-settlement attributes to isolate true human footprints
settlement_types <- c("Settlement_Boma", "Settlement_Campsite", "Town",
                      "Ranger Post","School" ,"Health Facility" ,"Ranger Station" ,
                      "Police Post" ,"Settlement_Safari Camp", "Settlement_Hotel/Lodge", "Village")
true_settlements_utm <- subset(settlements_utm, settlements_utm$type %in% settlement_types)

# Create a template raster based on the filtered settlements extent
# Note: Adding a buffer (e.g., ext(true_settlements_utm) + 5000) prevents harsh clipping at the edges
template_raster <- rast(ext(true_settlements_utm), resolution = 30)
crs(template_raster) <- crs(true_settlements_utm)

# Calculate Euclidean distance from every cell to the nearest true settlement
distance_raster <- distance(template_raster, true_settlements_utm)

# Convert distance values from meters to kilometers for easier modeling interpretation
#distance_raster_km <- distance_raster / 1000

# Review the result visually
plot(distance_raster, main = "Distance to settlements (Kilometers)")
plot(true_settlements_utm, add = TRUE, col = "red", cex = 0.5)

# Save the final GeoTIFF to your working directory
writeRaster(distance_raster, "data/rasters/distance_settlements.tif", overwrite = TRUE)


### DEM (AWS Alternative Source)

# Load raster template
gme_aoi <- rast("data/rasters/dem_30m.tif")

# Project template raster to UTM Zone 36S
gme_aoi_utm <- project(gme_aoi, "EPSG:32736")

# Convert spatial extent object to sf format for the API query
gme_sf <- st_as_sf(as.polygons(ext(gme_aoi_utm), crs = crs(gme_aoi_utm)))

# Download global seamless DEM matching your extent bounds
# z = 11 provides an excellent native ~30m landscape resolution scale
dem_raw <- get_elev_raster(locations = gme_sf, z = 11, src = "aws")
dem_wgs84 <- rast(dem_raw)

# Crop and project new DEM data onto your 30m UTM 36S template grid
mara_dem_30m <- project(dem_wgs84, gme_aoi_utm, method = "bilinear")
mara_dem_30m <- mask(mara_dem_30m, gme_aoi_utm)

# Save clean DEM raster
writeRaster(mara_dem_30m, "data/rasters/mara_elevation_30m.tif", overwrite = TRUE)

# Generate TPI with 3x3 window
mara_tpi_3x3 <- terrain(mara_dem_30m, v = "TPI")

# Save TPI raster
writeRaster(mara_tpi_3x3, "data/rasters/mara_tpi_3x3.tif", overwrite = TRUE)

# Plot outputs to confirm the Mara Triangle is fully covered
par(mfrow = c(1, 2))
plot(mara_dem_30m, main = "Elevation (m)")
plot(mara_tpi_3x3, main = "TPI (3x3 Window)")

# reset plotting
par(mfrow = c(1, 1))
