# Packages
library(sf)
library(raster)
library(data.table)
library(ctmm)
library(tidyverse)


# Analysis
input_dir <- "data/collar_data"
output_dir <- "output"

# list all csv files to process
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

# create an empty list to store home range size summaries
hr_summary_list <- list()

# loop through all the files
for (file_path in csv_files) {
  
  # get animal id from the filename
  animal_id <- tools::file_path_sans_ext(basename(file_path))
  message("processing: ", animal_id)
  
  # load data and create telemetry object
  df <- fread(file_path)
  l <- as.telemetry(df)
  
  # fit the movement model
  svf <- variogram(l)
  
  # save variogram as png and overwrite if it exists
  var_img_path <- file.path(output_dir, paste0(animal_id, "_variogram.png"))
  if (file.exists(var_img_path)) file.remove(var_img_path)
  
  png(var_img_path, width = 800, height = 600)
  plot(svf)
  title(paste("variogram -", animal_id))
  dev.off()
  
  guess <- ctmm.guess(l, variogram=svf, interactive=FALSE)
  fit <- ctmm.select(l, guess, trace=2)
  
  # calculate the ud
  ud <- akde(l, fit, weights=TRUE)
  
  # save home range plot as png and overwrite if it exists
  hr_img_path <- file.path(output_dir, paste0(animal_id, "_homerange_plot.png"))
  if (file.exists(hr_img_path)) file.remove(hr_img_path)
  
  png(hr_img_path, width = 800, height = 600)
  plot(ud)
  title(paste("akde home range -", animal_id))
  dev.off()
  
  # extract home range size metrics
  ud_sum <- summary(ud, units = TRUE)
  
  # find the row index containing the word area
  area_row_index <- grep("area", rownames(ud_sum$CI))
  area_row <- ud_sum$CI[area_row_index, ]
  
  # convert values to text before stripping units
  clean_low  <- as.numeric(gsub("[^0-9.-]", "", as.character(area_row[1])))
  clean_est  <- as.numeric(gsub("[^0-9.-]", "", as.character(area_row[2])))
  clean_high <- as.numeric(gsub("[^0-9.-]", "", as.character(area_row[3])))
  
  # save metrics into a data frame row
  hr_metrics <- data.frame(
    animal_id = animal_id,
    hr_size_sq_km = clean_est,
    hr_low_ci_sq_km = clean_low,
    hr_high_ci_sq_km = clean_high,
    stringsAsFactors = FALSE
  )
  
  hr_summary_list[[animal_id]] <- hr_metrics
  
  # convert ud object to sf
  ud_sp <- SpatialPolygonsDataFrame.UD(ud)
  ud_sf <- st_as_sf(ud_sp)
  
  # isolate only the official est home range polygon row
  ud_est_sf <- ud_sf[grepl("95% est", ud_sf$name), ]
  
  # export the est polygon directly to the output folder as a geopackage
  gpkg_filename <- file.path(output_dir, paste0(animal_id, "_AKDE_HomeRange_est.gpkg"))
  st_write(ud_est_sf, gpkg_filename, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
  
  message("finished: ", animal_id, "\n")
}

# compile and export csv to the output folder
final_hr_table <- bind_rows(hr_summary_list)
write.csv(final_hr_table, file.path(output_dir, "all_lions_hr_size_summary.csv"), row.names = FALSE)

message("all processing complete!")
