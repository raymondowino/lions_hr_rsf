# packages 
library(sf) 
library(terra) 
library(data.table) 
library(tidyverse) 

# analysis 
input_dir <- "data/collar_data" 
output_dir <- "output" 

# load raw maps 
roads_raw <- rast("data/rasters/distance_roads.tif") 
rivers_raw <- rast("data/rasters/distance_rivers.tif") 
dem_raw <- rast("data/rasters/mara_tpi_3x3.tif")
settlement_raw <- rast("data/rasters/distance_settlements.tif") 

# establish the master baseline grid reference 
master_crs <- crs(rivers_raw) 
master_grid <- rivers_raw 

message("aligning and projecting rasters...") 
roads_aligned <- terra::project(roads_raw, master_grid, method = "bilinear") 
dem_aligned <- terra::project(dem_raw, master_grid, method = "bilinear") 
settlement_aligned <- terra::project(settlement_raw, master_grid, method = "bilinear") 

# bundle into a single stack 
raster_stack <- c(roads_aligned, rivers_raw, dem_aligned, settlement_aligned) 
names(raster_stack) <- c("roads", "rivers", "dem", "settlement") 

# list tracking datasets 
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE) 
all_used_list <- list() 
all_available_list <- list() 
rsf_results_list <- list() 

# set seed for reproducible pseudo-random background background sampling
set.seed(42) 

# loop through each animal file 
for (file_path in csv_files) { 
  animal_id <- tools::file_path_sans_ext(basename(file_path)) 
  message("processing rsf for: ", animal_id) 
  
  gpkg_path <- file.path(output_dir, paste0(animal_id, "_AKDE_HomeRange_est.gpkg")) 
  if (!file.exists(gpkg_path)) { 
    warning("missing geopackage for ", animal_id, ". skipping.") 
    next 
  } 
  
  # load home range layers safely 
  lion_hr_sf <- read_sf(gpkg_path) 
  lion_hr_sf_est <- lion_hr_sf[grepl("95% est", lion_hr_sf$name), ] 
  lion_hr_vect <- vect(lion_hr_sf_est) 
  
  if (crs(lion_hr_vect) != master_crs) { 
    lion_hr_vect <- project(lion_hr_vect, master_crs) 
  } 
  
  # load point data 
  df <- fread(file_path) 
  lion_points <- vect(df, geom = c("location.long", "location.lat"), crs = "+proj=longlat +datum=WGS84") 
  lion_points <- project(lion_points, master_crs) 
  
  # strip points outside the polygon shape boundary 
  points_inside <- terra::is.related(lion_points, lion_hr_vect, "intersects") 
  lion_points_clean <- lion_points[points_inside, ] 
  
  if (length(lion_points_clean) == 0) { 
    warning("zero tracking points found inside home range for ", animal_id, ". skipping.") 
    next 
  } 
  
  # used data extractions 
  used_env <- terra::extract(raster_stack, lion_points_clean) 
  used_env$ID <- NULL 
  used_env$animal_id <- animal_id 
  used_env$case <- 1 
  used_env$w <- 1 
  
  # background available point generation 
  raster_cropped <- terra::crop(raster_stack, lion_hr_vect) 
  raster_clipped <- terra::mask(raster_cropped, lion_hr_vect) 
  available_env <- spatSample(lion_hr_vect, size = nrow(used_env) * 20, method = "random") 
  available_values <- terra::extract(raster_clipped, available_env) 
  available_values$ID <- NULL 
  available_values$animal_id <- animal_id 
  available_values$case <- 0 
  available_values$w <- 5000 
  
  # export single combined used/available file per lion ID
  individual_combined <- rbind(used_env, available_values)
  write.csv(individual_combined, file.path(output_dir, paste0(animal_id, "_used_available_points.csv")), row.names = FALSE)
  
  # bundle individual vectors together 
  all_used_list[[animal_id]] <- used_env 
  all_available_list[[animal_id]] <- available_values 
  
  rsf_data <- rbind(used_env, available_values) 
  rsf_data_clean <- na.omit(rsf_data) 
  
  if (nrow(rsf_data_clean) == 0) { 
    warning("empty data cells for ", animal_id, ". skipping.") 
    next 
  } 
  
  # generate z-scored metrics 
  rsf_data_clean$z_roads <- scale(rsf_data_clean$roads) 
  rsf_data_clean$z_rivers <- scale(rsf_data_clean$rivers) 
  rsf_data_clean$z_dem <- scale(rsf_data_clean$dem) 
  rsf_data_clean$z_settlement <- scale(rsf_data_clean$settlement) 
  
  # model fitting 
  rsf_model <- glm(case ~ z_roads + z_rivers + z_dem + z_settlement, data = rsf_data_clean, family = binomial, weights = w) 
  m_summary <- summary(rsf_model)$coefficients 
  
  # calculate confidence intervals and round metrics to 3 digits
  covariates <- c("roads", "rivers", "dem", "settlement") 
  
  beta_vals <- sapply(covariates, function(x) as.numeric(m_summary[paste0("z_", x), "Estimate"]))
  se_vals   <- sapply(covariates, function(x) as.numeric(m_summary[paste0("z_", x), "Std. Error"]))
  p_vals    <- sapply(covariates, function(x) as.numeric(m_summary[paste0("z_", x), "Pr(>|z|)"]))
  
  beta_lwr  <- beta_vals - (1.96 * se_vals)
  beta_upr  <- beta_vals + (1.96 * se_vals)
  
  lion_coefs <- data.frame( 
    animal_id   = rep(animal_id, length(covariates)), 
    covariate   = covariates, 
    beta        = round(beta_vals, 3), 
    beta_lwr    = round(beta_lwr, 3), 
    beta_upr    = round(beta_upr, 3), 
    se          = round(se_vals, 3), 
    p_value     = round(p_vals, 3), 
    odds_ratio  = round(exp(beta_vals), 3), 
    or_lwr      = round(exp(beta_lwr), 3), 
    or_upr      = round(exp(beta_upr), 3), 
    row.names   = NULL, 
    stringsAsFactors = FALSE 
  ) 
  rsf_results_list[[animal_id]] <- lion_coefs 
} 

# export master population datasets to disk 
final_used <- bind_rows(all_used_list) 
final_available <- bind_rows(all_available_list) 
all_combined_points <- rbind(final_used, final_available) %>% na.omit() 

# write final outputs
write.csv(all_combined_points, file.path(output_dir, "all_combined_points.csv"), row.names = FALSE) 

final_summary_table <- bind_rows(rsf_results_list) 
write.csv(final_summary_table, file.path(output_dir, "all_lions_rsf_coefficients_summary.csv"), row.names = FALSE)
message("All processing runs complete and files saved successfully!")

# Visualize

ggplot(final_summary_table, aes(x = odds_ratio, y = covariate, color = animal_id)) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") + 
  geom_errorbarh(aes(xmin = or_lwr, xmax = or_upr), height = 0.2, position = position_dodge(width = 0.4)) + 
  geom_point(position = position_dodge(width = 0.4), size = 2)+
  theme_bw(base_size = 12) + 
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5), 
    axis.text = element_text(face = "bold", size = 14, color = "black"), 
    axis.title = element_text(size = 14, face = "bold"), 
    axis.ticks = element_line(colour = "black"), 
    legend.position = "top"
  ) + 
  labs(title = "", x = "Odds Ratio", y = "", color = "Lion Identity")

# Two stage population level rsf
# Load coefficients from the individual results
individual <- read.csv("output/all_lions_rsf_coefficients_summary.csv")

#  Compute two-stage summary with directional logic
population_summary <- individual %>%
  group_by(covariate) %>%
  summarize(
    n_animals  = n_distinct(animal_id),
    mean_beta  = mean(beta),
    sd_beta    = sd(beta),
    se_beta    = sd_beta / sqrt(n_animals),
    t_value    = mean_beta / se_beta,
    p_value    = 2 * (1 - pt(abs(t_value), df = n_animals - 1)),
    ci_lwr     = mean_beta - (qt(0.975, df = n_animals - 1) * se_beta),
    ci_upr     = mean_beta + (qt(0.975, df = n_animals - 1) * se_beta),
    odds_ratio = exp(mean_beta)
  ) %>%
  mutate(
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "NS"
    )
  )

# Format population metrics back to Odds Ratio scale for plotting
plot_population <- population_summary %>%
  mutate(
    or_lwr = exp(ci_lwr),
    or_upr = exp(ci_upr)
  )

# Define vibrant colors matching your template aesthetic
my_colors <- c("#D35400", "#6A4A3C", "#0F65A1", "#2ECC71")

# Clear graphic device beforehand to avoid internal state errors
if (dev.cur() > 1) dev.off()

#  Generate the Raincloud Forest Plot
p <- ggplot() +
  # Null reference line at Odds Ratio = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  
  #  Density slab (half-eye) - Shifted up to clear the boxplot
  ggdist::stat_halfeye(
    data = individual,
    aes(x = odds_ratio, y = covariate, fill = covariate),
    adjust = 1, 
    scale = 0.35,         
    color = NA, 
    position = position_nudge(y = 0.35),
    show.legend = FALSE
  ) +
  
  #  Boxplot - Nudged upward to sit right beside the point data
  geom_boxplot(
    data = individual,
    aes(x = odds_ratio, y = covariate, fill = covariate),
    width = 0.08, 
    color = "black",
    fill = "white",
    linewidth = 0.8, 
    outlier.shape = NA,
    position = position_nudge(y = 0.15),
    show.legend = FALSE
  ) +
  
  #  Jittered individual lion points - Centered exactly on the main line (y = 0)
  ggbeeswarm::geom_quasirandom(
    data = individual,
    aes(x = odds_ratio, y = covariate, color = covariate),
    orientation = "y", 
    width = 0.08, 
    alpha = 0.4, 
    size = 2.5,
    show.legend = FALSE
  ) +
  
  # Population Mean Confidence Interval - Centered exactly on the main line (y = 0)
  geom_errorbarh(
    data = plot_population, 
    aes(xmin = or_lwr, xmax = or_upr, y = covariate), 
    height = 0.08, 
    color = "black", 
    linewidth = 1.0
  ) +
  
  # Population Mean Point 
  geom_point(
    data = plot_population, 
    aes(x = odds_ratio, y = covariate), 
    color = "black", 
    fill = "black",   # Made solid black instead of white
    shape = 21, 
    size = 5.0,       
    stroke = 1.0
  ) +
  
  # Color scaling and layout configuration
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_log10(breaks = c(0.4, 0.5, 1.0, 1.5, 2.0, 2.5), limits = c(0.3, 2.7)) +
  labs(x = "Odds Ratio", y = NULL) +
  
  # Clean, minimalist journal theme layout
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )

# plot
print(p)



