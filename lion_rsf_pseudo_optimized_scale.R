library(sf)
library(terra)
library(data.table)
library(tidyverse)
library(ggdist)
library(ggbeeswarm)

input_dir  <- "data/collar_data"
output_dir <- "output"
scales_to_test <- seq(100, 2000, by = 100)

z1 <- function(x) as.numeric(scale(x))

make_focal <- function(r, sc, res_m) {
  k <- round(sc / res_m)
  if (k %% 2 == 0) k <- k + 1
  focal(r, w = matrix(1, k, k), fun = mean, na.rm = TRUE)
}

fit_aic_by_scale <- function(df, cov, scales_to_test) {
  map_dfr(scales_to_test, function(sc) {
    nm <- paste0(cov, "_", sc)
    m <- glm(case ~ z1(df[[nm]]), data = df, family = binomial, weights = w)
    tibble(scale = sc, AIC = AIC(m))
  }) %>%
    mutate(dAIC = AIC - min(AIC), covariate = cov)
}

roads_raw      <- rast("data/rasters/distance_roads.tif")
rivers_raw     <- rast("data/rasters/distance_rivers.tif")
dem_raw        <- rast("data/rasters/mara_tpi_3x3.tif")
settlement_raw <- rast("data/rasters/distance_settlements.tif")

master_grid <- rivers_raw
master_crs  <- crs(master_grid)
res_m       <- res(master_grid)[1]

roads_base      <- project(roads_raw, master_grid, method = "bilinear")
rivers_base     <- rivers_raw
dem_base        <- project(dem_raw, master_grid, method = "bilinear")
settlement_base <- project(settlement_raw, master_grid, method = "bilinear")

raster_scale_list <- c(
  setNames(lapply(scales_to_test, \(sc) make_focal(roads_base, sc, res_m)), paste0("roads_", scales_to_test)),
  setNames(lapply(scales_to_test, \(sc) make_focal(rivers_base, sc, res_m)), paste0("rivers_", scales_to_test)),
  setNames(lapply(scales_to_test, \(sc) make_focal(settlement_base, sc, res_m)), paste0("settlement_", scales_to_test)),
  list(dem_base = dem_base)
)

scale_stack <- rast(raster_scale_list)

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
scale_selection_results <- list()
optimized_rsf_results <- list()

set.seed(42)

for (file_path in csv_files) {
  animal_id <- tools::file_path_sans_ext(basename(file_path))
  message("--- Processing animal: ", animal_id, " ---")
  
  gpkg_path <- file.path(output_dir, paste0(animal_id, "_AKDE_HomeRange_est.gpkg"))
  if (!file.exists(gpkg_path)) {
    warning("missing geopackage. skipping.")
    next
  }
  
  hr_sf <- read_sf(gpkg_path) %>% filter(grepl("95% est", name))
  hr_vect <- vect(hr_sf)
  if (crs(hr_vect) != master_crs) hr_vect <- project(hr_vect, master_crs)
  
  df <- fread(file_path)
  pts <- vect(df, geom = c("location.long", "location.lat"), crs = "EPSG:4326")
  pts <- project(pts, master_crs)
  
  inside <- is.related(pts, hr_vect, "intersects")
  pts <- pts[inside, ]
  if (nrow(pts) == 0) {
    warning("zero points inside HR. skipping.")
    next
  }
  
  used_env <- extract(scale_stack, pts)
  used_env$case <- 1
  used_env$w <- 1
  
  hr_clip <- mask(crop(scale_stack, hr_vect), hr_vect)
  
  avail_pts <- spatSample(hr_vect, size = nrow(used_env) * 20, method = "random")
  avail_env <- extract(hr_clip, avail_pts)
  avail_env$case <- 0
  avail_env$w <- 5000
  
  rsf_all <- bind_rows(used_env, avail_env) %>%
    select(-ID) %>%
    na.omit()
  
  best_scales <- list()
  for (cov in c("roads", "rivers", "settlement")) {
    tbl <- fit_aic_by_scale(rsf_all, cov, scales_to_test) %>% mutate(animal_id = animal_id)
    scale_selection_results[[paste0(animal_id, "_", cov)]] <- tbl
    best_scales[[cov]] <- tbl$scale[which.min(tbl$AIC)]
  }
  
  optimized_df <- tibble(
    case = rsf_all$case,
    w = rsf_all$w,
    z_roads = z1(rsf_all[[paste0("roads_", best_scales$roads)]]),
    z_rivers = z1(rsf_all[[paste0("rivers_", best_scales$rivers)]]),
    z_settlement = z1(rsf_all[[paste0("settlement_", best_scales$settlement)]]),
    z_dem = z1(rsf_all$dem_base)
  )
  
  final_model <- glm(case ~ z_roads + z_rivers + z_settlement + z_dem,
                     data = optimized_df, family = binomial, weights = w)
  
  coefs <- summary(final_model)$coefficients
  terms <- c("z_roads", "z_rivers", "z_settlement", "z_dem")
  vars  <- c("roads", "rivers", "settlement", "dem")
  
  beta <- sapply(terms, \(x) coefs[x, "Estimate"])
  se   <- sapply(terms, \(x) coefs[x, "Std. Error"])
  pval <- sapply(terms, \(x) coefs[x, "Pr(>|z|)"])
  
  lwr <- beta - 1.96 * se
  upr <- beta + 1.96 * se
  
  optimized_rsf_results[[animal_id]] <- tibble(
    animal_id = animal_id,
    covariate = vars,
    selected_scale = c(best_scales$roads, best_scales$rivers, best_scales$settlement, 30),
    beta = round(beta, 3),
    beta_lwr = round(lwr, 3),
    beta_upr = round(upr, 3),
    odds_ratio = round(exp(beta), 3),
    or_lwr = round(exp(lwr), 3),
    or_upr = round(exp(upr), 3),
    p_value = round(pval, 3)
  )
}

final_scale_profiles <- bind_rows(scale_selection_results)
write.csv(final_scale_profiles, file.path(output_dir, "raster_scale_aic_profiles.csv"), row.names = FALSE)

final_coef_summary <- bind_rows(optimized_rsf_results)
write.csv(final_coef_summary, file.path(output_dir, "optimized_multi_scale_coefficients.csv"), row.names = FALSE)

ggplot(final_scale_profiles, aes(x = factor(scale), y = dAIC, group = animal_id, color = animal_id)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~covariate, scales = "free_y") +
  labs(title = "Scale Optimization Profiles Across Covariates",
       x = "Neighborhood Smoothing Buffer Scale (meters)",
       y = "Delta AIC",
       color = "Lion ID") +
  theme_bw(base_size = 12)

ggplot(final_coef_summary, aes(x = odds_ratio, y = covariate, color = animal_id)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = or_lwr, xmax = or_upr), height = 0.2,
                 position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        axis.text = element_text(face = "bold", size = 14, color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(colour = "black"),
        legend.position = "top") +
  labs(title = "", x = "Odds Ratio", y = "", color = "Lion Identity")

individual <- read.csv("output/optimized_multi_scale_coefficients.csv")

population_summary <- individual %>%
  group_by(covariate) %>%
  summarize(
    n_animals = n_distinct(animal_id),
    mean_beta = mean(beta),
    sd_beta = sd(beta),
    se_beta = sd_beta / sqrt(n_animals),
    t_value = mean_beta / se_beta,
    p_value = 2 * (1 - pt(abs(t_value), df = n_animals - 1)),
    ci_lwr = mean_beta - qt(0.975, df = n_animals - 1) * se_beta,
    ci_upr = mean_beta + qt(0.975, df = n_animals - 1) * se_beta,
    odds_ratio = exp(mean_beta),
    .groups = "drop"
  ) %>%
  mutate(
    Significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    or_lwr = exp(ci_lwr),
    or_upr = exp(ci_upr)
  )

my_colors <- c("#D35400", "#6A4A3C", "#0F65A1", "#2ECC71")

if (dev.cur() > 1) dev.off()

p <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  ggdist::stat_halfeye(
    data = individual,
    aes(x = odds_ratio, y = covariate, fill = covariate),
    adjust = 1,
    scale = 0.35,
    color = NA,
    position = position_nudge(y = 0.35),
    show.legend = FALSE
  ) +
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
  ggbeeswarm::geom_quasirandom(
    data = individual,
    aes(x = odds_ratio, y = covariate, color = covariate),
    orientation = "y",
    width = 0.08,
    alpha = 0.4,
    size = 2.5,
    show.legend = FALSE
  ) +
  geom_errorbarh(
    data = population_summary,
    aes(xmin = or_lwr, xmax = or_upr, y = covariate),
    height = 0.08,
    color = "black",
    linewidth = 1.0
  ) +
  geom_point(
    data = population_summary,
    aes(x = odds_ratio, y = covariate),
    color = "black",
    fill = "black",
    shape = 21,
    size = 3.0,
    stroke = 1.0
  ) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_log10(breaks = c(0.4, 0.5, 1.0, 1.5, 2.0, 2.5), limits = c(0.3, 2.7)) +
  labs(x = "Odds Ratio", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
  )

print(p)
