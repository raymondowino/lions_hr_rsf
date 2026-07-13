# libraries
library(tidyverse)
library(ggdist)
library(ggbeeswarm)
library(glmmTMB)

# data
all_data <- read.csv("output/all_combined_points.csv")
head(all_data)
table(all_data$case,all_data$animal_id)

summary(all_data)


# z-score standardize elevation, slope, ET, HLI, and transformed 
# aspect

all_data$zroads <- scale(all_data$roads)
all_data$zrivers <- scale(all_data$rivers)
all_data$zdem <- scale(all_data$dem)
all_data$zsettlement <- scale(all_data$settlement)

# check again
head(all_data)

# Next, lets add some large weights to our available points. Lets try weights
# of 10000 just because we have relatively few available points.

all_data$w <- ifelse(all_data$case==1,1,10000)

table(all_data$case,all_data$w)

# fit rsf still with pseudoreplication
lion_bin_GLM <- glm(case ~ zroads + zrivers + zdem + zsettlement,
                    weights = w,
                    family = binomial,
                    data = all_data)

summary(lion_bin_GLM)

# fit rsf still with pseudoreplication (control)
lion_rsf <- glmmTMB(case ~ zroads + zrivers + zdem + zsettlement + (1 | animal_id), 
                    weights = w, 
                    family = binomial, 
                    data = all_data)

# View the summary
summary(lion_rsf)

# pull random effect
(lion_rsf_ranef <- as.data.frame(coef(lion_rsf)$cond$animal_id))

# comapare mean and sd to guassian distribution
(mean_RI <- mean(lion_rsf_ranef$`(Intercept)`))
(sd_RI <- sd(lion_rsf_ranef$`(Intercept)`))

# distibution form actual rndom effect
(est_mean_RI <- fixef(lion_rsf)$cond['(Intercept)'])
(est_sd_RI <- summary(lion_rsf)$varcor$cond$animal_id[1]^0.5)

# plot a histogram of our individual random intercepts
hist(lion_rsf_ranef$`(Intercept)`,
     breaks = seq(-13,-8,by=0.1),
     freq = FALSE)

# Add 2 density curves: 
# Red = Gaussian distribution estimated by glmmTMB random intercepts
# Blue = Gaussian distribution based on sample mean and SD of random intercepts
curve(dnorm(x, mean_RI, sd_RI), from = -13, to = 1, 
      add = T, col = "blue", lwd = 2)
curve(dnorm(x, est_mean_RI, est_sd_RI), from = -13, to = 1, 
      add = T, col = "red", lwd = 2)
legend("topleft",c("Model Estimates","Derived"),
       lwd = 2, col = c("red","blue"))


# Muff et al. (2020) recommends making the distribution of the 
# random intercept to be very large
lion_rsf_fixed <- glmmTMB(case ~ zroads + zrivers + zdem + zsettlement + (1 | animal_id), 
                    weights = w, 
                    family = binomial, 
                    data = all_data,
                    doFit=FALSE)

# fix the standard deviation of the first random term, which is 
# the (1|id) component in the above model equation
lion_rsf_fixed$parameters$theta[1] = log(1e3) 


#
curve(dnorm(x, 0, log(1e3)), from = -20, to = 20, 
      col = "green", lwd = 2, add = T)

#
curve(dnorm(x, 0, log(1e3)), from = -20, to = 20, 
      col = "green", lwd = 2)


#

lion_rsf_fixed$mapArg = list(theta=factor(c(NA)))

system.time(RI_only_fixed_mod <- glmmTMB:::fitTMB(lion_rsf_fixed))

# Lets look at our summary output

summary(RI_only_fixed_mod)



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


