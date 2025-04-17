library(tidyverse)
library(sf)
setwd("~/Documents/git/BIOSTAT696")
birth <- readRDS('data/input/birth_imp.RDS')
lat_long <- read_csv('data/input/cluster_locations.csv')[,-1] %>% 
  select(DHSCLUST, LATNUM, LONGNUM)
colnames(lat_long) <- c("clusterid", "lat", "long")
birth <- left_join(birth, lat_long, by = 'clusterid')

clusters <- read_sf("data/input/bdg/bdg_cluster.shp")
country <- read_sf("data/input/bdg/BGD_adm0.shp")


## Plot of clusters
ggplot() +
  geom_sf(data = country, fill = 'gray90', color = 'black') +
  geom_sf(data = clusters, color = "red", size = 1) +
  theme_minimal()



library(spBayes)
library(gstat)
library(sp)

set.seed(2025)
jitter_sd <- 0.03

birth$lat <- birth$lat + rnorm(n = nrow(birth), mean = 0, sd = jitter_sd)
birth$long <- birth$long + rnorm(n = nrow(birth), mean = 0, sd = jitter_sd)


sf_data <- st_as_sf(birth, coords = c("long", "lat"), crs = 4326)  # EPSG:4326 = WGS84 (lat/lon)
sf_proj <- st_transform(sf_data, crs = 32646)

variog <- variogram(birth_weight ~ 1, sf_proj)

# Convert distance to kilometers
variog$dist_km <- variog$dist / 1000

# Variogram
ggplot(variog, aes(x = dist_km, y = gamma)) +
  geom_point(color = "steelblue") +
  geom_line(color = "steelblue") +
  labs(x = "Distance (km)", y = "Semivariance") +
  theme_minimal(base_size = 14)



####### Variogram for the residuals

variog_res <- variogram(birth_weight ~ mother_bmi + sex_of_child + wealth_index, sf_proj)

# Convert distance to kilometers
variog_res$dist_km <- variog$dist / 1000

# Variogram
ggplot(variog_res, aes(x = dist_km, y = gamma)) +
  geom_point(color = "steelblue") +
  geom_line(color = "steelblue") +
  labs(x = "Distance (km)", y = "Semivariance") +
  ggtitle("Semivariogram of Residuals") + 
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))




variog_fitted <- fit.variogram(variog_res,
                           model = vgm(psill = 100000, "Exp", range = 20000, nugget = 250000))

# Generate fitted model values
dist_seq <- seq(0, max(variog_res$dist), length.out = 200)
model_line <- variogramLine(variog_fitted, dist_vector = dist_seq)
model_line$dist_km <- model_line$dist / 1000

# Remove the zero-distance point
model_line <- model_line[model_line$dist_km > 0, ]

# Step 5: Plot using ggplot2
ggplot(variog_res, aes(x = dist_km, y = gamma)) +
  geom_point(color = "steelblue", size = 2) +
  geom_line(color = 'steelblue') + 
  geom_line(data = model_line, aes(x = dist_km, y = gamma), color = "firebrick", linewidth = 1) +
  labs(x = "Distance (km)", y = "Semivariance", title = "Semivariogram of Residuals with Fitted Model") +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))







# Create a grid directly in sf
sf_grid_proj <- st_make_grid(sf_proj, cellsize = 2000, what = "centers")  # 1000 meter resolution

# Convert to Spatial for kriging
grid_sp <- as(st_sf(geometry = sf_grid_proj), "Spatial")
gridded(grid_sp) <- TRUE

sp_proj <- as(sf_proj, "Spatial")

# Compute the median of mother_bmi
median_bmi <- median(birth$mother_bmi, na.rm = TRUE)


# Then use their most frequent levels and copy those to the grid
# Get most frequent level
most_common_sex <- names(which.max(table(birth$sex_of_child)))

# Assign it properly to grid_sp
grid_sp$sex_of_child <- factor(most_common_sex, levels = levels(birth$sex_of_child))
grid_sp$wealth_index <- 1
grid_sp$mother_bmi <- median_bmi

# Perform kriging (GP prediction)
kriged <- krige(
  formula = birth_weight ~ mother_bmi + sex_of_child + wealth_index,
  locations = sp_proj,
  newdata = grid_sp,
  model = variog_fitted
)

# Convert to sf and clean
kriged_sf <- st_as_sf(kriged) %>%
  filter(!is.na(var1.pred), is.finite(var1.pred)) %>%
  st_transform(st_crs(country))  # Align CRS

# Extract coordinates and drop geometry
kriged_df <- kriged_sf %>%
  mutate(x = st_coordinates(.)[, 1],
         y = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()


ggplot() +
  geom_tile(data = kriged_df, aes(x = x, y = y, fill = var1.pred), width = 0.05, height = 0.05) +
  geom_sf(data = country, fill = NA, color = "black") +
  scale_fill_viridis_c(option = "plasma", name = "Predicted\nBirth Weight") +
  coord_sf(expand = FALSE) +
  labs(
    title = "Universal Kriging Birth Weight Predictions",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5)
  )


#### Compiling results from all models
library(ggplot2)
library(dplyr)
library(patchwork)
library(forcats)


vecchia <- read_csv("vecchia_parameter_estimates.csv")

nn_beta <- read_csv("nngp_beta_samples.csv") %>% 
  select(-iter) %>%
  summarise(across(everything(), list(
    estimate = median,
    Lower_CI = ~quantile(.x, 0.025),
    Upper_CI = ~quantile(.x, 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("param", ".value"),
    names_pattern = "^(.*)_(estimate|Lower_CI|Upper_CI)$"
  )
  

nn_theta <- read_csv("nngp_theta_samples.csv") %>%
  select(-iter) %>%
  summarise(across(everything(), list(
    estimate = median,
    Lower_CI = ~quantile(.x, 0.025),
    Upper_CI = ~quantile(.x, 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("param", ".value"),
    names_pattern = "^(.*)_(estimate|Lower_CI|Upper_CI)$"
  )

mesh_beta <- read_csv('beta_1_mgp.csv')%>%
  select(-iter) %>%
  summarise(across(everything(), list(
    estimate = median,
    Lower_CI = ~quantile(.x, 0.025),
    Upper_CI = ~quantile(.x, 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("param", ".value"),
    names_pattern = "^(.*)_(estimate|Lower_CI|Upper_CI)$"
  )

mesh_tau <- read_csv('tausq_1_mgp.csv')%>%
  select(-iter) %>%
  summarise(across(everything(), list(
    estimate = median,
    Lower_CI = ~quantile(.x, 0.025),
    Upper_CI = ~quantile(.x, 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("param", ".value"),
    names_pattern = "^(.*)_(estimate|Lower_CI|Upper_CI)$"
  )

mesh_theta <- read_csv('theta_1_mgp.csv')%>%
  select(-iter) %>%
  summarise(across(everything(), list(
    estimate = median,
    Lower_CI = ~quantile(.x, 0.025),
    Upper_CI = ~quantile(.x, 0.975)
  ), .names = "{.col}_{.fn}")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("param", ".value"),
    names_pattern = "^(.*)_(estimate|Lower_CI|Upper_CI)$"
  )


lm_krig <- lm(birth_weight ~ mother_bmi + sex_of_child + wealth_index, sf_proj)
  
# Get coefficients
coefs <- summary(lm_krig)$coefficients

# Convert to tibble and compute CIs
lm_summary <- as_tibble(coefs, rownames = "param") %>%
  transmute(
    param,
    estimate = Estimate,
    Lower_CI = Estimate - 1.96 * `Std. Error`,
    Upper_CI = Estimate + 1.96 * `Std. Error`
  )


mesh_beta$param <- c("(Intercept)","mother_bmi","sex_of_childmale", "wealth_index")
nn_beta$param<- c("(Intercept)","mother_bmi","sex_of_childmale", "wealth_index")
mesh_beta[1:4, 2:4] <- mesh_beta[1:4, 2:4]*1000



# Add model identifiers
mesh_beta$model <- "Mesh"
nn_beta$model <- "NNGP"
vecchia_model <- vecchia[4:7, ] %>% rename(param = param_name)
vecchia_model$model <- "Vecchia"
lm_summary$model <- "Full GP"

# Combine all models into one dataframe
combined <- bind_rows(mesh_beta, nn_beta, vecchia_model, lm_summary)

combined <- combined %>%
  mutate(
    model = factor(model, levels = c("Full GP", "Vecchia", "NNGP", "Mesh"))
  )

# Add pretty math labels for each parameter
combined <- combined %>%
  mutate(param = recode(param,
                        "(Intercept)" = "beta[Intercept]",
                        "mother_bmi" = "beta[BMI]",
                        "sex_of_childmale" = "beta[Male]",
                        "wealth_index" = "beta[Wealth]"
  ))

combined <- combined %>%
  mutate(param_label = recode(param,
                              "beta[BMI]" = "beta[BMI]",
                              "beta[Intercept]" = "beta[Intercept]",
                              "beta[Male]" = "beta[Male]",
                              "beta[Wealth]" = "beta[Wealth]"
  ))
combined <- combined %>%
  mutate(
    model = factor(model, levels = c("Full GP", "Vecchia", "NNGP", "Mesh"))
  )

# Function to make plot for each parameter
make_plot <- function(param_expr) {
  combined %>%
    filter(param == param_expr) %>%
    ggplot(aes(x = estimate, y = model, color = model)) +
    geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI),
                    position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    labs(x = NULL, y = parse(text = param_expr)) +
    scale_y_discrete(limits = c("Full GP", "Vecchia", "NNGP", "Mesh")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      axis.title.y = element_text(
        angle = 0,
        hjust = 0.5,        # Center horizontally
        vjust = 0.5,        # Optional: tweak vertically
        size = 14,
        margin = margin(r = 10)  # Add space between label and axis
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# Create plots
plots <- c("beta[Intercept]","beta[BMI]", "beta[Male]", "beta[Wealth]") %>%
  lapply(make_plot)

# Combine with patchwork
wrap_plots(plots, ncol = 1, guides = "collect") &
  plot_annotation(title = "Coefficient Estimates by Model",
                  theme = theme(plot.title = element_text(hjust = 0.5))) &
  theme(legend.position = "right")




######## Covariance parameters 

vec <- vecchia[1:3, ] %>% rename(param = param_name)
vec$param <- c("sigmasq", 'phi', 'tausq')
nn_theta$param <- c('phi', 'sigmasq', 'tausq')
mesh <- rbind(mesh_tau, mesh_theta)
full_gp <- data_frame(param = c("tausq", "phi", "sigmasq"),
                      estimate = c(variog_fitted$psill[1], variog_fitted$range[2], variog_fitted$psill[2]),
                      Lower_CI = c(NA, NA, NA),
                      Upper_CI = c(NA, NA, NA))


# Add model labels and adjust as needed
full_gp <- full_gp %>%
  mutate(model = "Full GP",
         estimate = ifelse(param == "phi", 1/estimate, estimate))

vec <- vec %>%
  mutate(model = "Vecchia")

nn_theta <- nn_theta %>%
  mutate(model = "NNGP",
         estimate = ifelse(param == 'phi', estimate/100, estimate),
         Lower_CI = ifelse(param =='phi', Lower_CI/100, Lower_CI),
         Upper_CI = ifelse(param == 'phi', Upper_CI/100, Upper_CI))

mesh <- mesh %>%
  mutate(model = "Mesh",
         estimate = ifelse(param %in% c("tausq", "sigmasq"), estimate * 1e6, estimate),
         Lower_CI = ifelse(param %in% c("tausq", "sigmasq"), Lower_CI * 1e6, Lower_CI),
         Upper_CI = ifelse(param %in% c("tausq", "sigmasq"), Upper_CI * 1e6, Upper_CI))

# Combine all into one tidy dataframe
cov_params <- bind_rows(full_gp, vec, nn_theta, mesh) %>%
  mutate(model = factor(model, levels = c("Full GP", "Vecchia", "NNGP", "Mesh")))


make_cov_plot <- function(param_expr, nice_label) {
  cov_params %>%
    filter(param == param_expr) %>%
    ggplot(aes(x = estimate, y = model, color = model)) +
    geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI),
                    position = position_dodge(width = 0.6)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(x = NULL, y = nice_label) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title.y = element_text(angle = 0, hjust = 0.5, margin = margin(r = 10)),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    )
}

# Create each plot
tau_plot <- make_cov_plot("tausq", expression(tau^2))
sigma_plot <- make_cov_plot("sigmasq", expression(sigma^2))
phi_plot <- make_cov_plot("phi", expression(phi))


# Create individual plots
plots <- list(
  make_cov_plot("tausq", expression(tau^2)),
  make_cov_plot("sigmasq", expression(sigma^2)),
  make_cov_plot("phi", expression(phi))
)

# Combine and share legend
wrap_plots(plots, ncol = 1, guides = "collect") &
  plot_annotation(
    title = "Covariance Parameter Estimates",
    theme = theme(plot.title = element_text(hjust = 0.5))
  ) &
  theme(legend.position = "right")
