### 0. SETUP
# Install required packages
install.packages(c('sf', 'dplyr', 'readr', 'httr', 'janitor', 'ggplot2', 'rdrobust'))

# Load required packages
library(sf)
library(dplyr)
library(readr)
library(httr)
library(janitor)
library(ggplot2)
library(rdrobust)

### 1. DOWNLOAD SHAPEFILES

## 1.1. Download shapefiles for TIF districts

tif_url <- "https://data.cityofchicago.org/resource/eejr-xtfb.geojson"
tif <- st_read(tif_url)

## 1.2. Download shapefiles for Sanitary and Ship Canal 

waterways_url <- "https://data.cityofchicago.org/resource/knfe-65pw.geojson"
waterways <- st_read(waterways_url)
canal <- waterways %>%
    filter(name == "SANITARY AND SHIP CANAL")

## 1.3. Dissolve canal segments into one shapefile

canal_dissolved <- canal %>%
    st_union() %>%
    st_as_sf()

## 1.4. Create a buffer to define the canal corridor
canal_buffer <- canal_dissolved %>%
 
  # Convert from lat/long to projected CRS for Chicago area: NAD83 UTM Zone 16N (EPSG 26916)
  st_transform(26916) %>%
  # Build a 1 kilometer buffer around the canal
  st_buffer(1000) %>%
  # Convert back to the orignal CRS
  st_transform(st_crs(canal_dissolved))

## 1.5. Drop invalid TIF geometries
tif <- tif[st_is_valid(tif),]

## 1.6. Create dummy variable for TIFs that touch corridor
tif <- tif %>%
  mutate(
    canal_tif = as.integer(lengths(st_intersects(., canal_buffer)) > 0)
  )

### 2. DOWNLOAD OUTCOME VARIABLE
## 2.1. Load parcel data from Cook County Assessor's Office
valuation <- read.csv("https://github.com/neildpatel/canaldistricts/blob/main/Assessor_-_Assessed_Values_20251119.csv")
parceluni <- read.csv("https://github.com/neildpatel/canaldistricts/blob/main/Assessor_-_Parcel_Universe_(Current_Year_Only)_20251119.csv")

# Join the CSVs using parcel pin-ID
parcelmerged <- inner_join(valuation, parceluni, by = "pin")

# A little data cleaning to convert the valuation field from string to numeric
parcelmerged$mailed_tot <- parcelmerged$mailed_tot %>% 
  gsub("\\$", "", x = _) %>%   
  gsub(",", "", x = _) %>%     
  as.numeric()

# A little more data cleaning - filter for parcel values over $10,000 and with valid lat/long coordinates
parcelvalues <- parcelmerged %>%
  filter(mailed_tot > 10000, 
         !is.na(latitude), 
         !is.na(longitude),
  )

# 2.2 Convert dataframe to shapefile
parcel_sf <- st_as_sf(
  parcelvalues,
  coords = c("longitude", "latitude"),
  crs = st_crs(tif)   
)

# 2.3. Transform to projected CRS to use distances in meters
tif_utm    <- st_transform(tif, 26916)
parcel_utm <- st_transform(parcel_sf, 26916)

# 2.4. Calculate distance to nearest TIF border for all TIFs
tif_union    <- st_union(tif_utm)
tif_boundary <- st_boundary(tif_union)

dist_mat <- st_distance(parcel_utm, tif_boundary)  
parcel_sf$dist_to_tif_border_m <- as.numeric(dist_mat)

# 2.5. Create a signed distance (inside TIF = positive, outside = negative)
inside_any_tif <- lengths(st_within(parcel_utm, tif_utm)) > 0

parcel_sf$signed_dist_to_tif_border_m <- ifelse(
  inside_any_tif,
  parcel_sf$dist_to_tif_border_m,
  -parcel_sf$dist_to_tif_border_m
)

# 2.6. Nearest TIF type (sanitary vs other) for each parcel
nearest_tif_idx <- st_nearest_feature(parcel_utm, tif_utm)

parcel_sf$nearest_canal_tif <- tif_utm$canal_tif[nearest_tif_idx]

parcel_sf$nearest_tif_type <- ifelse(
  parcel_sf$nearest_canal_tif == 1,
  "Sanitary Canal TIFs",
  "Non-Canal TIFs"
)

df <- parcel_sf %>%
  st_drop_geometry() %>%
  filter(
    !is.na(signed_dist_to_tif_border_m),
    !is.na(certified_tot),
    !is.na(nearest_tif_type),
    certified_tot > 0
  ) %>%
  mutate(
    log_certified = log(certified_tot)
  )

df_canal <- df %>% filter(nearest_tif_type == "Sanitary Canal TIFs")
df_other <- df %>% filter(nearest_tif_type == "Non-Canal TIFs")

nrow(df_canal)
nrow(df_other)

### 3. RDD ANALYSIS

## 3.1. Sanitary canal TIF borders
rd_canal <- rdrobust(
  y = df_canal$log_certified,
  x = df_canal$signed_dist_to_tif_border_m,
  c = 0
)

summary(rd_canal)

## 3.2. Other TIF borders
rd_other <- rdrobust(
  y = df_other$log_certified,
  x = df_other$signed_dist_to_tif_border_m,
  c = 0
)

summary(rd_other)

### 4. PLOTS

## 4.1 RDD Plot

# Extract optimal bandwidths from rdrobust objects
h_canal <- rd_canal$bws[1]   # est. bandwidth for canal TIF borders
h_other <- rd_other$bws[1]   # est. bandwidth for other TIF borders

h_canal
h_other

# Subset data to within each bandwidth
canal_band <- df_canal %>%
  filter(abs(signed_dist_to_tif_border_m) <= h_canal)

other_band <- df_other %>%
  filter(abs(signed_dist_to_tif_border_m) <= h_other)

# Split left/right and build triangular kernel weights

# Canal: left and right
canal_L <- canal_band %>%
  filter(signed_dist_to_tif_border_m < 0)
canal_R <- canal_band %>%
  filter(signed_dist_to_tif_border_m >= 0)

w_canal_L <- (h_canal - abs(canal_L$signed_dist_to_tif_border_m)) / h_canal
w_canal_R <- (h_canal - abs(canal_R$signed_dist_to_tif_border_m)) / h_canal

# Other TIFs: left and right
other_L <- other_band %>%
  filter(signed_dist_to_tif_border_m < 0)
other_R <- other_band %>%
  filter(signed_dist_to_tif_border_m >= 0)

w_other_L <- (h_other - abs(other_L$signed_dist_to_tif_border_m)) / h_other
w_other_R <- (h_other - abs(other_R$signed_dist_to_tif_border_m)) / h_other

# Local linear regressions (same structure as rdrobust)
fit_L_canal <- lm(log_certified ~ signed_dist_to_tif_border_m,
                  data = canal_L, weights = w_canal_L)
fit_R_canal <- lm(log_certified ~ signed_dist_to_tif_border_m,
                  data = canal_R, weights = w_canal_R)

fit_L_other <- lm(log_certified ~ signed_dist_to_tif_border_m,
                  data = other_L, weights = w_other_L)
fit_R_other <- lm(log_certified ~ signed_dist_to_tif_border_m,
                  data = other_R, weights = w_other_R)

# Quick check: the jump from these fits should match rdrobust
coef(fit_R_canal)[1] - coef(fit_L_canal)[1]
coef(fit_R_other)[1] - coef(fit_L_other)[1]

# Sequences of x-values within each bandwidth
xL_canal  <- seq(-h_canal, 0, length.out = 200)
xR_canal  <- seq(0, h_canal, length.out = 200)

xL_other  <- seq(-h_other, 0, length.out = 200)
xR_other  <- seq(0, h_other, length.out = 200)

# Prediction data frames

pred_canal_L <- data.frame(
  signed_dist_to_tif_border_m = xL_canal
) %>%
  mutate(
    log_fit = predict(fit_L_canal, newdata = .),
    tif_type = "Canal TIF boundary",
    side = "Left"
  )

pred_canal_R <- data.frame(
  signed_dist_to_tif_border_m = xR_canal
) %>%
  mutate(
    log_fit = predict(fit_R_canal, newdata = .),
    tif_type = "Canal TIF boundary",
    side = "Right"
  )

pred_other_L <- data.frame(
  signed_dist_to_tif_border_m = xL_other
) %>%
  mutate(
    log_fit = predict(fit_L_other, newdata = .),
    tif_type = "Non-Canal TIF boundary",
    side = "Left"
  )

pred_other_R <- data.frame(
  signed_dist_to_tif_border_m = xR_other
) %>%
  mutate(
    log_fit = predict(fit_R_other, newdata = .),
    tif_type = "Non-Canal TIF boundary",
    side = "Right"
  )

# Combine all fitted lines
pred_lines <- bind_rows(
  pred_canal_L, pred_canal_R,
  pred_other_L, pred_other_R
)

# Combine band data with a tif_type label
canal_band_plot <- canal_band %>%
  mutate(tif_type = "Canal TIF boundary")

other_band_plot <- other_band %>%
  mutate(tif_type = "Non-Canal TIF boundary")

plot_band <- bind_rows(canal_band_plot, other_band_plot)

# Set axis limits
x_min <- -100
x_max <- 100

y_min <- 9
y_max <- 17

ggplot() +
  # Gap at cutoff
  annotate("rect",
           xmin = -2, xmax = 2,
           ymin = -Inf, ymax = Inf,
           fill = "white", alpha = 0.5) +
  
  geom_point(
    data = plot_band %>% 
      filter(between(signed_dist_to_tif_border_m, -100, 100)),
    aes(x = signed_dist_to_tif_border_m, y = log_certified),
    alpha = 0.05,
    size = 0.6,
    color = "darkturquoise"
  ) +
  
  geom_line(
    data = pred_lines %>% 
      filter(between(signed_dist_to_tif_border_m, -100, 100)),
    aes(x = signed_dist_to_tif_border_m, y = log_fit, color = side),
    size = 1.6
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  coord_cartesian(
    xlim = c(-100, 100),
    ylim = c(9, 17.0)   # <<--- zoomed y-axis
  ) +
  
  labs(
    x = "Distance to nearest TIF boundary (m)\nnegative = outside TIF, positive = inside TIF, ",
    y = "log(parcel value)",
    title = "RDD results: Canal and Non-Canal TIFs",
    color = "Side"
  ) +
  facet_wrap(~ tif_type) +
  scale_color_manual(values = c("Left" = "blue", "Right" = "red")) +
  theme_minimal(base_size = 14)

## 4.2 Boxplot comparison of log parcel values for canal and non-canal TIFs

ggplot(df, aes(x = nearest_tif_type, y = log_certified)) +
  geom_boxplot(width = 0.5, fill = "turquoise", color = "grey20") +
  labs(
    title = "Log Parcel Values for Canal and Non-Canal ",
    x = "",
    y = "log(Parcel Value)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )

## 4.3 Visualizing the location and treatment assignment of parcels
parcel_plot_df <- parcel_sf %>%
  mutate(
    tif_type = case_when(
      nearest_canal_tif == 1 & inside_any_tif == 1 ~ "Inside Canal TIF",
      nearest_canal_tif == 1 & inside_any_tif == 0 ~ "Outside Canal TIF",
      nearest_canal_tif == 0 & inside_any_tif == 1 ~ "Inside Non-Canal TIF",
      nearest_canal_tif == 0 & inside_any_tif == 0 ~ "Outside Non-Canal TIF"
    )
  )

fill_colors <- c(
  "Inside Canal TIF"        = "#b30000",  # dark red
  "Outside Canal TIF"       = "#fc9272",  # light red
  "Inside Non-Canal TIF"    = "#08519c",  # dark blue
  "Outside Non-Canal TIF"   = "#9ecae1"   # light blue
)

tif_shaded <- tif %>%
  mutate(
    tif_category = case_when(
      canal_tif == 1 ~ "Inside Canal TIF",        
      canal_tif == 0 ~ "Inside Non-Canal TIF"
    )
  )

# Latitude bounds
lat_range <- range(st_coordinates(parcel_sf)[, "Y"], na.rm = TRUE)

# Plot with zoom applied
ggplot() +
  
# TIF polygons
  geom_sf(
    data = tif_shaded,
    aes(fill = tif_category),
    alpha = 0.25,
    color = "white",
    size = 0.3
  ) +
  
# Parcel points
  geom_sf(
    data = parcel_plot_df,
    aes(color = tif_type),
    alpha = 0.6,
    size = 0.4
  ) +
  
  scale_fill_manual(
    values = fill_colors,
    name = "TIF Type"
  ) +
  scale_color_manual(
    values = fill_colors,
    name = "Parcel Type"
  ) +
  
  labs(
    title = "Parcels Inside and Outside of Canal vs Non-Canal TIF Districts",
    x = "",
    y = ""
  ) +
  
  coord_sf(
    ylim = lat_range,
    expand = FALSE
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16)
  )

## 4.4. Single TIF

single_tif <- tif_shaded %>% 
  filter(tif_category == "Inside Canal TIF") %>%
  slice(3)

single_tif <- st_transform(single_tif, st_crs(parcel_plot_df))

single_boundary <- st_boundary(single_tif)

# Distance of all parcels to that boundary
d_to_boundary <- st_distance(parcel_plot_df, single_boundary)
d_to_boundary <- as.numeric(d_to_boundary)   # drop units

# Select parcels within 100m of the boundary
parcel_nearby <- parcel_plot_df %>%
  mutate(
    dist_to_single = d_to_boundary,
    inside_single  = st_intersects(geometry, single_tif, sparse = FALSE)
  ) %>%
  filter(dist_to_single <= 100)

single_ring <- st_buffer(single_boundary, dist = 100)

ggplot() +
# TIF polygon
  geom_sf(
    data = single_tif,
    fill = "grey80",
    color = "black",
    size = 0.6,
    alpha = 0.3
  ) +
  
# 100-meter ring
  geom_sf(
    data = single_ring,
    fill = "yellow",
    alpha = 0.15,
    color = NA
  ) +
  
# Parcels inside TIF
  geom_sf(
    data = parcel_nearby %>% filter(inside_single),
    color = "blue",     # dark red
    size = 1,
    alpha = 0.8
  ) +
  
# Parcels outside TIF
  geom_sf(
    data = parcel_nearby %>% filter(!inside_single),
    color = "turquoise",     # light red
    size = 1,
    alpha = 0.8
  ) +
  
# Zoom to the selected TIF
  coord_sf(
    xlim = st_bbox(single_tif)[c("xmin","xmax")],
    ylim = st_bbox(single_tif)[c("ymin","ymax")],
    expand = FALSE
  ) +
  
  labs(
    title = "Parcels within 100m of the Little Village East TIF boundary",
    x = "",
    y = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## 4.5 Canal corridor
# Subset canal TIFs
canal_tifs <- tif %>%
  filter(canal_tif == 1)

# Build a unified plotting object with identifiers for the legend
plot_canals <- canal %>%
  mutate(type = "Canal")

plot_buffer <- canal_buffer %>%
  mutate(type = "1 km Canal Buffer")

plot_tifs <- canal_tifs %>%
  mutate(type = "Canal TIF District")

# Combined bounding box for nice zooming
bbox <- st_bbox(canal_buffer)

# Colors for legend
fill_cols <- c(
  "1 km Canal Buffer" = "lightblue",
  "Canal TIF District" = "pink",
  "Canal (centerline)" = NA
)

line_cols <- c(
  "Canal" = "blue",
  "1 km Canal Buffer" = NA,
)

ggplot() +
  # BUFFER (polygon)
  geom_sf(
    data = plot_buffer,
    aes(fill = type),
    color = NA,
    alpha = 0.25
  ) +
  
  # TIF polygons
  geom_sf(
    data = plot_tifs,
    aes(fill = type, color = type),
    alpha = 0.5,
    linewidth = 0.5
  ) +
  
  # Canal line
  geom_sf(
    data = plot_canals,
    aes(color = type),
    linewidth = 1
  ) +
  
  scale_fill_manual(values = fill_cols, name = "Layers") +
  scale_color_manual(values = line_cols, name = "Layers") +
  
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  
  labs(
    title = "Sanitary & Ship Canal Buffer and Canal TIF Districts",
    x = "", y = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

