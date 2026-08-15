# 2025.11.11 Max Lindmark
# Script to read and prepare SST and SAL data

# 1. Reference for data product:

# Centre National de Recherches Météorologiques;
# Centre Européen de Recherche et Formation Avancée en Calcul Scientifique.
# WCRP CMIP6: the CNRM-CERFACS team CNRM-CM6-1-HR model output for the ‘historical’ experiment.
# Preprint at https://catalogue.ceda.ac.uk/uuid/f22d74c09e72480da1f56f29e0f80039 (2020).

# 2. Find the correct layer:
# https://catalogue.ceda.ac.uk/uuid/f22d74c09e72480da1f56f29e0f80039/

# Click download -> r1i1p1f2 (the core ensemble) -> Omon (Ocean monthly) -> tos (SST), sos (SAL)
# -> gn -> latest (takes you to a .nc file). These are around 4 gb!

# Final path is:
# https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1-HR/historical/r1i1p1f2/Omon
#
# ... And for forecasts:
# https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/CNRM-CERFACS/CNRM-CM6-1-HR/ssp245/r1i1p1f2/Omon/sos/gn/v20191202

# -----------------------------
# Load packages
library(tidyverse)
library(ncdf4)
library(lubridate)
library(tidylog)
library(tictoc)
library(data.table)
library(viridis)
library(raster)
library(ggsidekick)
theme_set(theme_sleek())

# -----------------------------
# Memory-safe function to read one time step at a time
read_nc_memory_safe_ncdf4 <- function(path, varname,
                                      lat_bounds = c(35, 80),
                                      lon_bounds = c(-30, 30),
                                      start_date = as.Date("1940-01-01"),
                                      end_date = NULL) {
  nc <- nc_open(path)

  lat <- ncvar_get(nc, "lat") # [y,x]
  lon <- ncvar_get(nc, "lon") # [y,x]
  time <- ncvar_get(nc, "time") # days since 1850-01-01
  origin <- as.Date("1850-01-01")
  time_dates <- origin + time
  if (is.null(end_date)) end_date <- max(time_dates)

  # Lat/lon mask
  mask <- lat >= lat_bounds[1] & lat <= lat_bounds[2] &
    lon >= lon_bounds[1] & lon <= lon_bounds[2]

  # Time indices in desired range
  idx_time <- which(time_dates >= start_date & time_dates <= end_date)
  df_list <- vector("list", length(idx_time))

  for (i in seq_along(idx_time)) {
    t_idx <- idx_time[i]
    message("Reading time step: ", i, "/", length(idx_time), " (", time_dates[t_idx], ")")

    # Read one time step [y,x] at t_idx
    slice <- ncvar_get(nc, varid = varname, start = c(1, 1, t_idx), count = c(-1, -1, 1))

    # Flatten to tidy
    yx_idx <- which(mask, arr.ind = TRUE)
    df <- tibble(
      y = yx_idx[, 1],
      x = yx_idx[, 2],
      latitude = lat[mask],
      longitude = lon[mask],
      date = time_dates[t_idx],
      year = year(date),
      month = month(date),
      value = slice[mask]
    )

    df_list[[i]] <- df
  }

  nc_close(nc)
  bind_rows(df_list)
}

# -----------------------------
# Paths to CMIP6 data
home <- here::here()

# Historical (1940-2014)
path_tos_hist <- paste0(home, "/data/env/tos_Omon_CNRM-CM6-1-HR_historical_r1i1p1f2_gn_185001-201412.nc")
path_sos_hist <- paste0(home, "/data/env/sos_Omon_CNRM-CM6-1-HR_historical_r1i1p1f2_gn_185001-201412.nc")

tos_hist <- read_nc_memory_safe_ncdf4(path_tos_hist, "tos", start_date = as.Date("1940-01-01"))
sos_hist <- read_nc_memory_safe_ncdf4(path_sos_hist, "sos", start_date = as.Date("1940-01-01"))

d_old <- bind_rows(
  tos_hist |> mutate(var = "sst"),
  sos_hist |> mutate(var = "sal")
)

# Forecast (2015-2025)
path_tos_f <- paste0(home, "/data/env/tos_Omon_CNRM-CM6-1-HR_ssp245_r1i1p1f2_gn_201501-210012.nc")
path_sos_f <- paste0(home, "/data/env/sos_Omon_CNRM-CM6-1-HR_ssp245_r1i1p1f2_gn_201501-210012.nc")

tos_f <- read_nc_memory_safe_ncdf4(path_tos_f, "tos", start_date = as.Date("2015-01-01"), end_date = as.Date("2025-12-31"))
sos_f <- read_nc_memory_safe_ncdf4(path_sos_f, "sos", start_date = as.Date("2015-01-01"), end_date = as.Date("2025-12-31"))

d_new <- bind_rows(
  tos_f |> mutate(var = "sst"),
  sos_f |> mutate(var = "sal")
)

# Average by month and combine historical + forecast
tic()
dt_old <- as.data.table(d_old)
dt_old_annual <- dt_old[, .(value = mean(value, na.rm = TRUE)),
  by = .(y, x, latitude, longitude, year, var)
]
d_old_sum <- as_tibble(dt_old_annual)
toc()

tic()
dt_new <- as.data.table(d_new)
dt_new_annual <- dt_new[, .(value = mean(value, na.rm = TRUE)),
  by = .(y, x, latitude, longitude, year, var)
]
d_new_sum <- as_tibble(dt_new_annual)
toc()

d <- bind_rows(d_old_sum, d_new_sum)

# -----------------------------
# Add ICES sub divisions

shape <- shapefile(paste0(home, "/data/shapefiles/StatRec_map_Areas_Full_20170124.shp"))

pts <- SpatialPoints(cbind(d$longitude, d$latitude),
  proj4string = CRS(proj4string(shape))
)

d$sub_div <- over(pts, shape)$Area_27

d <- d |> drop_na(sub_div)

# -----------------------------
# Make some plots

d <- d |> drop_na(value)

# d |>
#   filter(year == min(year)) |>
#   ggplot(aes(longitude, latitude, color = value)) +
#   facet_wrap(~var, scales = "free") +
#   geom_point(size = 0.3) +
#   scale_color_viridis() +
#   coord_cartesian(expand = 0) +
#   theme(legend.position = "bottom")
#
# ggsave(paste0(home, "/output/first_year_map.png"), width = 21, height = 11, unit = "cm")
#
# # Plot annual means by sub-div
# d |>
#   summarise(value = mean(value),
#             .by = c(year, var, sub_div)) |>
#   ggplot(aes(year, value, color = sub_div)) +
#   geom_line(alpha = 0.5) +
#   facet_wrap(~var, scales = "free") +
#   theme(legend.position = "bottom",
#         legend.title = element_blank()) +
#   guides(color = guide_legend(ncol = 12))
#
# ggsave(paste0(home, "/output/annual_means_by_sd.png"), width = 22, height = 16, unit = "cm")
#
# # Now do z-scored trends by sd
# d |>
#   summarise(value = mean(value),
#             .by = c(year, var, sub_div)) |>
#   mutate(value_sc = value - mean(value), .by = c(var, sub_div)) |>
#   ggplot(aes(year, value_sc, color = sub_div)) +
#   geom_line(alpha = 0.5) +
#   facet_wrap(~var, scales = "free") +
#   theme(legend.position = "bottom",
#         legend.title = element_blank()) +
#   guides(color = guide_legend(ncol = 12))
#
# ggsave(paste0(home, "/output/annual_means_by_sd_sc.png"), width = 22, height = 16, unit = "cm")


d <- d |>
  summarise(value = mean(value), .by = c(sub_div, year, var))

# Save annual data!
write_csv(d, paste0(home, "/output/sst_sal_sd_year_1940_2020.csv"))
