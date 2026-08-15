library(readr)
library(tidyr)
library(ggplot2)
library(dplyr)

home <- here::here()

d <- read_csv(paste0(home, "/output/sst_sal_sd_year_1940_2020.csv")) |>
  pivot_wider(names_from = "var", values_from = "value")

key <- read_csv(paste0(home, "/data/stocks_per_area2025.csv")) |>
  dplyr::mutate(sub_div = sub("^27\\.", "", ICES_Area)) |>
  dplyr::mutate(
    StockKeyLabel = ifelse(StockKeyLabel == "cod.21.1.isc", "cod.21.1", StockKeyLabel),
    StockKeyLabel = ifelse(StockKeyLabel == "cod.27.1-2.coastN", "cod.27.1-2coastN", StockKeyLabel)
  ) |>
  dplyr::mutate(sub_div = ifelse(StockKeyLabel %in% c("cod.21.1", "cod.27.1-2coastN"),
    "1.b", sub_div
  ))

# load(paste0(home, "/data/Final dataset productivity.RData"))

# b0.rds is the base stock-year grid: sid/year/B0/species, all sourced
# directly from icesdata in 02-b0.R. proddatafin.RData is no longer used at
# all -- abi, ffmsy, ssb_msy, yield_msy were dropped from the model.
B0 <- readRDS(paste0(home, "/output/b0.rds"))

b0 <- B0 |>
  dplyr::rename(
    .id = sid,
    b0 = B0,
    Species = species
  )

# With a list of keys for which stocks are in which ices sub div...
# Loop through each stock, find which sub divs, take mean, and end up with 3 col df:
# stock, sal, temp
# save for joining with FLR data

dl <- list()

for (i in unique(b0$.id)) {
  prod_sub <- b0 |>
    dplyr::filter(.id == i)

  key_sub <- key |> dplyr::filter(StockKeyLabel == i)

  sort(unique(d$sub_div))

  key_sub$sub_div

  var_sub <- d |>
    dplyr::filter(sub_div %in% key_sub$sub_div) |>
    dplyr::filter(year %in% prod_sub$year) |>
    dplyr::summarise(
      salinity = mean(sal),
      sst = mean(sst), .by = year
    )

  prod_sub <- prod_sub |> dplyr::left_join(var_sub, by = "year")

  dl[[i]] <- prod_sub
}

dat <- dplyr::bind_rows(dl)

# Now calculate mean temperature and sst in domain
domain_sub_divs <- b0 |>
  dplyr::distinct(.id) |>
  dplyr::left_join(key |> dplyr::rename(.id = StockKeyLabel), by = ".id")

var_annual_mean <- d |>
  dplyr::filter(sub_div %in% domain_sub_divs$sub_div) |>
  dplyr::summarise(
    annual_salinity = mean(sal),
    annual_sst = mean(sst), .by = year
  )

dat <- dat |>
  dplyr::left_join(var_annual_mean, by = "year")

dat <- dat |>
  dplyr::rename(stock = .id)

write_csv(dat, paste0(home, "/output/dat.csv"))

# Some plots
ggplot(dat, aes(year, sst, color = stock)) +
  geom_line() +
  theme(legend.position = "bottom") +
  guides(color = "none")

ggplot(dat, aes(year, salinity, color = stock)) +
  geom_line() +
  theme(legend.position = "bottom") +
  guides(color = "none")

dat |>
  dplyr::mutate(sst_sc = sst - mean(sst), .by = stock) |>
  ggplot(aes(year, sst_sc, color = stock)) +
  geom_line() +
  theme(legend.position = "bottom") +
  guides(color = "none")

dat |>
  dplyr::mutate(salinity_sc = salinity - mean(salinity), .by = stock) |>
  ggplot(aes(year, salinity_sc, color = stock)) +
  geom_line() +
  theme(legend.position = "bottom") +
  guides(color = "none")
