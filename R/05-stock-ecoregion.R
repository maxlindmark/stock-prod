# Assign each stock to an ICES ecoregion, used to group stocks in the 1995-2019 trend figure in 06-glmm-hgam.qmd (and in the exploratory ecoregion HGAMs in scrap/).
# Each stock's ICES areas (data/stocks_per_area2025.csv, same key and patches as 04-vbge-stock.qmd and 03-join-env-flr.R) are overlaid on the ICES ecoregion polygons, and the stock is assigned to the ecoregion covering the largest share of its area. Stocks where no ecoregion covers at least 50% of the area are "Widely distributed".
# Small ecoregions are merged with neighbours so that every group has enough stocks to estimate a regional trend: Faroes, Greenland and Icelandic Waters -> "Iceland, Faroes & Greenland"; Oceanic Northeast Atlantic (only Rockall haddock, 6b) -> Celtic Seas.
# Shapefiles are from the condition-hgam project.

library(sf)
library(dplyr)
library(readr)
library(tidyr)

home <- here::here()
shp_home <- "/Users/maxlindmark/Dropbox/max-work/R/condition-hgam/data/shapefiles/"

sf_use_s2(FALSE)

# Equal-area projection (ETRS89-LAEA) so overlap areas are comparable; polygons simplified to 2 km to speed up the intersection
areas <- st_read(paste0(shp_home, "ICES_areas/ICES_Areas_20160601_cut_dense_3857.shp"), quiet = TRUE) |>
  st_transform(3035) |>
  st_make_valid() |>
  st_simplify(dTolerance = 2000) |>
  st_make_valid() |>
  select(Area_27)

eco <- st_read(paste0(shp_home, "ICES_ecoregions/ICES_ecoregions_20171207_erase_ESRI.shp"), quiet = TRUE) |>
  st_transform(3035) |>
  st_make_valid() |>
  st_simplify(dTolerance = 2000) |>
  st_make_valid() |>
  select(Ecoregion)

# Area (m2) of each ICES area polygon falling in each ecoregion
overlap <- suppressWarnings(st_intersection(areas, eco)) |>
  mutate(area = as.numeric(st_area(geometry))) |>
  st_drop_geometry()

stocks <- unique(read_csv(paste0(home, "/output/dat.csv"), show_col_types = FALSE)$stock)

key <- read_csv(paste0(home, "/data/stocks_per_area2025.csv"), show_col_types = FALSE) |>
  mutate(
    sub_div = sub("^27\\.", "", ICES_Area),
    StockKeyLabel = ifelse(StockKeyLabel == "cod.21.1.isc", "cod.21.1", StockKeyLabel),
    StockKeyLabel = ifelse(StockKeyLabel == "cod.27.1-2.coastN", "cod.27.1-2coastN", StockKeyLabel)
  ) |>
  mutate(sub_div = ifelse(StockKeyLabel %in% c("cod.21.1", "cod.27.1-2coastN"), "1.b", sub_div)) |>
  filter(StockKeyLabel %in% stocks) |>
  distinct(stock = StockKeyLabel, sub_div)

stopifnot(all(stocks %in% key$stock))

# Key areas are given at different levels (e.g. "4" or "4.a"), so match each to all shapefile polygons at that level or finer
area_ids <- unique(overlap$Area_27)
key_poly <- key |>
  reframe(Area_27 = area_ids[area_ids == sub_div | startsWith(area_ids, paste0(sub_div, "."))], .by = c(stock, sub_div)) |>
  distinct(stock, Area_27)

shares <- key_poly |>
  left_join(overlap, by = "Area_27", relationship = "many-to-many") |>
  summarise(area = sum(area), .by = c(stock, Ecoregion)) |>
  mutate(share = area / sum(area), .by = stock)

stock_eco <- shares |>
  slice_max(share, n = 1, by = stock, with_ties = FALSE) |>
  transmute(stock, ecoregion_raw = Ecoregion, share) |>
  mutate(
    ecoregion = case_when(
      # West Greenland cod is in NAFO 1 (no ICES ecoregion), and the area key lists East Greenland cod's NAFO 1F part as ICES area 1 (Barents Sea), so both are set by hand
      stock %in% c("cod.21.1", "cod.2127.1f14") ~ "Greenland",
      # Icelandic stocks in 5a + 14; area 14 (East Greenland) makes the Greenland Sea the largest by area
      stock %in% c("aru.27.5a14", "usk.27.5a14") ~ "Icelandic Waters",
      # areas 5, 6, 12 and 14, only 51% in one ecoregion
      stock == "reg.27.561214" ~ "Widely distributed",
      share < 0.5 ~ "Widely distributed",
      TRUE ~ ecoregion_raw
    ),
    # merge small ecoregions with neighbours
    ecoregion = case_when(
      ecoregion %in% c("Faroes", "Greenland", "Icelandic Waters") ~ "Iceland, Faroes & Greenland",
      ecoregion == "Oceanic Northeast Atlantic" ~ "Celtic Seas",
      ecoregion == "Bay of Biscay and the Iberian Coast" ~ "Bay of Biscay & Iberian Coast",
      TRUE ~ ecoregion
    )
  ) |>
  arrange(ecoregion, stock)

stopifnot(setequal(stock_eco$stock, stocks))

count(stock_eco, ecoregion)

write_csv(stock_eco, paste0(home, "/output/stock_ecoregion.csv"))
