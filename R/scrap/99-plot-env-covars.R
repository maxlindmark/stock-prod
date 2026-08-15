library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(ggsidekick); theme_set(theme_sleek())

home <- here::here()

## Plot environmental data
#env_dat <- read_csv(paste0(home, "/output/sst_sal_sd_year_1940_2020.csv"))

mid_year <- round(median(1959:2022))

dat <- read_csv(paste0(home, "/output/dat.csv"))

d <- dat |>
  mutate(
    species_f = as.factor(Species),
    stock_f = as.factor(stock),
    year_ct = year - mid_year
  ) |>
  group_by(stock_f, species_f) |>
  arrange(year, .by_group = TRUE) |>
  mutate(
    sst_lag1      = lag(sst, 1),
    salinity_lag1 = lag(salinity, 1)
  ) |>
  ungroup() |>
  # Filter to analysis period BEFORE scaling
  filter(year >= 1959) |>
  left_join(nao, by = "year") |>
  # Scale within stock (mean and SD), matching how the response is scaled
  mutate(
    sst_sc_lag1       = (sst_lag1 - mean(sst_lag1, na.rm = TRUE)) / sd(sst_lag1, na.rm = TRUE),
    salinity_sc_lag1  = (salinity_lag1 - mean(salinity_lag1, na.rm = TRUE)) / sd(salinity_lag1, na.rm = TRUE),
    .by = stock_f
  ) |>
  # Scale response within stock
  mutate(
    b0_sc = (b0 - mean(b0, na.rm = TRUE)) / sd(b0, na.rm = TRUE),
    .by = stock_f
  ) |>
  tidylog::drop_na(b0_sc) |>
  tidylog::drop_na(sst_sc_lag1) |>
  tidylog::drop_na(salinity_sc_lag1) |>
  droplevels()

plotdat <- d |>
  ungroup() |>
  select(year, stock_f, sst_sc_lag1, salinity_sc_lag1) |>
  pivot_longer(c(sst_sc_lag1, salinity_sc_lag1), names_to = "var", values_to = "value_sc") |>
  mutate(var = factor(var,
                      levels = c("sst_sc_lag1", "salinity_sc_lag1"),
                      labels = c("SST", "Salinity")))

plotdat |>
  ggplot(aes(year, value_sc, color = stock_f)) +
  geom_line(alpha = 0.5) +
  geom_smooth(aes(year, value_sc), inherit.aes = FALSE,
              method = "gam", formula = y ~ s(x, k = 7), color = "grey30") +
  coord_cartesian(expand = 0) +
  scale_color_viridis_d(option = "mako") +
  facet_wrap(~var, nrow = 1, scales = "free_y") +
  geom_text(
    data = data.frame(
      var = factor(c("SST", "Salinity"), levels = c("SST", "Salinity")),
      label = c("A", "B")
    ),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.5, vjust = 1.5, size = 4, fontface = "plain", inherit.aes = FALSE
  ) +
  scale_x_continuous(breaks = seq(1965, 2015, by = 10)) +
  guides(color = "none") +
  labs(x = "Year", y = "z-scored (1-year lag, within-stock scaled)") +
  theme(legend.position = "none")

ggsave(paste0(home, "/figures/supp/sst_sal_trends.pdf"),
       width = 20, height = 11, units = "cm"
)
