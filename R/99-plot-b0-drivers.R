library(FLCore)
library(icesdata)
library(FLasher)
library(FLSRTMB)
library(FLRef)
library(ggplotFL)
library(FLRebuild)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape)
library(ggsidekick)
theme_set(theme_sleek())

home <- here::here()

# ML: Reuse the filtered icesdata + fitted Beverton-Holt SRRs saved by 02-b0.R,
# ML: rather than re-fitting ~80 stock-recruit relationships via TMB again.
srr_fits <- readRDS(paste0(home, "/output/srr_fits.rds"))
icesdata <- srr_fits$icesdata
bhs <- srr_fits$bhs

th <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.minor = element_blank())

# ML: Fix an FLQuant slot at its age-specific time-series mean, removing year-to-year variation while preserving the age structure.
fix_to_mean <- function(flq) {
  means <- apply(flq, c(1, 3, 4, 5, 6), mean, na.rm = TRUE)
  for (y in seq(dim(flq)[2])) flq[, y, , , , ] <- means
  flq
}

# ML: Test the function
# grab a stock with time-varying weight-at-age
# x <- icesdata[["cod.27.1-2"]]
# # look at the raw weight-at-age - a matrix of age x year
# stock.wt(x)
# # fix it
# fix_to_mean(stock.wt(x))

# ML: For each stock, fix *one* parameter at its age-specific time-series mean and
# calculate how much B0 variance drops vs the baseline.
# If a WG assumed a parameter constant by year, that parameter's contribution is automatically zero.
# See example at the end...
dat <- ldply(names(bhs), function(id) {
  x <- icesdata[[id]]
  sr <- bhs[[id]]

  # ML: one copy per parameter, each with that parameter fixed at its mean
  x_fix_wt <- x
  stock.wt(x_fix_wt) <- fix_to_mean(stock.wt(x))
  x_fix_mat <- x
  mat(x_fix_mat) <- fix_to_mean(mat(x))
  x_fix_M <- x
  m(x_fix_M) <- fix_to_mean(m(x))

  # ML: b0dyn() projects the stock at F=0 using time-varying biology and the fitted SR;
  # ssb() extracts the resulting unfished SSB time series.
  # baseline has everything varying (that the stock has)
  # Each counterfactual scenario fixes one parameter.
  # The contribution of a parameter is the baseline variance minus variance when that
  # parameter is fixed.
  # SR residuals: same logic, but removed via residuals=0 argument in b0dyn.
  b0_base <- var(log(as.numeric(ssb(b0dyn(x, sr)))), na.rm = TRUE)

  V <- c(
    wt    = b0_base - var(log(as.numeric(ssb(b0dyn(x_fix_wt, sr)))), na.rm = TRUE),
    mat   = b0_base - var(log(as.numeric(ssb(b0dyn(x_fix_mat, sr)))), na.rm = TRUE),
    M     = b0_base - var(log(as.numeric(ssb(b0dyn(x_fix_M, sr)))), na.rm = TRUE),
    resid = b0_base - var(log(as.numeric(ssb(b0dyn(x, sr, residuals = 0)))), na.rm = TRUE)
  )

  V <- abs(V)

  if (is.na(sum(V)) || sum(V) == 0) {
    return(NULL)
  }

  data.frame(sid = id, driver = names(V), var = V, var_pct = V / sum(V) * 100)
}, .id = NULL)

dat$driver <- factor(dat$driver,
  levels = c("wt", "mat", "M", "resid"),
  labels = c("Weight-at-age", "Maturity", "Natural mortality", "SR residuals")
)
resid_order <- dat$sid[dat$driver == "SR residuals"][order(dat$var_pct[dat$driver == "SR residuals"])]
dat$sid <- factor(dat$sid, levels = resid_order)

total_order <- dat |>
  summarise(total = sum(var), .by = sid) |>
  arrange(total) |>
  pull(sid)

dat$sid <- factor(dat$sid, levels = total_order)

# dat_long <- dat |>
#   pivot_longer(c(var, var_pct), names_to = "type", values_to = "value") |>
#   mutate(type = ifelse(type == "var_pct", "Relative (% of total)", "Absolute (variance of log B0)"))
#
# p1 <- ggplot(dat_long, aes(value, sid, fill = driver)) +
#   geom_col() +
#   facet_wrap(~type, scales = "free_x") +
#   coord_cartesian(expand = 0) +
#   scale_fill_brewer(palette = "Set2") +
#   guides(fill = guide_legend(position = "inside")) +
#   theme(
#     legend.position.inside = c(0.35, 0.1),
#     axis.text.y = element_text(size = 8),
#     plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")
#     ) +
#   labs(x = NULL, y = NULL, fill = NULL)
#
# length(unique(dat$sid))

dat_long <- dat |>
  pivot_longer(c(var, var_pct), names_to = "type", values_to = "value") |>
  mutate(type = ifelse(type == "var_pct",
    "Relative~('%'~of~total)",
    "Absolute~(variance~of~log~B[0])"
  ))

p1 <- ggplot(dat_long, aes(value, sid, fill = driver)) +
  geom_col() +
  facet_wrap(~type, scales = "free_x", labeller = label_parsed) +
  coord_cartesian(expand = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  guides(fill = guide_legend(position = "inside")) +
  theme(
    legend.position.inside = c(0.35, 0.1),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 12, 5.5, 5.5, "pt")
  ) +
  labs(x = NULL, y = NULL, fill = NULL)

p1

# Plot mean and median % per driver across all stocks
p2 <- dat |>
  summarise(Mean = mean(var_pct), Median = median(var_pct), .by = driver) |>
  pivot_longer(c(Mean, Median), names_to = "stat", values_to = "value") |>
  filter(stat == "Mean") |>
  ggplot(aes(value, reorder(driver, value), fill = driver)) +
  geom_col() +
  facet_wrap(~stat, nrow = 2) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none") +
  coord_cartesian(expand = 0) +
  labs(
    y = NULL,
    x = expression(
      paste(
        "% of ",
        B[0],
        " variance"
      )
    ),
    fill = NULL
  )

(p1 / p2) +
  plot_layout(heights = c(1, 0.08)) +
  plot_annotation(tag_levels = "A")

ggsave(paste0(home, "/figures/b0_drivers_mean.pdf"),
  width = 19, height = 23, units = "cm"
)

# Plot some example b0 series to look at high vs low variance stocks over time
ids <- c("cod.27.47d20", "ple.27.21-23")

dat_ex <- ldply(ids, function(id) {
  cbind(sid = id, as.data.frame(ssb(b0dyn(icesdata[[id]], bhs[[id]])), drop = TRUE))
})

dat_ex <- dat_ex |>
  mutate(b0_sc = data / mean(data), .by = sid)

ggplot(dat_ex, aes(year, b0_sc)) +
  geom_line() +
  facet_wrap(~sid, scales = "fixed") +
  theme_minimal(10) +
  labs(x = NULL, y = "Dynamic B0 (relative to mean)")
