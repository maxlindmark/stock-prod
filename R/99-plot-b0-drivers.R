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
library(ggtext)
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

  # ML: whether each parameter is actually observed to vary by year for this
  # stock, so we want to calculate the average relative contribution not across all,
  # but by group (groups have same input time varying)
  wt_tv  <- !isTRUE(all.equal(as.numeric(stock.wt(x)), as.numeric(stock.wt(x_fix_wt))))
  mat_tv <- !isTRUE(all.equal(as.numeric(mat(x)), as.numeric(mat(x_fix_mat))))
  M_tv   <- !isTRUE(all.equal(as.numeric(m(x)), as.numeric(m(x_fix_M))))
  tv_labs <- c("Wt", "Mat", "M")[c(wt_tv, mat_tv, M_tv)]
  tv_group <- if (length(tv_labs) == 0) "None t.v." else paste(tv_labs, collapse = " + ")

  data.frame(sid = id, driver = names(V), var = V, var_pct = V / sum(V) * 100, tv_group = tv_group)
}, .id = NULL)

dat$driver <- factor(dat$driver,
  levels = c("wt", "M", "mat", "resid"),
  labels = c("Weight-at-age", "Natural mortality", "Maturity", "SR residuals")
)

driver_cols <- setNames(
  RColorBrewer::brewer.pal(nlevels(dat$driver), "Set2"),
  levels(dat$driver)
)

# ML: order tv_group facets by how many parameters are time-varying (fewest first)
tv_group_order <- dat |>
  distinct(sid, tv_group) |>
  mutate(n_tv = ifelse(tv_group == "None t.v.", 0, lengths(gregexpr("\\+", tv_group)) + 1)) |>
  distinct(tv_group, n_tv) |>
  arrange(n_tv, tv_group) |>
  pull(tv_group)
dat$tv_group <- factor(dat$tv_group, levels = tv_group_order)
resid_order <- dat$sid[dat$driver == "SR residuals"][order(dat$var_pct[dat$driver == "SR residuals"])]
dat$sid <- factor(dat$sid, levels = resid_order)

total_order <- dat |>
  summarise(total = sum(var), .by = sid) |>
  arrange(total) |>
  pull(sid)

dat$sid <- factor(dat$sid, levels = total_order)

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

# ML: color each tv_group abbreviation (Wt / Mat / M) to match the fill color
# of the driver it stands for, so the facet strip text doubles as a legend.
driver_cols_dark <- setNames(
  RColorBrewer::brewer.pal(nlevels(dat$driver), "Dark2"),
  levels(dat$driver)
)
abbrev_cols <- c(
  Wt  = driver_cols_dark[["Weight-at-age"]],
  Mat = driver_cols_dark[["Maturity"]],
  M   = driver_cols_dark[["Natural mortality"]]
)
colorize_tv_group <- function(x) {
  if (x == "None t.v.") {
    return(x)
  }
  parts <- strsplit(x, " \\+ ")[[1]]
  if (x == "Wt + Mat + M") {
    parts <- c("Wt", "M", "Mat")
  }
  paste(
    paste0("<span style='color:", abbrev_cols[parts], "'>", parts, "</span>"),
    collapse = " + "
  )
}

group_n <- dat |>
  distinct(sid, tv_group) |>
  count(tv_group)
group_labeller <- as_labeller(setNames(
  paste0(vapply(as.character(group_n$tv_group), colorize_tv_group, character(1)), " (n=", group_n$n, ")"),
  group_n$tv_group
))

p2 <- dat |>
  summarise(Mean = mean(var_pct), .by = c(driver, tv_group)) |>
  ggplot(aes(Mean, driver, fill = driver)) +
  geom_col() +
  facet_wrap(~factor(tv_group, levels = c("None t.v.",
                                          "Wt", 
                                          "Wt + M",
                                          "Wt + Mat",
                                          "Wt + Mat + M")),
             nrow = 1, labeller = group_labeller) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_discrete(limits = rev(levels(dat$driver))) +
  scale_x_continuous(breaks = c(25, 50, 75)) +
  theme(legend.position = "none", strip.text.x = ggtext::element_markdown()) +
  coord_cartesian(expand = 0) +
  labs(
    y = NULL,
    x = expression(
      paste(
        "Mean % of ",
        B[0],
        " variance"
      )
    ),
    fill = NULL
  )

p2

(p1 / p2) +
  plot_layout(heights = c(1, 0.1)) +
  plot_annotation(tag_levels = "A")

ggsave(paste0(home, "/figures/b0_drivers_mean.pdf"),
  width = 21, height = 24, units = "cm"
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
