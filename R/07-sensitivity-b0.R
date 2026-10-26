library(FLCore)
library(FLasher)
library(FLSRTMB)
library(FLRef)
library(mgcv)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(ggsidekick)
library(RColorBrewer)

theme_set(theme_sleek())

home <- here::here()

# numbers-at-age start: b0dyn() (icesdata package) with its doubling bug fixed
source(paste0(home, "/R/functions/b0dyn_cor.R"))

# Sensitivity of the main results (04-glmm-hgam.qmd) to how the B0 projection is started:
# 1) Main run vs the same B0 with the first 10 years of each stock's projection removed, to drop any remaining start-up transient.
# 2) Main run (recruitment-adjusted equilibrium start, b0dyn_eq_rec()) vs a start from the observed numbers-at-age with fishing removed: b0dyn() (icesdata package) with its doubling bug fixed, b0dyn_cor().
# For each B0 version the HGAM from 04 is refitted and the global smooth and the per-stock percent change in B0 between 1995 and 2019 are compared with the main run.

# Main data (B0 from 02-b0.R, SST from 03-join-env-flr.R)
dat <- read_csv(paste0(home, "/output/dat.csv"))

# B0 from the numbers-at-age start (b0dyn_cor()), with the same filtered icesdata and SR fits as 02-b0.R
srr_fits <- readRDS(paste0(home, "/output/srr_fits.rds"))
icesdata <- srr_fits$icesdata
bhs <- srr_fits$bhs

b0_naa <- bind_rows(lapply(names(bhs), function(id) {
  tibble(stock = id, as.data.frame(ssb(b0dyn_cor(icesdata[[id]], bhs[[id]])), drop = TRUE)) |>
    transmute(stock, year, b0_naa = data)
}))

# The three B0 versions, all on the same stock-year grid and with the same covariates
dat_versions <- list(
  "Main" = dat,
  "First 10 years removed" = dat |> filter(year >= min(year) + 10, .by = stock),
  "Numbers-at-age start" = dat |> left_join(b0_naa, by = c("stock", "year")) |> mutate(b0 = b0_naa) |> select(-b0_naa)
)

# Data preparation and model, as in 04-glmm-hgam.qmd
mid_year <- round(median(1959:2022))

prep_data <- function(dat) {
  dat |>
    mutate(
      species_f = as.factor(Species),
      stock_f = as.factor(stock),
      year_ct = year - mid_year
    ) |>
    group_by(stock_f, species_f) |>
    arrange(year, .by_group = TRUE) |>
    mutate(sst_lag1 = lag(sst, 1)) |>
    ungroup() |>
    filter(year >= 1959) |>
    mutate(
      sst_sc_lag1 = (sst_lag1 - mean(sst_lag1, na.rm = TRUE)) / sd(sst_lag1, na.rm = TRUE),
      .by = stock_f
    ) |>
    mutate(
      b0_mean_stock = mean(b0, na.rm = TRUE),
      b0_sd_stock = sd(b0, na.rm = TRUE),
      b0_sc = (b0 - b0_mean_stock) / b0_sd_stock,
      .by = stock_f
    ) |>
    tidyr::drop_na(b0_sc, sst_sc_lag1) |>
    droplevels()
}

fit_hgam <- function(d) {
  mgcv::gam(
    b0_sc ~ s(year_ct, k = 12, m = 2) +
      s(year_ct, stock_f, bs = "fs", k = 4, m = 1),
    data = d,
    family = scat(),
    method = "REML"
  )
}

# Global smooth with 95% CI, as in 04
get_global <- function(m, d) {
  nd <- d |>
    distinct(year_ct) |>
    mutate(stock_f = levels(d$stock_f)[1], year = year_ct + mid_year)
  gl <- predict(m, newdata = nd, exclude = "s(year_ct,stock_f)", se.fit = TRUE)
  nd |> mutate(global = gl$fit, global_lwr = gl$fit - 1.96 * gl$se.fit, global_upr = gl$fit + 1.96 * gl$se.fit)
}

# Per-stock percent change in B0 between 1995 and 2019, simulated from the coefficient posterior and back-transformed to raw B0, as in 04
get_pct_change <- function(m, d, nsim = 2000) {
  b0_stats <- d |> distinct(stock_f, b0_mean = b0_mean_stock, b0_sd = b0_sd_stock)
  b <- coef(m)
  V <- vcov(m, unconditional = TRUE)
  set.seed(42)
  out <- list()
  for (s in levels(d$stock_f)) {
    yrs <- range(d$year[d$stock_f == s])
    # stock must cover both 1995 and 2019
    if (yrs[1] > 1995 || yrs[2] < 2019) next
    sub <- tibble(stock_f = factor(s, levels = levels(d$stock_f)), year = c(1995, 2019), year_ct = year - mid_year)
    st <- b0_stats |> filter(stock_f == s)
    Xp <- predict(m, newdata = sub, type = "lpmatrix")
    B <- mgcv::rmvn(nsim, b, V)
    early <- as.numeric(B %*% Xp[1, ]) * st$b0_sd + st$b0_mean
    late <- as.numeric(B %*% Xp[2, ]) * st$b0_sd + st$b0_mean
    pct <- (late / early - 1) * 100
    out[[s]] <- tibble(
      stock_f = s,
      pct_change = median(pct),
      pct_lwr = unname(quantile(pct, 0.025)),
      pct_upr = unname(quantile(pct, 0.975)),
      prop_early_nonpos = mean(early <= 0)
    )
  }
  # drop stocks with simulated 1995 B0 <= 0, as in 04
  bind_rows(out) |> filter(prop_early_nonpos == 0)
}

res <- lapply(names(dat_versions), function(v) {
  print(v)
  d <- prep_data(dat_versions[[v]])
  m <- fit_hgam(d)
  list(
    global = get_global(m, d) |> mutate(version = v),
    pct = get_pct_change(m, d) |> mutate(version = v),
    n_stocks = nlevels(d$stock_f),
    dev_expl = summary(m)$dev.expl
  )
})
names(res) <- names(dat_versions)

global <- bind_rows(lapply(res, `[[`, "global"))
pct <- bind_rows(lapply(res, `[[`, "pct")) |>
  mutate(sig = case_when(pct_lwr > 0 ~ "Increasing", pct_upr < 0 ~ "Decreasing", TRUE ~ "No change"))

# Summary numbers for the SI text
pct |>
  summarise(
    n_stocks = n(),
    median_pct = round(median(pct_change), 1),
    n_decreasing = sum(sig == "Decreasing"),
    n_increasing = sum(sig == "Increasing"),
    .by = version
  ) |>
  left_join(tibble(version = names(res), n_stocks_model = sapply(res, `[[`, "n_stocks"), dev_expl = round(sapply(res, `[[`, "dev_expl"), 3)), by = "version") |>
  as.data.frame()

# Agreement with the main run: per-stock percent change (Spearman correlation, stocks in both runs) and global smooth (correlation over shared years)
for (v in names(res)[-1]) {
  pp <- inner_join(filter(pct, version == "Main"), filter(pct, version == v), by = "stock_f")
  gg <- inner_join(filter(global, version == "Main"), filter(global, version == v), by = "year")
  cat(v, ": per-stock % change n =", nrow(pp), ", Spearman r =", round(cor(pp$pct_change.x, pp$pct_change.y, method = "spearman"), 2),
    ", same direction:", sum(sign(pp$pct_change.x) == sign(pp$pct_change.y)),
    ", same significance class:", sum(pp$sig.x == pp$sig.y),
    ", median % change main / alt (shared stocks):", round(median(pp$pct_change.x), 1), "/", round(median(pp$pct_change.y), 1),
    ", decreasing main / alt:", sum(pp$sig.x == "Decreasing"), "/", sum(pp$sig.y == "Decreasing"),
    ", increasing main / alt:", sum(pp$sig.x == "Increasing"), "/", sum(pp$sig.y == "Increasing"),
    ", Spearman (3 dec.):", round(cor(pp$pct_change.x, pp$pct_change.y, method = "spearman"), 3),
    "; global smooth r =", round(cor(gg$global.x, gg$global.y), 2), "\n")
}

# Plots: one figure per sensitivity test, A = global smooth, B = per-stock percent change vs the main run
pal <- brewer.pal(n = 3, "Dark2")

plot_sens <- function(v) {
  cols <- setNames(pal[1:2], c("Main", v))

  p_global <- global |>
    filter(version %in% c("Main", v)) |>
    mutate(version = factor(version, levels = c("Main", v))) |>
    ggplot(aes(year, global, color = version, fill = version)) +
    geom_ribbon(aes(ymin = global_lwr, ymax = global_upr), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    coord_cartesian(expand = 0) +
    guides(color = guide_legend(position = "inside"), fill = guide_legend(position = "inside")) +
    theme(legend.position.inside = c(0.05, 0.05), legend.justification = c(0, 0)) +
    labs(x = "Year", y = expression(paste("Global smooth of ", B[0]^{sc})), color = NULL, fill = NULL)

  pp <- inner_join(filter(pct, version == "Main"), filter(pct, version == v), by = "stock_f", suffix = c("_main", "_alt"))

  p_pct <- ggplot(pp, aes(pct_change_main, pct_change_alt)) +
    geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
    geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
    geom_abline(color = "grey40") +
    geom_errorbar(aes(ymin = pct_lwr_alt, ymax = pct_upr_alt), width = 0, alpha = 0.3, color = "grey30") +
    geom_errorbar(aes(xmin = pct_lwr_main, xmax = pct_upr_main), width = 0, alpha = 0.3, color = "grey30") +
    geom_point(size = 1.2, alpha = 0.8) +
    labs(
      x = expression(paste("Percent change in ", B[0], ", main")),
      y = bquote("Percent change in" ~ B[0] * "," ~ .(tolower(v)))
    )

  (p_global | p_pct) + plot_annotation(tag_levels = "A")
}

plot_sens("First 10 years removed")
ggsave(paste0(home, "/figures/supp/sens_b0_trim10.pdf"), width = 21, height = 9.5, units = "cm")

plot_sens("Numbers-at-age start")
ggsave(paste0(home, "/figures/supp/sens_b0_naa_start.pdf"), width = 21, height = 9.5, units = "cm")
