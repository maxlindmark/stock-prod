library(FLCore)
library(FLasher)
library(FLSRTMB)
library(dplyr)
library(ggplot2)
library(ggsidekick)
library(RColorBrewer)

theme_set(theme_sleek())

home <- here::here()

# Restart test of the starting numbers-at-age of the B0 projection (supporting figures for the choice of start, see 99-b0dyn-selftest.qmd for the full exploration).
# Each stock's B0 is projected over its full series, and then restarted 30, 35, ..., 55 years into the series with the same start method. By then the full-length run no longer depends on its start, so the closer the restarted run starts to it, and the faster it converges on it, the smaller the transient a start causes.
# Three starts are compared: the observed numbers-at-age with fishing removed (b0dyn() with its doubling bug fixed), the unfished equilibrium, and the unfished equilibrium with the SR curve scaled by the mean recruitment deviation (main method).

source(paste0(home, "/R/functions/b0dyn_cor.R"))
source(paste0(home, "/R/functions/b0dyn_eq.R"))
source(paste0(home, "/R/functions/b0dyn_eq_rec.R"))

# filtered icesdata + fitted Beverton-Holt SRRs saved by 02-b0.R
srr_fits <- readRDS(paste0(home, "/output/srr_fits.rds"))
icesdata <- srr_fits$icesdata
bhs <- srr_fits$bhs

versions <- list(
  "Numbers-at-age" = b0dyn_cor,
  "Equilibrium" = b0dyn_eq,
  "Equilibrium, rec-dev adjusted" = b0dyn_eq_rec
)
version_cols <- setNames(brewer.pal(3, "Dark2"), names(versions))

offsets <- seq(30, 55, by = 5)

# log difference between restarted and full-length B0 in every year from the restart onwards; each method is compared with its own full-length run
st_traj <- list()

for (id in names(bhs)) {
  x <- icesdata[[id]]
  sr <- bhs[[id]]

  for (v in names(versions)) {
    ref <- ssb(versions[[v]](x, sr))

    for (off in offsets) {
      T0 <- dims(x)$minyear + off
      # only restarts with at least 10 years left to project
      if (dims(x)$maxyear - T0 < 9) next

      run <- ssb(versions[[v]](window(x, start = T0), sr))
      yrs <- T0:dims(x)$maxyear
      diff <- log(c(run[, as.character(yrs)])) - log(c(ref[, as.character(yrs)]))
      st_traj[[length(st_traj) + 1]] <- tibble(sid = id, version = v, offset = off, since = yrs - T0, diff = diff)
    }
  }
}

st_traj <- bind_rows(st_traj) |>
  mutate(version = factor(version, levels = names(versions)))

# Summary numbers for the SI text
# restart year: percent difference between restarted and full-length B0
st <- filter(st_traj, since == 0)
st |>
  summarise(
    median_pct = round(100 * (exp(median(diff)) - 1)),
    q25_pct = round(100 * (exp(quantile(diff, 0.25)) - 1)),
    q75_pct = round(100 * (exp(quantile(diff, 0.75)) - 1)),
    n = n(),
    n_stocks = n_distinct(sid),
    .by = version
  ) |>
  as.data.frame()

# years until within 10% (|log difference| < 0.1) for the rest of the series; NA if still off in the last year
conv_yrs <- st_traj |>
  summarise(last_off = max(c(-1, since[abs(diff) >= 0.1])), last_yr = max(since), .by = c(sid, version, offset)) |>
  mutate(yrs = ifelse(last_off == last_yr, NA, last_off + 1))

conv_yrs |>
  summarise(
    median_years = median(yrs, na.rm = TRUE),
    within_10pct_from_restart = round(100 * mean(!is.na(yrs) & yrs == 0)),
    converged_within_10yrs = round(100 * mean(!is.na(yrs) & yrs <= 10)),
    not_converged = round(100 * mean(is.na(yrs))),
    n = n(),
    .by = version
  ) |>
  as.data.frame()

# label a log ratio as the percent difference it corresponds to
pct_lab <- function(x) ifelse(abs(x) < 1e-9, "0%", sprintf("%+.0f%%", 100 * (exp(x) - 1)))

# Restart year: restarted vs full-length B0 by restart offset
n_stocks <- st |>
  distinct(sid, offset) |>
  count(offset)

ggplot(st, aes(factor(offset), diff, color = version)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.4, size = 0.8) +
  geom_boxplot(aes(fill = version), position = position_dodge(width = 0.8), width = 0.15, outlier.shape = NA, alpha = 0.3, linewidth = 0.4) +
  geom_text(data = n_stocks, aes(factor(offset), -Inf, label = paste0("n=", n)), inherit.aes = FALSE, vjust = -0.5, size = 3, color = "grey40") +
  scale_color_manual(values = version_cols) +
  scale_fill_manual(values = version_cols) +
  scale_y_continuous(breaks = log(c(0.05, 0.1, 0.25, 0.5, 1, 2, 4)), labels = pct_lab) +
  labs(x = "Restart, years into series", y = expression(paste("Restarted vs full-length ", B[0], " (log scale)")), color = NULL, fill = NULL) +
  theme(legend.position = "bottom")

ggsave(paste0(home, "/figures/supp/b0_restart_start.pdf"), width = 18, height = 12, units = "cm")

# Convergence after the restart: median and interquartile range of |log difference| by year since restart, pooled over stocks and offsets, up to 20 years (from year 10 the number of stock-offset combinations declines)
conv <- st_traj |>
  filter(since <= 20) |>
  summarise(med = median(abs(diff)), q25 = quantile(abs(diff), 0.25), q75 = quantile(abs(diff), 0.75), .by = c(version, since))

ggplot(conv, aes(since, med, color = version, fill = version)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  geom_line() +
  geom_point(size = 1) +
  scale_color_manual(values = version_cols) +
  scale_fill_manual(values = version_cols) +
  scale_y_continuous(breaks = log(1 + c(0, 0.1, 0.25, 0.5, 1, 2)), labels = \(x) sprintf("%.0f%%", 100 * (exp(x) - 1))) +
  labs(x = "Years since restart", y = "Absolute difference (log scale)", color = NULL, fill = NULL) +
  theme(legend.position = "bottom")

ggsave(paste0(home, "/figures/supp/b0_restart_convergence.pdf"), width = 14, height = 11, units = "cm")
