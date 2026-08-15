library(dplyr)
library(readr)
library(tidyr)
library(slider)
library(purrr)
library(ggplot2)
library(ggrepel)
library(brms)
library(tidybayes)
library(icesdata)
library(FLCore)
library(rfishbase)
library(ggsidekick)
theme_set(theme_sleek())

home <- here::here()

# ML: Rebuild the same model data (d) and predictions used in 04-glmm-brms.qmd.
# Not saved as an artifact by that script, so reconstructed here from scratch.
dat <- read_csv(paste0(home, "/output/dat.csv"))

nao <- read_csv(paste0(home, "/data/env/NAO data - extended time series.csv")) |>
  dplyr::rename(year = Year)

nao <- nao |>
  arrange(year) |>
  mutate(nao = slide_dbl(winter_NAO, mean, .before = 4, .complete = TRUE))

mid_year <- round(median(1959:2022))

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
  filter(year >= 1959) |>
  left_join(nao, by = "year") |>
  mutate(
    sst_sc_lag1       = (sst_lag1 - mean(sst_lag1, na.rm = TRUE)) / sd(sst_lag1, na.rm = TRUE),
    salinity_sc_lag1  = (salinity_lag1 - mean(salinity_lag1, na.rm = TRUE)) / sd(salinity_lag1, na.rm = TRUE),
    .by = stock_f
  ) |>
  mutate(nao_sc = as.numeric(scale(nao))) |>
  mutate(
    b0_sc = (b0 - mean(b0, na.rm = TRUE)) / sd(b0, na.rm = TRUE),
    .by = stock_f
  ) |>
  tidyr::drop_na(b0_sc, sst_sc_lag1, salinity_sc_lag1, nao_sc) |>
  droplevels()

m2b <- readRDS(paste0(home, "/output/m2b.rds"))

# ML: 1996 anchor, not 1960: most stocks don't have predictions back to 1960
# (see 04-glmm-brms.qmd for the year-coverage check). Restrict to the 75 of 79
# stocks with data by 1996.
pm_smooth_with_re <- d |>
  tidyr::drop_na() |>
  filter(year >= 1959) |>
  add_epred_draws(m2b, ndraws = 5000, seed = 42)

stock_scale <- d |>
  summarise(
    mean_b0 = mean(b0, na.rm = TRUE),
    sd_b0 = sd(b0, na.rm = TRUE),
    .by = stock_f
  )

stock_min_year <- d |> summarise(min_year = min(year), .by = stock_f)
stocks_since_1996 <- stock_min_year |> filter(min_year <= 1996) |> pull(stock_f)

pct_data_1996 <- pm_smooth_with_re |>
  filter(stock_f %in% stocks_since_1996) |>
  left_join(stock_scale, by = "stock_f") |>
  mutate(b0_raw = .epred * sd_b0 + mean_b0) |>
  group_by(stock_f, .draw) |>
  summarise(
    b0_ref = b0_raw[which.min(abs(year - 1996))],
    b0_new = b0_raw[which.min(abs(year - 2020))],
    pct_change = (b0_new - b0_ref) / b0_ref * 100,
    .groups = "drop"
  ) |>
  group_by(stock_f) |>
  median_qi(pct_change, .width = 0.95) |>
  mutate(stock_f = as.character(stock_f))

# ML: Mean and cumulative F per stock, F-at-age averaged across ages each year,
# over each stock's full available history. Crude (not selectivity-weighted),
# fine for an exploratory pass, not a precise fishing-pressure metric.
f_summary <- names(icesdata) |>
  map(function(id) {
    stk <- icesdata[[id]]
    f_by_year <- apply(harvest(stk)@.Data, 2, mean, na.rm = TRUE)
    tibble(
      stock_f = id,
      mean_F = mean(f_by_year, na.rm = TRUE),
      cumulative_F = sum(f_by_year, na.rm = TRUE)
    )
  }) |>
  list_rbind()

# ML: Life-history traits + species, same source as the "Anything common..."
# section of 04-glmm-brms.qmd.
traits <- names(icesdata) |>
  map(function(stock_name) {
    stk <- icesdata[[stock_name]]
    fl <- attributes(stk)$fishlife
    if (is.null(fl)) {
      return(NULL)
    }
    tibble(
      stock_f = stock_name,
      species = attributes(stk)$species,
      linf = unname(fl["linf"]),
      k = unname(fl["k"]),
      age_maturity = mean(as.numeric(mat50(stk)), na.rm = TRUE)
    )
  }) |>
  list_rbind()

# ML: WoRMS species-to-family lookup from the stomachr package data.
load("/Users/maxlindmark/Dropbox/max-work/R/stomachr/data/worms_lookup.rda")
fam_lookup <- worms_lookup |> distinct(scientific_name, family)

# ML: A handful of stock species aren't in worms_lookup. Patched manually here
# (verified against standard fish taxonomy) rather than editing the stomachr
# package's data file, which belongs to a different project.
fam_manual <- tibble(
  scientific_name = c(
    "Lophius budegassa", "Dicentrarchus labrax", "Lepidorhombus boscii",
    "Sardina pilchardus", "Sebastes mentella", "Sebastes norvegicus",
    "Brosme brosme"
  ),
  family = c(
    "Lophiidae", "Moronidae", "Scophthalmidae",
    "Clupeidae", "Sebastidae", "Sebastidae",
    "Lotidae"
  )
)

fam_lookup <- bind_rows(fam_lookup, fam_manual)

pca_dat <- pct_data_1996 |>
  select(stock_f, pct_change) |>
  left_join(f_summary, by = "stock_f") |>
  left_join(traits, by = "stock_f") |>
  left_join(fam_lookup, by = c("species" = "scientific_name")) |>
  drop_na(mean_F, cumulative_F, linf, k, age_maturity)

cat("n stocks in analysis:", nrow(pca_dat), "\n")
cat("n stocks with family matched:", sum(!is.na(pca_dat$family)), "\n")
cat("Unmatched species (add to worms_lookup if family breakdown matters):\n")
print(unique(pca_dat$species[is.na(pca_dat$family)]))

# PCA on life-history / fishing-pressure traits
pca_fit <- pca_dat |>
  select(mean_F, cumulative_F, linf, k, age_maturity) |>
  scale() |>
  prcomp()

summary(pca_fit)
pca_fit$rotation

pca_scores <- as_tibble(pca_fit$x) |>
  bind_cols(pca_dat |> select(stock_f, species, pct_change)) |>
  mutate(direction = ifelse(pct_change < 0, "Decline", "Increase"))

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = direction)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = stock_f), size = 2, max.overlaps = 15) +
  scale_color_manual(values = c(Decline = "tomato2", Increase = "steelblue")) +
  labs(color = NULL)

p_pca

ggsave(paste0(home, "/figures/supp/pct_change_pca.png"), width = 20, height = 16, units = "cm")

# Direct correlations with pct_change
cor_mat <- pca_dat |>
  select(pct_change, mean_F, cumulative_F, linf, k, age_maturity) |>
  cor(use = "pairwise.complete.obs") |>
  round(2)

cor_mat

cor_tests <- map_dfr(c("mean_F", "cumulative_F", "linf", "k", "age_maturity"), function(v) {
  ct <- cor.test(pca_dat$pct_change, pca_dat[[v]])
  tibble(variable = v, r = unname(ct$estimate), p_value = ct$p.value)
})

cor_tests

# By family
p_family <- pca_dat |>
  filter(!is.na(family)) |>
  ggplot(aes(reorder(family, pct_change, median), pct_change)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  coord_flip() +
  labs(x = NULL, y = "% change in B0 (1996-2020)")

p_family

ggsave(paste0(home, "/figures/supp/pct_change_by_family.png"), width = 18, height = 14, units = "cm")

family_summary <- pca_dat |>
  filter(!is.na(family)) |>
  summarise(median_pct = median(pct_change), mean_pct = mean(pct_change), n = dplyr::n(), .by = family) |>
  arrange(median_pct)

family_summary

# ML: n per family is small (mostly 1-4 stocks; Gadidae n=26 is the exception),
# so treat these as suggestive, not confirmatory.
aov(pct_change ~ family, data = pca_dat |> filter(!is.na(family))) |> summary()
kruskal.test(pct_change ~ family, data = pca_dat |> filter(!is.na(family)))

## ---------------------------------------------------------------------------
## Cluster stocks on their actual response (not on predictor traits), then
## test whether cluster membership lines up with habitat guild, family, or F.
## ---------------------------------------------------------------------------

# ML: Total (fixed + random) SST and salinity slopes per stock, same source
# as the "Plot random slopes" section of 04-glmm-brms.qmd.
re_long <- m2b |>
  spread_draws(
    b_sst_sc_lag1, b_salinity_sc_lag1,
    r_stock_f[stock_f, term]
  ) |>
  filter(term %in% c("sst_sc_lag1", "salinity_sc_lag1")) |>
  ungroup() |>
  mutate(
    fixed_effect = case_when(
      term == "sst_sc_lag1" ~ b_sst_sc_lag1,
      term == "salinity_sc_lag1" ~ b_salinity_sc_lag1
    ),
    total_slope = fixed_effect + r_stock_f
  ) |>
  group_by(stock_f, term) |>
  median_qi(total_slope, .width = 0.95) |>
  ungroup() |>
  mutate(stock_f = as.character(stock_f)) |>
  select(stock_f, term, total_slope) |>
  pivot_wider(names_from = term, values_from = total_slope)

# ML: Demersal/pelagic/benthic habitat guild from FishBase. Hits a live web
# API on first call for any species not already locally cached by rfishbase,
# so this step needs network access and can be slow the first time it runs.
species_list <- unique(pca_dat$species)
fb <- species(species_list) |>
  as_tibble() |>
  select(species = Species, DemersPelag)

response_dat <- pct_data_1996 |>
  select(stock_f, pct_change) |>
  left_join(re_long, by = "stock_f") |>
  left_join(traits |> select(stock_f, species), by = "stock_f") |>
  left_join(fam_lookup, by = c("species" = "scientific_name")) |>
  left_join(fb, by = "species") |>
  left_join(f_summary, by = "stock_f") |>
  drop_na(pct_change, sst_sc_lag1, salinity_sc_lag1)

cat("\nn stocks with habitat guild matched:", sum(!is.na(response_dat$DemersPelag)), "of", nrow(response_dat), "\n")
cat("Unmatched species (FishBase name mismatch, check manually if needed):\n")
print(unique(response_dat$species[is.na(response_dat$DemersPelag)]))

# Hierarchical clustering on the response itself: how much a stock declined
# or increased, and how strongly it responds to SST and salinity.
clust_vars <- response_dat |>
  select(pct_change, sst_sc_lag1, salinity_sc_lag1) |>
  scale()

hc <- hclust(dist(clust_vars), method = "ward.D2")

k <- 3
response_dat$cluster <- factor(cutree(hc, k = k))

plot(hc, labels = response_dat$stock_f, main = paste0("Hierarchical clustering (k = ", k, ") of stock response"))
rect.hclust(hc, k = k)

# ML: Publication-style version, minimal ink, colorblind-safe palette
# (Dark2), 68% normal-density ellipses per cluster instead of point labels.
pal_cluster <- c("Cluster 1" = "#1b9e77", "Cluster 2" = "#d95f02", "Cluster 3" = "#7570b3")

p_cluster_nature <- response_dat |>
  mutate(cluster_lab = paste("Cluster", cluster)) |>
  ggplot(aes(sst_sc_lag1, pct_change, color = cluster_lab, fill = cluster_lab)) +
  stat_ellipse(geom = "polygon", alpha = 0.10, level = 0.68, linewidth = 0.4, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_color_manual(values = pal_cluster) +
  scale_fill_manual(values = pal_cluster) +
  labs(
    x = "Total SST slope (fixed + random)",
    y = "% change in B0 (1996 to 2020)",
    color = NULL, fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = c(0.84, 0.13),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 6, 4, 4)
  )

p_cluster_nature

ggsave(paste0(home, "/figures/supp/pct_change_clusters_nature.pdf"), width = 89, height = 80, units = "mm")
ggsave(paste0(home, "/figures/supp/pct_change_clusters_nature.png"), width = 89, height = 80, units = "mm", dpi = 600)

# Does cluster membership line up with habitat guild or family?
cat("\n=== Cluster x habitat guild ===\n")
print(table(response_dat$cluster, response_dat$DemersPelag, useNA = "ifany"))

cat("\n=== Cluster x family (top few) ===\n")
print(table(response_dat$cluster, response_dat$family, useNA = "ifany"))

# Direct test: does pct_change differ by habitat guild?
p_habitat <- response_dat |>
  filter(!is.na(DemersPelag)) |>
  ggplot(aes(reorder(DemersPelag, pct_change, median), pct_change)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  coord_flip() +
  labs(x = NULL, y = "% change in B0 (1996-2020)")

p_habitat

ggsave(paste0(home, "/figures/supp/pct_change_by_habitat.png"), width = 18, height = 10, units = "cm")

habitat_summary <- response_dat |>
  filter(!is.na(DemersPelag)) |>
  summarise(median_pct = median(pct_change), mean_pct = mean(pct_change), n = dplyr::n(), .by = DemersPelag) |>
  arrange(median_pct)

habitat_summary

cat("\n=== ANOVA / Kruskal-Wallis: pct_change ~ habitat guild ===\n")
aov(pct_change ~ DemersPelag, data = response_dat |> filter(!is.na(DemersPelag))) |> summary()
kruskal.test(pct_change ~ DemersPelag, data = response_dat |> filter(!is.na(DemersPelag)))

## ---------------------------------------------------------------------------
## What explains the decline/increase itself (pct_change), not the SST slope.
## Joint model: does habitat guild or fishing history add anything once each
## stock's own SST slope is accounted for?
## ---------------------------------------------------------------------------

response_dat2 <- response_dat |>
  mutate(bottom_group = ifelse(
    DemersPelag %in% c("benthopelagic", "bathydemersal", "demersal"),
    "Bottom-associated", "Pelagic"
  ))

fit_full <- lm(pct_change ~ mean_F + bottom_group + sst_sc_lag1, data = response_dat2)
cat("\n=== pct_change ~ mean_F + bottom_group + SST slope ===\n")
print(summary(fit_full))

r2_full <- summary(fit_full)$r.squared

p_explain <- response_dat2 |>
  ggplot(aes(sst_sc_lag1, pct_change, color = bottom_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey50") +
  geom_smooth(
    aes(group = 1), method = "lm", color = "grey30", fill = "grey70",
    linewidth = 0.6, alpha = 0.25
  ) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(values = c("Bottom-associated" = "#d95f02", "Pelagic" = "#1b9e77")) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = paste0("R² = ", round(r2_full, 2)),
    hjust = -0.15, vjust = 1.5, size = 3.2, color = "grey20"
  ) +
  labs(
    x = "Total SST slope (fixed + random)",
    y = "% change in B0 (1996 to 2020)",
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = c(0.82, 0.14),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.ticks = element_line(linewidth = 0.3, color = "black"),
    panel.grid = element_blank()
  )

p_explain

ggsave(paste0(home, "/figures/supp/pct_change_explained.pdf"), width = 12, height = 10, units = "cm")
ggsave(paste0(home, "/figures/supp/pct_change_explained.png"), width = 12, height = 10, units = "cm", dpi = 600)
