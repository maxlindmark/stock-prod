## Declining productivity of Northeast Atlantic fish stocks over six decades
In this paper we fit hierarchical generalized additive models to time series of dynamic unfished biomass ($B_0$) across 79 ICES-assessed stocks, to estimate each stock's trajectory relative to a shared temporal trend. Next we correlated each stock's 1995–2019 trend with its historical fishing pressure and life-history traits, to see if declines are steeper in more heavily fished or particular life-history stocks.

### Reproducing Results

To reproduce our results you can either:

1. Fork the repository, clone it, open a new RStudio project with version control, and paste the repo url

2. Download a zip and work locally on your computer

We use [`renv`](https://rstudio.github.io/renv/articles/renv.html) to manage package versions. Once you've downloaded the project, run `renv::restore()` in your current working directory. This will install the package versions we used when this repository was archived. Note that packages are installed in a stand-alone project library for this paper, and will not affect your installed R packages anywhere else! `renv` does *not* help with different versions of R. We used R version 4.6.1, and ran the analysis on a 24 GB Apple M2 Sequoia 15.6.1 laptop.

### Repository structure

`R`: code to prepare data, fit models, and make figures

`data`: stock-area csv; most data are generated from loaded R packages. `data/env` (CMIP6 SST and salinity from CNRM-CM6-1-HR, about 4 GB per file) is not included; download the files as described at the top of `R/01-sal-sst.R`. Only `01-sal-sst.R` needs them, and its output (`output/sst_sal_sd_year_1940_2020.csv`) is included, so you can skip step 01 and start from `02-b0.R`

`figures`: figures for paper are saved here

`output`: .rds files generated from FLR packages ready for analysis