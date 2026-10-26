library(remotes)
library(FLCore)
library(icesdata)
library(FLasher)
library(FLSRTMB)
library(FLRef)
library(ggplotFL)
library(FLRebuild)
library(plyr)
library(dplyr)
library(reshape)

home <- here::here()

# dynamic B0 from a recruitment-adjusted unfished equilibrium start (see 08-b0-restart-test.R)
source(paste0(home, "/R/functions/b0dyn_eq_rec.R"))

# ML: icesdata is a named list of FLStock objects - one per ICES stock assessment.
# ML: Each FLStock holds the full assessment output: numbers-at-age, weight-at-age,
# ML: maturity, natural mortality (M), and fishing mortality (F) by age and year.
data(icesdata)

th <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.minor = element_blank())

# ML: Keep only stocks where the equilibrium simulation (eqsim) can estimate R0
# ML: (unfished recruitment). R0 is needed to set the stock-recruit model.
icesdata <- icesdata[laply(icesdata, function(x) "r0" %in% dimnames(eqsim(x))[[1]])]
icesdata <- icesdata[laply(icesdata, function(x) !is.na(eqsim(x)["r0"]))]

# add harvest and spawn.m
# ML: harvest.spwn and m.spwn control what fraction of annual mortality has already
# ML: occurred by the time fish spawn. Setting a floor of 1/12 avoids numerical
# ML: issues in the SPR (spawning potential ratio) calculations below.
nms <- names(icesdata)
icesdata <- mlply(names(icesdata), function(id) {
  print(id)
  x <- icesdata[[id]]
  harvest.spwn(x) <- qmax(harvest.spwn(x), 1 / 12)
  m.spwn(x) <- qmax(m.spwn(x), 1 / 12)
  x
})
names(icesdata) <- nms

# Assume r0 prior with CV = 0.3 (time series average of sel, m, mat and waa)
# ML: Fit a Beverton-Holt stock-recruit (SR) relationship per stock using TMB.
# ML: spr0y() computes SPR0 (expected SSB per recruit at zero fishing) for each year;
# ML: its time-series mean is used to scale the SR curve.
# ML: s.logitsd = 5 gives a weakly informative prior on steepness.
bhs <- (llply(names(icesdata), function(id) {
  print(id)
  x <- icesdata[[id]]
  r0 <- eqsim(x)["r0"]
  tryIt(srrTMB(as.FLSR(x, model = bevholt),
    spr0 = mean((spr0y(x))),
    r0.pr = c(r0, 0.3), s.logitsd = 5
  ))
}))
names(bhs) <- nms
bhs <- FLSRs(bhs[!laply(bhs, is.null)])

# ML: Save the filtered icesdata + fitted SRRs so 05-plot-b0-drivers.R can reuse
# ML: them without re-fitting ~80 Beverton-Holt models via TMB.
saveRDS(list(icesdata = icesdata, bhs = bhs), file = paste0(home, "/output/srr_fits.rds"))

# ML: b0dyn_eq_rec() projects the stock forward with F = 0 using the time-varying biological
# ML: parameters (wt, mat, M) and the fitted SR model. The result is the SSB the stock
# ML: would have had each year if it had never been fished ("dynamic B0").
# ML: It starts from the unfished equilibrium numbers-at-age (first-year biology, SR curve scaled by the mean recruitment deviation) instead of the observed numbers-at-age used by b0dyn() (icesdata package); this start was closest to the long-run B0 in the restart test (08-b0-restart-test.R).
# ML: B0.2 is the same projection but with SR residuals set to zero, giving the
# ML: deterministic (mean-curve) unfished trajectory.
dat <- mdply(names(bhs), function(id) {
  print(id)
  cbind(sid = id, as.data.frame(FLQuants(
    "SSB" = ssb(icesdata[[id]]),
    "B0" = ssb(b0dyn_eq_rec(icesdata[[id]], bhs[[id]])),
    "B0.2" = ssb(b0dyn_eq_rec(icesdata[[id]], bhs[[id]], residuals = 0)),
    "Rec" = rec(icesdata[[id]]),
    "Residuals" = residuals(bhs[[id]])
  ), drop = TRUE))
})

dat <- cast(dat, year + sid ~ qname, value = "data")

# ML: Species is available directly as an attribute on each icesdata FLStock
# ML: object, same as the stock id (`name`) -- no need to depend on
# ML: proddatafin.RData for it.
species_lookup <- ldply(names(icesdata), function(id) {
  data.frame(sid = id, species = attributes(icesdata[[id]])$species)
})

# ML: Save the object 03-join-env-flr.R expects: named `B0`, with columns
# ML: sid, year, B0 (absolute dynamic B0, tonnes), and species. B0 is not scaled by a static reference level (e.g. virgin SSB) because the HGAM standardizes B0 within stock, which cancels any per-stock constant.
B0 <- transmute(dat, sid = sid, year = year, B0 = B0)
B0 <- merge(B0, species_lookup, by = "sid")
saveRDS(B0, file = paste0(home, "/output/b0.rds"))

# 3 stocks are dropped because eqsim() couldn't produce a usable virgin-recruitment (R0) estimate for them — missing entirely for anf.27.3a46, and NA for cod.27.24-32 and her.27.6a7bc — and R0 is required to set up the stock-recruit model.



