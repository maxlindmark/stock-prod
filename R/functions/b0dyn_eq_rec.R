# Dynamic B0 from a recruitment-adjusted unfished equilibrium start. Same as b0dyn() (icesdata package), except for the starting numbers-at-age: instead of the observed numbers with fishing removed, the first year(s) get the unfished equilibrium numbers-at-age for the first year's biology, with the SR curve scaled by the mean recruitment deviation. Tested against b0dyn() and other starts in 99-b0dyn-selftest.qmd
b0dyn_eq_rec = function(object, sr, residuals = 1) {
  # which years to set: the first year if youngest age is 0 or 1, the first 2 years if it is 2, etc.
  yrs = seq(max(1, dims(object)$min))
  # first-year biology by age: natural mortality, stock weight, maturity, and fraction of M before spawning
  m1 = c(m(object)[, 1])
  w1 = c(stock.wt(object)[, 1])
  mat1 = c(mat(object)[, 1])
  ms1 = c(m.spwn(object)[, 1])
  na = length(m1)
  # proportion of a recruit surviving to each age from natural mortality only; the plus group is arrivals / fraction dying each year
  l = cumprod(c(1, exp(-m1[-na])))
  if (!is.na(range(object)[["plusgroup"]])) l[na] = l[na] / (1 - exp(-m1[na]))
  # SSB per recruit of an unfished stock at equilibrium
  spr0 = sum(l * w1 * mat1 * exp(-m1 * ms1))
  # mean recruitment deviation; residuals = 0 sets it to 1, i.e. the plain fitted curve
  k = if (residuals == 1) mean(exp(c(residuals(sr))), na.rm = TRUE) else 1
  # unfished equilibrium recruitment from the Beverton-Holt curve with a scaled by k: solves R = k * a * SSB / (b + SSB) with SSB = R * spr0
  r_eq = k * c(params(sr)["a"]) - c(params(sr)["b"]) / spr0
  # starting numbers-at-age: equilibrium recruitment times survival to each age
  for (y in yrs) stock.n(object)[, y] = r_eq * l
  # no catches or fishing in the starting years
  catch.n(object)[, yrs] = 0
  landings.n(object)[, yrs] = 0
  discards.n(object)[, yrs] = 0
  harvest(object)[, yrs] = 0
  # project forward with F = 0, recruitment from the SR curve times the recruitment deviations (residuals = 0: none)
  rsdl = window(exp(residuals(sr) * residuals), start = dims(object)$minyear)
  stk0 = FLCore::ffwd(object, fbar = fbar(object)[, -yrs] %=% 0, sr = sr, deviances = rsdl)
  stk0
}
