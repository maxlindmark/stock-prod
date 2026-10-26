# Dynamic B0 from an unfished equilibrium start, without the recruitment-deviation adjustment: same as b0dyn_eq_rec(), but the equilibrium recruitment comes from the fitted SR curve itself (a - b / spr0), not the curve scaled by the mean recruitment deviation. Used only for comparison in 08-b0-restart-test.R
b0dyn_eq = function(object, sr, residuals = 1) {
  yrs = seq(max(1, dims(object)$min))
  m1 = c(m(object)[, 1])
  w1 = c(stock.wt(object)[, 1])
  mat1 = c(mat(object)[, 1])
  ms1 = c(m.spwn(object)[, 1])
  na = length(m1)
  l = cumprod(c(1, exp(-m1[-na])))
  if (!is.na(range(object)[["plusgroup"]])) l[na] = l[na] / (1 - exp(-m1[na]))
  spr0 = sum(l * w1 * mat1 * exp(-m1 * ms1))
  # equilibrium recruitment from the unscaled SR curve
  r_eq = c(params(sr)["a"]) - c(params(sr)["b"]) / spr0
  for (y in yrs) stock.n(object)[, y] = r_eq * l
  catch.n(object)[, yrs] = 0
  landings.n(object)[, yrs] = 0
  discards.n(object)[, yrs] = 0
  harvest(object)[, yrs] = 0
  rsdl = window(exp(residuals(sr) * residuals), start = dims(object)$minyear)
  stk0 = FLCore::ffwd(object, fbar = fbar(object)[, -yrs] %=% 0, sr = sr, deviances = rsdl)
  stk0
}
