# Dynamic B0 from a numbers-at-age start: b0dyn() (icesdata package) with its doubling bug fixed. b0dyn() sets the starting numbers to stock.n + stock.n %*% srvAdj, i.e. the observed numbers plus the observed numbers with fishing removed, which counts them twice; here they are stock.n %*% srvAdj only. Used as the numbers-at-age alternative to b0dyn_eq_rec() (see 99-b0dyn-selftest.qmd and 07-sensitivity-b0.R)
b0dyn_cor = function(object, sr, residuals = 1) {
  # which years to set: the first year if youngest age is 0 or 1, the first 2 years if it is 2, etc.
  yrs = seq(max(1, dims(object)$min))
  # survival lost to fishing: at the start of the year a fish of age a has been fished at ages 1 to a-1, so age 1 gets 1, age 2 exp(F1), age 3 exp(F1 + F2), etc., using the starting year's F
  srvAdj = exp(apply(apply(harvest(object)[, yrs], c(1, 3:6), mean), 3:6, cumsum))
  srvAdj[-1] = srvAdj[-dim(srvAdj)[1]]
  srvAdj[1] = 1
  # bug fix: observed numbers with fishing removed (b0dyn() has stock.n(object)[, yrs] + stock.n(object)[, yrs] %*% srvAdj)
  stock.n(object)[, yrs] = stock.n(object)[, yrs] %*% srvAdj
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
