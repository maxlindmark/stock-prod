# Calculate ABI by stock and year, see 10.1111/faf.12789
library(readr)
library(remotes)
library(FLCore)
library(FLBRP)
library(FLRebuild)
library(doParallel)
library(foreach)
library(plyr)
library(dplyr)
library(icesdata)

# icesdata = icesdata
ids <- names(icesdata)

ncores <- detectCores() - 4
cl <- makeCluster(ncores)
registerDoParallel(cl)

## Estimating a Beverton and Holt for all stocks
bh1 <- foreach(id = ids, .packages = c("FLCore", "FLBRP", "icesdata")) %dopar% {
  landings.n(icesdata[[id]]) <- catch.n(icesdata[[id]])
  landings.wt(icesdata[[id]]) <- catch.wt(icesdata[[id]])
  discards.n(icesdata[[id]])[] <- 0
  discards.wt(icesdata[[id]])[] <- 0
  catch(icesdata[[id]]) <- computeCatch(icesdata[[id]], slot = "all")
  icesdata::eql(icesdata[[id]], model = "bevholtSV")
}

bh1 <- FLBRPs(bh1)
names(bh1) <- ids
stopCluster(cl)

rfs <- ldply(bh1, function(x) refpts(x)["msy", "ssb", drop = TRUE])
colnames(rfs)[2] <- "Bmsy"

ts <- mdply(names(icesdata), function(id) {
  cbind(
    .id = id,
    model.frame(FLQuants(ABImsy = abi(icesdata[[id]], bh1[[id]], ref = "msy")), drop = TRUE)
  )
})[, -c(1, 12:14, 19)]

ts_abi <- ts[c(1:3)]


write_csv(ts_abi, paste0(here::here(), "/output/ts_abi.csv"))
