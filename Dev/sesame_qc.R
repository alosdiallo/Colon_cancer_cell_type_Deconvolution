library(sesame)
sesameDataCache()
library(minfi)
library(ggplot2)

idat <- "idats"
sdfs = openSesame(idat, func = NULL) 
saveRDS(sdfs, "sesame_variables/sdfs.RDS")
betas = openSesame(idat, func = getBetas) 
saveRDS(betas, "sesame_variables/betas.RDS")
qcs = openSesame(idat, prep="QCDPB", func=sesameQC_calcStats)
saveRDS(qcs, "sesame_variables/qcs.RDS")
w = sesameQC_getStats(qcs[[1]], "frac_dt")

targets <- read.metharray.sheet(idat)
saveRDS(targets, "sesame_variables/targets.RDS")

RGset <- read.metharray.exp(targets=targets, extended= TRUE) 
saveRDS(RGset, "sesame_variables/RGset.RDS")

library(parallel)
v = mclapply(searchIDATprefixes("."),readIDATpair,mc.cores = 2)
saveRDS(v, "sesame_variables/v.RDS")
