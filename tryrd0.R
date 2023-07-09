#!/usr/bin/env Rscript
# this script does what? try stuff with readRDS()
library(ggplot2)
library(Cairo)


bcdir <- "/mnt/sdb1/ncrepros/mayo-brain-mets/data/norm-counts/091/batch-corrected/"
renv <- load(paste0(bcdir, "bm-log2-TMMCPM-protein-batch-corrected.RData"))
alt <- readRDS(paste0(bcdir, "bm.log2CPMTMM.filt.pc.corrected.rds"))

saveRDS(bm.log2CPMTMM.pc.corrected, paste0(bcdir, "bmcpmcorr.rds"), compress=F)
saveRDS(alt, paste0(bcdir, "alt_bmcpmcorr.rds"), compress=F)

# test sized:
a1 <- bm.log2CPMTMM.pc.corrected[1:4,1:4]
a2 <- alt[1:4,1:4]
saveRDS(a1, paste0(bcdir, "a1.rds"), compress=F)
a2[1,1] <- 0.125
a2[1,2] <- 0.125
# a2[1,3] <- 0.125
saveRDS(a2, paste0(bcdir, "a2.rds"), compress=F)

write.csv(a1, paste0(bcdir, "a1.csv"), quote=F)
write.csv(a2, paste0(bcdir, "a2.csv"), quote=F)
write.csv(bm.log2CPMTMM.pc.corrected, paste0(bcdir, "bmcpmcorr.csv"), quote=F)
write.csv(alt, paste0(bcdir, "alt_bmcpmcorr.csv"), quote=F)

write.csv(matrix(sprintf("%a", a1),byrow=T, nrow=nrow(a1)), paste0(bcdir, "a1_p.csv"), quote=F)
write.csv(matrix(sprintf("%a", a2),byrow=T, nrow=nrow(a2)), paste0(bcdir, "a2_p.csv"), quote=F)

# Cairo image template
# CairoPNG("fname.png", 800, 800)
# put plot command here
# dev.off()
