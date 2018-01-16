# The main part of this script is derived from Nikaido-san

# argsments: ifile, oprefix, id, binsize

timestamp()

#############################################
library(data.table)

print(ifile)
print(oprefix)
print(id)
print(binsize)

#########################################
# Calculate whether RSS is significant
#########################################
# Load data: assuming[position (Int), depth (Float), RSS (Bool)]
dt = fread(sprintf("gzip -dc %s", ifile), header=TRUE)
# dt[, depth.sm := runmed(depth, k = 17)]
# df = as.data.frame(dt)

# Binning
breaks = seq(dt[, min(bp)], dt[, max(bp)] + binsize, by = binsize)
bpbin = .bincode(dt$bp, breaks, right=FALSE)
dt[, bpmed := breaks[bpbin] + floor(binsize/2)]

dt2 = copy(dt)
dt2 = dt2[, depthMean := mean(depth), by = bpmed]
dt2 = unique(dt2[, .(bpmed, depthMean)], by=NULL)
dim(dt2)
posRSS = dt[RSS == TRUE, max(bp)]
dt2[, RSS := (bpmed < posRSS)]

# Encode relative distance (the distance from the closest upstream RS site candidate)
if(strand == "+"){
        dt2[bpmed < posRSS, bpmed_relative := bpmed-min(bpmed)]
        dt2[bpmed >= posRSS, bpmed_relative := bpmed-min(bpmed)]
}else{
        dt2[bpmed < posRSS, bpmed_relative := max(bpmed)-bpmed]
        dt2[bpmed >= posRSS, bpmed_relative := max(bpmed)-bpmed]      
}

setnames(dt2, "bpmed", "bp")
setnames(dt2, "bpmed_relative", "bp_relative")
setnames(dt2, "depthMean", "depth")
dt2[, depth.sm := runmed(depth, k = 17)]

df = as.data.frame(dt2)


# Perform linear regression
df.lm    = lm(depth ~ bp, data = df)
df.rss.lm = lm(depth ~ bp + RSS, data = df)

df.rel.lm = lm(depth ~ bp_relative, data = df)

df.sm.lm = lm(depth.sm ~ bp, data = df)
df.sm.rss.lm = lm(depth.sm ~ bp + RSS, data = df)

# Obtain fitted values from the estimated paramter values
y.pred    = fitted(df.lm)
y.rss.pred = fitted(df.rss.lm)

y.rel.pred = fitted(df.rel.lm)

y.sm.pred = fitted(df.sm.lm)
y.sm.rss.pred = fitted(df.sm.rss.lm)

# Peform significance test for the slope
lm.anv = anova(df.lm, df.rss.lm)
print(lm.anv)

aic.lm = AIC(df.lm)
aic.lm.rel = AIC(df.rel.lm)

# https://stats.stackexchange.com/questions/8513/test-equivalence-of-non-nested-models/8519#8519
## The Cox Test for Comparing Non-Nested Models
require(lmtest)
pval_cox_2_over_1 = coxtest(df.lm, df.rel.lm)[[4]][1]
pval_cox_1_over_2 = coxtest(df.lm, df.rel.lm)[[4]][2]
print(pval_cox_2_over_1)
print(pval_cox_1_over_2)


lm.sm.anv = anova(df.sm.lm, df.sm.rss.lm)
print(lm.sm.anv)

pval.raw = lm.anv[2,6]
pval.sm = lm.sm.anv[2,6]


##############################
# Plot
##############################

labels3 = c("raw", "baseline (raw)", "RSS (raw)", "rel_dist (raw)")
n3 = length(labels3)
title.txt3 = sprintf("%s\np_raw=%.3e p_cox=%.3e dAIC(baseline-rel_dist)=%.1f", id, pval.raw, pval_cox_2_over_1, aic.lm - aic.lm.rel)

labels2 = c("smooth", "predict (smooth)", "RSS (smooth)")
n2 = length(labels2)

labels = c("raw", "smooth", "predict (raw)", "predict (smooth)", "RSS (raw)", "RSS (smooth)")
title.txt = sprintf("%s\np_raw=%.3e p_smooth=%.3e", id, pval.raw, pval.sm)
n = length(labels)
ymax = max(c(df$depth.sm, df$depth)) * 2

###
pdfname = sprintf("%s.pdf", oprefix)
pdf(pdfname)

matplot(df$bp, cbind(df$depth.sm, y.sm.pred, y.sm.rss.pred),
        type =  c("s", rep("l",(n2-1))), 
        lty = n2:1, col = 1:n2, 
        main = title.txt,
        lwd = c(1, rep(2, n2-1)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels2, col = 1:n2, lty = n2:1, lwd = c(1, rep(2, n2-1)))
rug(posRSS, lwd = 1)

matplot(df$bp, cbind(df$depth, df$depth.sm, y.pred, y.sm.pred, y.rss.pred, y.sm.rss.pred),
        type =  c("s", "s", rep("l",(n-2))), 
        lty = n2:1, col = 1:n, 
        main = title.txt,
        lwd = c(1, 1, rep(2, n-2)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels, col = 1:n, lty = n:1, lwd = c(1, 1, rep(2, n-2)))
rug(posRSS, lwd = 1)

matplot(df$bp, cbind(df$depth, y.pred, y.rss.pred, y.rel.pred),
        type =  c("s", "s", rep("l",(n3-1))), 
        lty = n3:1, col = 1:n3, 
        main = title.txt3,
        lwd = c(1, rep(2, n3-1)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels3, col = 1:n3, lty = n3:1, lwd = c(1, rep(2, n3-1)))
rug(posRSS, lwd = 1)

dev.off()

###
pngname = sprintf("%s.smooth.png", oprefix)
png(pngname)

matplot(df$bp, cbind(df$depth.sm, y.sm.pred, y.sm.rss.pred),
        type =  c("s", rep("l",(n2-1))), 
        lty = n2:1, col = 1:n2, 
        main = title.txt,
        lwd = c(1, rep(2, n2-1)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels2, col = 1:n2, lty = n2:1, lwd = c(1, rep(2, n2-1)))
rug(posRSS, lwd = 1)

dev.off()

###
pngname = sprintf("%s.rawSmooth.png", oprefix)
png(pngname)

matplot(df$bp, cbind(df$depth, df$depth.sm, y.pred, y.sm.pred, y.rss.pred, y.sm.rss.pred),
        type =  c("s", "s", rep("l",(n-2))), 
        lty = n2:1, col = 1:n, 
        main = title.txt,
        lwd = c(1, 1, rep(2, n-2)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels, col = 1:n, lty = n:1, lwd = c(1, 1, rep(2, n-2)))
rug(posRSS, lwd = 1)

dev.off()

###
pngname = sprintf("%s.raw_rel_dist.png", oprefix)
png(pngname)

matplot(df$bp, cbind(df$depth, y.pred, y.rss.pred, y.rel.pred),
        type =  c("s", "s", rep("l",(n3-1))), 
        lty = n3:1, col = 1:n3, 
        main = title.txt3,
        lwd = c(1, rep(2, n3-1)),
        ylim = c(0, ymax),
        ylab = "Depth of reads", xlab = "bp")
legend("topleft", legend = labels3, col = 1:n3, lty = n3:1, lwd = c(1, rep(2, n3-1)))
rug(posRSS, lwd = 1)

dev.off()

#############################
# Save
#############################
dt2[, y.pred := y.pred]
dt2[, y.rss.pred := y.rss.pred]
dt2[, y.rel.pred := y.rel.pred]
dt2[, y.sm.pred := y.sm.pred]
dt2[, y.sm.rss.pred := y.sm.rss.pred]

ofile = sprintf("%s.res.txt.gz", oprefix)
gzHandler = gzfile(ofile, "w")
write.table(dt2, gzHandler, sep = "\t", quote=FALSE, row.names=FALSE)
close(gzHandler)

ofile = sprintf("%s.pval.txt", oprefix)
write.table(data.table(
						id = id, 
						pval.raw = pval.raw,
                                                pval.raw.rel = pval_cox_2_over_1,
                                                pval.raw.rel_reverse = pval_cox_1_over_2,                                                
                                                aic_lm = aic.lm,
                                                aic_lm_rel = aic.lm.rel,
						pval.sm = pval.sm, 
						df.lm_Intercept = coef(df.lm)[1],
						df.lm_bp = coef(df.lm)[2],
						df.rss.lm_Intercept = coef(df.rss.lm)[1],
						df.rss.lm_bp = coef(df.rss.lm)[2],
						df.rss.lm_RSSTRUE = coef(df.rss.lm)[3],
                                                df.rel.lm_Intercept = coef(df.rel.lm)[1],
                                                df.rel.lm_bp = coef(df.rel.lm)[2],
						df.sm.lm_Intercept = coef(df.sm.lm)[1],
						df.sm.lm_bp = coef(df.sm.lm)[2],
						df.sm.rss.lm_Intercept = coef(df.sm.rss.lm)[1],
						df.sm.rss.lm_bp = coef(df.sm.rss.lm)[2],
						df.sm.rss.lm_RSSTRUE = coef(df.sm.rss.lm)[3]
					),
	ofile, sep = "\t", quote=FALSE, row.names=FALSE)



##########################################
sessionInfo()

timestamp()

