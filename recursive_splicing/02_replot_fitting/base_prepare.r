
###############

library(data.table)
library(ggplot2)
library(dplyr)
library(dtplyr)
library(magrittr)

#############

if(! file.exists(dirname(ofile))) dir.create(dirname(ofile))

################
# Functions
checkIfSawtooth = function(Strands, Intercepts, Bps, Rsss, BpsBefore){
  # Added new condition: after/before slope ratio > 1
  sapply(seq_along(Strands), function(i){
    Strand = Strands[i]
    Intercept = Intercepts[i]
    Bp = Bps[i]
    Rss = Rsss[i]
    BpBefore = BpsBefore[i]
    if(is.na(Rss)){
      return(NA)
    }
  
      
    if(Strand == "+"){
#     if(Intercept > 0 & Bp < 0 & Rss < 0 & Bp/BpBefore > 1){
      if(Intercept > 0 & Bp < 0 & Rss < 0 & Bp < BpBefore){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
    if(Strand == "-"){
      # if(Intercept < 0 & Bp > 0 & Rss > 0 & Bp/BpBefore > 1){
      if(Intercept < 0 & Bp > 0 & Rss > 0 & Bp > BpBefore){
        return(TRUE)
      }else{
        return(FALSE)
      }
      
    }
    return(NA)
  })
}

#################
# Load data and preprocess
dt = rbindlist(
  lapply(ifiles, function(ifile){
    dt = fread(sprintf("gzip -dc %s", ifile), header=TRUE)
    dt[, path := ifile]
    dt
  })
)
dt[, group := basename(dirname(path))]


dtpval = rbindlist(
  lapply(pvalfiles, function(ifile){
    dt = fread(ifile, header=TRUE)
    dt[, path := ifile]
    dt
  })
)
dtpval[, group := basename(dirname(path))]

dtpval[, strand := strand]
dtpval[, sawtooth_raw := checkIfSawtooth(strand, df.rss.lm_Intercept, df.rss.lm_bp, df.rss.lm_RSSTRUE, df.lm_bp)] %>% invisible()
dtpval[, pval.raw_new := pval.raw]
dtpval[sawtooth_raw==FALSE, pval.raw_new := 1]

dtpval[, isRs := pval.raw_new < thPval]

# Select high probability samples
dtProb = fread(probabilityFile, header=TRUE, sep="\t")
dtProb = dtProb[, .(label, rsDetectionProb, label_sum_down, label_sum_up)]
dtpval = dtpval[group %in% dtProb[rsDetectionProb > thProbability & label_sum_up > 0 & label_sum_down > 0, label]]
print(str(dtpval))

save(dt, dtpval, file=path_rdata)

########################