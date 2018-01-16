# idir_base = "products/_12_ngsplot_NEBsensAndStrand/all/Tier3"
# idir_base2= "products/_12_ngsplot_NEBsensAndStrand_bg/all/Tier3"
# odir_base = "products/_13_aggregationPlot/NEBsensAndStrand/all/Tier3"
# fileLabel = "data/files_to_analyze.txt"

#########################
timestamp()


args=(commandArgs(TRUE))

idir_base = args[1]
idir_base2= args[2]
odir_base = args[3]
fileLabel = args[4]

##########


library(data.table); library(dplyr); library(magrittr); library(dtplyr)
library(scales)


if(!file.exists(odir_base)){
  dir.create(odir_base, recursive = T)
}

# source("scripts/_functions.r")

collectDeepToolsWithBackground = function(
            idir_base, 
            idir_base2,
            odir_base,
            fileLabel
            ){

  labels <- fread(fileLabel, header = T) %>% select(label) %>% unlist
  labels2 <- paste0(labels, "_bg")
  # groups <- gsub("_.+", "", gsub("_mcf_1", "", gsub("RamDA_", "", labels)))
  # groups_unique <- unique(groups)

  ifiles <- sprintf("%s/%s/%s.%s", idir_base, labels, labels, "txt.gz")
  ifiles2 <- sprintf("%s/%s/%s.%s", idir_base2, labels, labels, "txt.gz")

  for(i in seq_along(ifiles)){
    ifile <- ifiles[i]
    print(sprintf("%s, %s", labels[i], ifiles[i]))
    dt <- fread(sprintf("gzip -dc %s", ifile), header = FALSE, skip = 1)
    if(i == 1){
      arr1 <- array(0, dim = c(nrow(dt), ncol(dt)-6, length(ifiles)))
      dimnames(arr1)[[1]] = dt[, V4]
    }
    ma = dt %>% select(-V1,-V2,-V3,-V4,-V5,-V6) %>% as.matrix()
    ma[is.nan(ma)] = 0
    arr1[,,i] <- ma
  }
  dimnames(arr1)[[3]] <- labels


  for(i in seq_along(ifiles2)){
    ifile <- ifiles2[i]
    print(sprintf("%s, %s", labels2[i], ifiles2[i]))
    dt <- fread(sprintf("gzip -dc %s", ifile), header = FALSE, skip = 1)
    if(i == 1){
      arr2 <- array(0, dim = c(nrow(dt), ncol(dt)-6, length(ifiles2)))
      dimnames(arr2)[[1]] = dt[, V4]
    }
    ma = dt %>% select(-V1,-V2,-V3,-V4,-V5,-V6) %>% as.matrix()
    ma[is.nan(ma)] = 0
    arr2[,,i] <- ma
  }
  dimnames(arr2)[[3]] <- labels2


  path_rdata = sprintf("%s/%s.Rdata", odir_base, "rdata")
  save(arr1, arr2, file = path_rdata)
}


collectDeepToolsWithBackground(idir_base, idir_base2, odir_base, fileLabel)


sessionInfo()

timestamp()
