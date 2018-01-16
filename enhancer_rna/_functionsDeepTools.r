library(data.table); library(dplyr); library(magrittr); library(dtplyr)
library(scales)


collectNgsplotWithBackground = function(
            idir_base, 
            idir_base2,
            odir_base,
            fileLabel
            ){

  labels <- fread(fileLabel, header = T) %>% select(label) %>% unlist
  labels2 <- paste0(labels, "_bg")
  # groups <- gsub("_.+", "", gsub("_mcf_1", "", gsub("RamDA_", "", labels)))
  # groups_unique <- unique(groups)

  ifiles <- sprintf("%s/%s/%s/%s.%s", idir_base, labels, "ngsplot", "hm1", "txt")
  ifiles2 <- sprintf("%s/%s/%s/%s.%s", idir_base2, labels, "ngsplot", "hm1", "txt")

  for(i in seq_along(ifiles)){
    ifile <- ifiles[i]
    print(sprintf("%s, %s", labels[i], ifiles[i]))
    dt <- fread(ifile, header = T)
    if(i == 1){
      arr1 <- array(0, dim = c(nrow(dt), ncol(dt)-4, length(ifiles)))
      dimnames(arr1)[[1]] = dt[, gname]
    }
    arr1[,,i] <- dt %>% select(-gid, -gname , -tid, -strand) %>% as.matrix()
  }
  dimnames(arr1)[[3]] <- labels


  for(i in seq_along(ifiles2)){
    ifile <- ifiles2[i]
    print(sprintf("%s, %s", labels2[i], ifiles2[i]))
    dt <- fread(ifile, header = T)
    if(i == 1){
      arr2 <- array(0, dim = c(nrow(dt), ncol(dt)-4, length(ifiles2)))
      dimnames(arr2)[[1]] = dt[, gname]
    }
    arr2[,,i] <- dt %>% select(-gid, -gname , -tid, -strand) %>% as.matrix()
  }
  dimnames(arr2)[[3]] <- labels2


  path_rdata = sprintf("%s/%s.Rdata", odir_base, "rdata")
  save(arr1, arr2, file = path_rdata)
}


colTrimmedMeans <- function(mat, top){
  prob <- 1-top/100
  apply(mat, 2, function(x){
      mean(x[x <= quantile(x, prob)], na.rm = T)
    })
}

colTrimmedSe <- function(mat, top){
  prob <- 1-top/100
  apply(mat, 2, function(x){
    sqrt(var(x[x <= quantile(x, prob)], na.rm = T)/length(x[x <= quantile(x, prob)]))
  })
}

countTrimmed <- function(mat, top){
  prob <- 1-top/100
  nrow(mat) # sum(mat[,1] <= quantile(mat[,1], prob))
}

layeredPlotTrim <- function(arr, vec_groups, metalabel, top = 1, pallet = "OrRd"){
  vec_groups = as.character(vec_groups)
     # cat("1\n")
  require(biovizBase)
     # cat("2\n")

  listOfVc <- list()
  listOfVc_se <- list()
  
  N <- 0

  xmax = dim(arr)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    ma <- apply(arr[,,group == groups, drop = FALSE], 1:2, sum)

    ma <- ma / sum(group == groups)
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
   # cat("18\n")
  mn <- sprintf("%s, sum of normalized counts\n(%s)", metalabel, experiment_name)
  sb <- paste0("top ", top, "% trimmed mean (N=", N, ")")

  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_avgprof.pdf")
  pdf(pdfname)

  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(listOfVc)+1)[2:(length(listOfVc)+1)]
  print(length(listOfVc))

  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]

    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)
  }
  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]

    lines(vc, lwd = 2, col = pal[i])
  }

  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/

  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   

  dev.off()

  par(xpd=FALSE,mar=c(5,4,4,2))
}


layeredCountPlotTrim <- function(arr, vec_groups, metalabel, top = 1, pallet = "OrRd"){
  vec_groups = as.character(vec_groups)

  require(biovizBase)
  listOfVc <- list()
  listOfVc_se <- list()
  N <- 0

  xmax = dim(arr)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    arr_tmp <- arr[,,group == groups, drop = FALSE]
    arr_tmp[arr_tmp > 0] <- 1
    ma <- apply(arr_tmp, 1:2, sum)
    
    ma <- ma / sum(group == groups)
    vc <- colTrimmedMeans(ma, top = top)
    vc_se <- colTrimmedSe(ma, top = top)
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
  mn <- sprintf("%s, number of detected expression\n(%s)", metalabel, experiment_name)
  sb <- paste0("top ", top, "% trimmed mean (N=", N, ")")
  
  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_count.pdf")
  pdf(pdfname)
  
  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(listOfVc)+1)[2:(length(listOfVc)+1)]
  print(length(listOfVc))
  
  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)
  
  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)
  }
  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    
    lines(vc, lwd = 2, col = pal[i])
  }
  
  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/
  
  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   
  
  dev.off()
  
  par(xpd=FALSE,mar=c(5,4,4,2))
}


layeredFracCountPlotTrim <- function(arr, vec_groups, metalabel, top = 1, frac = 0.1, pallet = "OrRd"){
  vec_groups = as.character(vec_groups)

  require(biovizBase)
  listOfVc <- list()
  listOfVc_se <- list()
  N <- 0

  xmax = dim(arr)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    arr_tmp <- arr[,,group == groups, drop = FALSE]
    arr_tmp[arr_tmp > 0] <- 1
    ma <- apply(arr_tmp, 1:2, sum)
    ma <- ma / sum(group == groups)
    ma[ma>=frac] <- 1
    ma[ma<frac] <- 0
    
    vc <- colTrimmedMeans(ma, top = top)
    vc_se <- colTrimmedSe(ma, top = top)
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
  mn <- sprintf("%s, fraction of enchancers (%f)\n(%s)", metalabel, frac, experiment_name)
  sb <- paste0("top ", top, "% trimmed mean (N=", N, ")")
  
  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_fracCount.pdf")
  pdf(pdfname)
  
  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(listOfVc)+1)[2:(length(listOfVc)+1)]
  print(length(listOfVc))
  
  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)
  
  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)
  }
  for(i in seq_along(listOfVc)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    
    lines(vc, lwd = 2, col = pal[i])
  }
  
  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/
  
  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   
  
  dev.off()
  
  par(xpd=FALSE,mar=c(5,4,4,2))
}

layeredFGPlotTrim <- function(arr1, arr2, vec_groups, metalabel, top = 1, pallet = "OrRd"){
  vec_groups = as.character(vec_groups)
     # cat("1\n")
  require(biovizBase)
     # cat("2\n")

  listOfVc <- list()
  listOfVc_se <- list()
  
  N <- 0
  N2 <- 0

  xmax = dim(arr1)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    ma <- apply(arr1[,,group == groups, drop = FALSE], 1:2, sum)

    ma <- ma / sum(group == groups)
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
  for(group in vec_groups){
    ma <- apply(arr2[,,group == groups, drop = FALSE], 1:2, sum)

    ma <- ma / sum(group == groups)
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N2 <- countTrimmed(ma, top)
  }
   # cat("18\n")
  mn <- sprintf("%s, sum of normalized counts\n(%s)", metalabel, experiment_name)
  sb <- paste0("top ", top, "% trimmed mean (N=", N, "," , N2, ")")

  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_foreground_background.pdf")
  pdf(pdfname)

  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(vec_groups)+1)[2:(length(vec_groups)+1)]
  print(length(listOfVc))

  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

  for(i in 1:length(vec_groups)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)

    j <- i + length(vec_groups)
    vc <- listOfVc[[j]]
    vc_se <- listOfVc_se[[j]]
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)

  }
  for(i in 1:length(vec_groups)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    lines(vc, lwd = 2, col = pal[i])

    j <- i + length(vec_groups)
    vc <- listOfVc[[j]]
    vc_se <- listOfVc_se[[j]]
    lines(vc, lwd = 2, col = pal[i], lty = 2)    
  }

  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=c(vec_groups, paste0(vec_groups, "bg")), col=rep(pal,2), lwd = 2, title=NA, lty = rep(1:2, each = length(vec_groups)))  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/

  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   

  dev.off()

  par(xpd=FALSE,mar=c(5,4,4,2))
}


layeredRatioPlotTrim <- function(arr1, arr2, vec_groups, metalabel, top = 1, pallet = "OrRd"){
  vec_groups = as.character(vec_groups)
  
  require(biovizBase)

  listOfVc <- list()

  N <- 0
  N2 <- 0

  xmax = dim(arr1)[2]
  cen = (xmax + 1) /2

  for(group in vec_groups){
    ma <- apply(arr1[,,group == groups, drop = FALSE], 1:2, sum)
    ma <- ma / sum(group == groups)
    vc <- colTrimmedMeans(ma, top = top)

    ma2 <- apply(arr2[,,group == groups, drop = FALSE], 1:2, sum)
    ma2 <- ma2 / sum(group == groups)
    vc2 <- colTrimmedMeans(ma2, top = top)

    vc <- log2((vc+10^(-7))/(vc2+10^(-7)))

    listOfVc <- append(listOfVc, list(vc))
    
    N <- countTrimmed(ma, top)
    N2 <- countTrimmed(ma2, top)
  }

  yl <- "log2((Enhancer+10e-7)/(Background+10e-7))"

  mn <- sprintf("%s, sum of normalized counts\n(%s)", metalabel, experiment_name)
  sb <- paste0("top ", top, "% trimmed mean (N=", N, "," , N2, ")")

  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_ratio.pdf")
  pdf(pdfname)

  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(vec_groups)+1)[2:(length(vec_groups)+1)]
  print(length(listOfVc))

  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

  for(i in 1:length(vec_groups)){
    vc <- listOfVc[[i]]
    lines(vc, lwd = 2, col = pal[i])
  }

  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/

  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   

  dev.off()

  par(xpd=FALSE,mar=c(5,4,4,2))
}

########################
# 2017/03/19
layeredFGPlotTrimNormalize <- function(arr1, arr2, vec_groups, metalabel, top = 1, pallet = "OrRd", 
                                        expTopRatio = 0.7, rngExp = c(35:38, 42:45)){
  vec_groups = as.character(vec_groups)
     # cat("1\n")
  require(biovizBase)
     # cat("2\n")

  listOfVc <- list()
  listOfVc_se <- list()
  
  N <- 0
  N2 <- 0

  xmax = dim(arr1)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    ma <- apply(arr1[,,group == groups, drop = FALSE], 1:2, sum)
    # Divide by the number of sampels
    ma <- ma / sum(group == groups)
    # Select top enhnacers with `floor(nrow(ma) * expTopRatio)` coverage (within rngExp) 
    rsma = order(rowSums(ma[, rngExp]), decreasing = TRUE)
    expTop = floor(nrow(ma) * expTopRatio)
    rng = rsma[1:expTop]
    ma = ma[rng, ]
    # Normalized to make the sum of each enhancer to be 1 before aggreations
    ma = t(apply(ma, 1, function(x){if(sum(x)==0){x}else{x/sum(x)}}))
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
  for(group in vec_groups){
    ma <- apply(arr2[,,group == groups, drop = FALSE], 1:2, sum)
    # Divide by the number of sampels
    ma <- ma / sum(group == groups)
    # # Select top enhnacers with `floor(nrow(ma) * expTopRatio)` coverage (within rngExp) 
    # rsma = order(rowSums(ma[, rngExp]), decreasing = TRUE)
    # expTop = floor(nrow(ma) * expTopRatio)
    # rng = rsma[1:expTop]
    # ma = ma[rng, ]
    # Normalized to make the sum of each enhancer to be 1 before aggreations
    ma = t(apply(ma, 1, function(x){if(sum(x)==0){x}else{x/sum(x)}}))
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N2 <- countTrimmed(ma, top)
  }
   # cat("18\n")
  mn <- sprintf("%s, sum of normalized counts\n(%s)", metalabel, experiment_name)
  sb <- sprintf("Top %.1f%% trimmed mean, top %.1f RamDA expression (N=%d,%d)", top, expTopRatio * 100, N, N2)

  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_foreground_background.pdf")
  pdf(pdfname)

  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(vec_groups)+1)[2:(length(vec_groups)+1)]
  print(length(listOfVc))

  par(xpd=TRUE,mar=c(5,4,4,10))
  plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

  for(i in 1:length(vec_groups)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)

    j <- i + length(vec_groups)
    vc <- listOfVc[[j]]
    vc_se <- listOfVc_se[[j]]
    polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)

  }
  for(i in 1:length(vec_groups)){
    vc <- listOfVc[[i]]
    vc_se <- listOfVc_se[[i]]
    lines(vc, lwd = 2, col = pal[i])

    j <- i + length(vec_groups)
    vc <- listOfVc[[j]]
    vc_se <- listOfVc_se[[j]]
    lines(vc, lwd = 2, col = pal[i], lty = 2)    
  }

  xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
  ycoord <- ymax
  legend(xcoord, ycoord, legend=c(vec_groups, paste0(vec_groups, "bg")), col=rep(pal,2), lwd = 2, title=NA, lty = rep(1:2, each = length(vec_groups)))  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/

  box()
  axis(1, at = at, labels = lb)
  axis(2)
  abline(v=cen, col = alpha("grey", 0.2))   

  dev.off()

  par(xpd=FALSE,mar=c(5,4,4,2))
}



layeredFGPlotTrimNormalizeSeparateFG <- function(arr1, arr2, vec_groups, metalabel, top = 1, pallet = "OrRd", 
                                        expTopRatio = 0.7, rngExp = c(35:38, 42:45)){
  vec_groups = as.character(vec_groups)
     # cat("1\n")
  require(biovizBase)
     # cat("2\n")

  listOfVc <- list()
  listOfVc_se <- list()
  
  N <- 0
  N2 <- 0

  xmax = dim(arr1)[2]
  cen = (xmax + 1)/2

  for(group in vec_groups){
    ma <- apply(arr1[,,group == groups, drop = FALSE], 1:2, sum)
    # Divide by the number of sampels
    ma <- ma / sum(group == groups)
    # Select top enhnacers with `floor(nrow(ma) * expTopRatio)` coverage (within rngExp) 
    rsma = order(rowSums(ma[, rngExp]), decreasing = TRUE)
    expTop = floor(nrow(ma) * expTopRatio)
    rng = rsma[1:expTop]
    ma = ma[rng, ]
    # Normalized to make the sum of each enhancer to be 1 before aggreations
    ma = t(apply(ma, 1, function(x){if(sum(x)==0){x}else{x/sum(x)}}))
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N <- countTrimmed(ma, top)
  }
  for(group in vec_groups){
    ma <- apply(arr2[,,group == groups, drop = FALSE], 1:2, sum)
    # Divide by the number of sampels
    ma <- ma / sum(group == groups)
    # # Select top enhnacers with `floor(nrow(ma) * expTopRatio)` coverage (within rngExp) 
    # rsma = order(rowSums(ma[, rngExp]), decreasing = TRUE)
    # expTop = floor(nrow(ma) * expTopRatio)
    # rng = rsma[1:expTop]
    # ma = ma[rng, ]
    # Normalized to make the sum of each enhancer to be 1 before aggreations
    ma = t(apply(ma, 1, function(x){if(sum(x)==0){x}else{x/sum(x)}}))
     # cat("14\n")
    vc <- colTrimmedMeans(ma, top = top)
     # cat("15\n")
    vc_se <- colTrimmedSe(ma, top = top)
       # cat("16\n")
    listOfVc <- append(listOfVc, list(vc))
    listOfVc_se <- append(listOfVc_se, list(vc_se))
    
    N2 <- countTrimmed(ma, top)
  }
   # cat("18\n")
  mn <- sprintf("%s, sum of normalized counts\n(%s)", metalabel, experiment_name)
  sb <- sprintf("Top %.1f%% trimmed mean, top %.1f RamDA expression (N=%d,%d)", top, expTopRatio * 100, N, N2)

  ymax <- extendrange(unlist(listOfVc))[2]
  pal <- colorBlindSafePal(pallet)(length(vec_groups)+1)[2:(length(vec_groups)+1)]
  print(length(listOfVc))

  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".layered_foreground_background.pdf")
  pdf(pdfname)

  par(xpd=TRUE,mar=c(5,4,4,10))

  # Foreground
  {
    plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

    for(i in 1:length(vec_groups)){
      vc <- listOfVc[[i]]
      vc_se <- listOfVc_se[[i]]
      polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)
    }
    for(i in 1:length(vec_groups)){
      vc <- listOfVc[[i]]
      vc_se <- listOfVc_se[[i]]
      lines(vc, lwd = 2, col = pal[i])
    }

    xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
    ycoord <- ymax
    legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA, lty = 1)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/
    box()
    axis(1, at = at, labels = lb)
    axis(2)
    abline(v=cen, col = alpha("grey", 0.2))
  }

  {
    plot(listOfVc[[1]], type = "n", xlab = xl, axes = F, ylim = c(0, ymax), ylab = yl,  main = mn, sub = sb)

    for(i in 1:length(vec_groups)){
      j <- i + length(vec_groups)
      vc <- listOfVc[[j]]
      vc_se <- listOfVc_se[[j]]
      polygon(x = c(1:xmax, xmax:1), y = c(vc + vc_se, rev(vc - vc_se)), col = alpha(pal[i], 0.3), border = NA)

    }
    for(i in 1:length(vec_groups)){
      j <- i + length(vec_groups)
      vc <- listOfVc[[j]]
      vc_se <- listOfVc_se[[j]]
      lines(vc, lwd = 2, col = pal[i])    
    }

    xcoord <- max(1:xmax)+((max(1:xmax)-min(1:xmax))/20)
    ycoord <- ymax
    legend(xcoord, ycoord, legend=vec_groups, col=pal, lwd = 2, title=NA, lty = 1)  # http://www.r-bloggers.com/adding-color-to-r-plot-a-function/

    box()
    axis(1, at = at, labels = lb)
    axis(2)
    abline(v=cen, col = alpha("grey", 0.2))   

  }


  dev.off()

  par(xpd=FALSE,mar=c(5,4,4,2))
}



