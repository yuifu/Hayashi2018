library(data.table); library(dplyr); library(magrittr); library(dtplyr)
library(scales)

##############

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
  nrow(mat) # sum(mat[,1] < quantile(mat[,1], prob))
}

sumMatrix <- function(sublabels){
  flag <- 1
  # cat("7\n")
  for(label in sublabels){
    # cat("8\n")
    ifile <- sprintf("%s/%s/%s/%s.%s", idir_base, label, "ngsplot", "hm1", "txt")
    print(label)
    # ifiles <- sprintf("%s/%s/%s/%s.%s", idir_base, sublabels, sublabels, "avgprof", "txt")
    # cat("9\n")
    dt <- fread(ifile, header = T)
    # cat("10\n")
    if(flag == 1){
      ma <- matrix(0, nrow = dt %>% nrow, ncol = (dt %>%ncol)-4)
      print("Reset ma")
      flag <- 0
    }    
    # if(!exists("ma")){
    #   ma <- matrix(0, nrow = dt %>% nrow, ncol = (dt %>%ncol)-4)
    #   print("Reset ma")
    # }    
    # cat("11\n")
    ma <- ma + dt %>% select(-gid, -gname , -tid, -strand) %>% as.matrix()
    # cat("12\n")
    #   vc <- vc + dt %>% unlist
  }
  return(ma)
}

meanMatrix <- function(sublabels){
  flag <- 1
  # cat("7\n")
  for(label in sublabels){
    # cat("8\n")
    ifile <- sprintf("%s/%s/%s/%s.%s", idir_base, label, "ngsplot", "hm1", "txt")
    print(label)
    # ifiles <- sprintf("%s/%s/%s/%s.%s", idir_base, sublabels, sublabels, "avgprof", "txt")
    # cat("9\n")
    dt <- fread(ifile, header = T)
    # cat("10\n")
    if(flag == 1){
      ma <- matrix(0, nrow = dt %>% nrow, ncol = (dt %>%ncol)-4)
      print("Reset ma")
      flag <- 0
    }    
    # if(!exists("ma")){
    #   ma <- matrix(0, nrow = dt %>% nrow, ncol = (dt %>%ncol)-4)
    #   print("Reset ma")
    # }    
    # cat("11\n")
    ma <- ma + dt %>% select(-gid, -gname , -tid, -strand) %>% as.matrix()
    # cat("12\n")
    #   vc <- vc + dt %>% unlist
  }
  ma <- ma/length(sublabels)
  return(ma)
}

heatmapTrim <- function(vec_groups, metalabel, top = 1, max_value = NA){
  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".heatmap.pdf")
  pdf(pdfname)
  
  for(group in vec_groups){
    sublabels <- labels[group == groups]
    print(group)
    
    # ma <- sumMatrix(sublabels)
    ma <- meanMatrix(sublabels)
    sel <- order(rowSums(ma, na.rm = T), decreasing = T)
    sel <- sel[seq_along(sel) > nrow(ma) * (top/100)]
    print(min(sel))
    print(max(ma[sel,]))

    N <- length(sel)

    mn <- sprintf("%s, %s (%s)", metalabel, group, experiment_name)
    sb <- paste0("top ", top, "% trimmed mean (N=", N, ")")
    plotHeatmap(ma[sel, ], main = mn, sub = sb, max_value = max_value, min_value = NA)
  }
  dev.off()
}

heatmapTrimFromArray <- function(arr, vec_groups, metalabel, top = 1, max_value = NA){
  pdfname <- sprintf("%s/%s%s", odir_base, metalabel, ".heatmap.pdf")
  pdf(pdfname)
  
  for(group in vec_groups){
    ma <- apply(arr[,,group == groups, drop = FALSE], 1:2, mean)
    print(group)
    
    sel <- order(rowSums(ma, na.rm = T), decreasing = T)
    sel <- sel[seq_along(sel) > nrow(ma) * (top/100)]
    print(min(sel))
    print(max(ma[sel,]))

    N <- length(sel)

    mn <- sprintf("%s, %s (%s)", metalabel, group, experiment_name)
    sb <- paste0("top ", top, "% trimmed mean (N=", N, ")")
    plotHeatmap(ma[sel, ], main = mn, sub = sb, max_value = max_value, min_value = NA)
  }
  dev.off()
}


cut2 <- function (x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE, 
                  dig.lab = 3L, ordered_result = FALSE, fixed_max, ...) 
{
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (length(breaks) == 1L) {
    if (is.na(breaks) || breaks < 2L) 
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    rx[2L] <- fixed_max
    if (dx == 0) {
      dx <- abs(rx[1L])
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
                        length.out = nb)
    }
    else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
                               dx/1000)
    }
  }
  else nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks)) 
    stop("'breaks' are not unique")
  codes.only <- FALSE
  if (is.null(labels)) {
    for (dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(0 + breaks, digits = dig, width = 1L)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
        break
    }
    labels <- if (ok) 
      paste0(if (right) 
        "("
        else "[", ch.br[-nb], ",", ch.br[-1L], if (right) 
          "]"
        else ")")
    else paste("Range", seq_len(nb - 1L), sep = "_")
    if (ok && include.lowest) {
      if (right) 
        substr(labels[1L], 1L, 1L) <- "["
      else substring(labels[nb - 1L], nchar(labels[nb - 
                                                     1L], "c")) <- "]"
    }
  }
  else if (is.logical(labels) && !labels) 
    codes.only <- TRUE
  else if (length(labels) != nb - 1L) 
    stop("lengths of 'breaks' and 'labels' differ")
  code <- .bincode(x, breaks, right, include.lowest)
  if (codes.only) 
    code
  else factor(code, seq_along(labels), labels, ordered = ordered_result)
}


plotHeatmap <- function(mat, main, sub, max_value = NA, min_value = NA){
  require(grid)
  
  breaks <- 100
  pal <- colorRampPalette(c('white','blue'))
  pal <- pal(breaks)
  
  if(is.na(max_value)){
    max_value <- max(mat, na.rm = T)  
  }
  if(is.na(min_value)){  
    min_value <- min(mat, na.rm = T)
  }

  roworder <- order(rowSums(mat, na.rm = T), decreasing = T)
  
  grid.newpage()
  pushViewport(plotViewport(c(1,1,1,1)))
  col_ratio <- unit(c(1,1,4,1), c("null","cm", "null", "null"))
  pushViewport(viewport(layout = grid.layout(5, length(col_ratio), heights = unit(c(1, 1, 2,12, 1), c("lines", "lines", "lines", "null", "line") ), widths=col_ratio)))
  
  pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 3))
  plotHeatmap2(mat[roworder,], breaks, pal, max_value)
  grid.rect(gp = gpar(fill = NA, col = "grey", alpha = 0.5))
  pushViewport(viewport(xscale = c(min(at), max(at))))
  grid.xaxis(gp = gpar(cex = 0.5), at = at, label = lb)
  popViewport()     
  
  popViewport()     
  
  pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3))
  plotColorKey(pal[1:breaks], min_value, max_value)
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
  grid.text(main)
  popViewport()   
  
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
  grid.text(sub, gp = gpar(cex = 0.8))
  popViewport()   
  
  popViewport()   
}


plotHeatmap2 <- function(mat, breaks, pal, fixed_max){
  nr <- nrow(mat)
  grid.raster(matrix(pal[as.numeric(cut2(mat, breaks, fixed_max = fixed_max))], nrow = nr),
              interpolate = F, width =1, height = 1)
}
plotColorKey <- function(colors, min_color, max_color){
  pushViewport(viewport(layout = grid.layout(
    1,2,
    widths = unit(c(4,1), c("null","null"))
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  
  
  grid.raster(matrix(colors,nrow=1), y = 0.25, interpolate = F, width = 1, height = 0.5, vjust = 0)
  grid.rect(x=0, y = 0.25, width=1, height = 0.5, gp = gpar(fill = NA, col="black"), vjust = 0, hjust = 0)
  
  grid.text(min_color, x=0, y=0.24, hjust = 0, vjust = 1, gp = gpar(cex = 0.5))
  grid.text(sprintf("%.3f", max_color), x=1, y=0.24, hjust = 1, vjust = 1, gp = gpar(cex = 0.5))
  
  popViewport()
  popViewport()
}





