

# path_rdata <- "products/_02_Rdata_NEB_es_enhancer_01_extend/rdata.Rdata"
# odir_base <- "products/_03_count_ngsplot_layered_NEB__es_enhancer_01"
# pathLabel 
# experiment_name <- "_01_ngsplot_NEB__es_enhancer_01"

args=(commandArgs(TRUE))

path_rdata = args[1]
odir_base = args[2]
pathLabel = args[3]
experiment_name = args[4]
metalabelPrefix = args[5]
pathSelectBed = args[6]
pathOutliers = args[7]

#########################
timestamp()

library(data.table); library(dplyr); library(magrittr); library(dtplyr)
library(scales)


if(!file.exists(odir_base)){
  dir.create(odir_base, recursive = T)
}

selCol = 41:120

xl <- "Position relative to enhancer center (bp)"
yl <- "Averaged coverage"
at <- c(1, 20, 40.5, 60, 80)
lb <- c(-2000, -1000, 0, 1000, 2000)

dtLabel =  fread(pathLabel, header = T)
dtOutlier = fread(pathOutliers, header=TRUE)
dtOutlier[group %in% c("G1", "S", "G2M"), label := paste0(label, "_mcf_1")]
dtLabel = dtLabel[! label %in% dtOutlier[, label]]

labels = dtLabel[, label]
groups <- dtLabel[, group]
groups_unique <- unique(groups)


##############
load(path_rdata)

##############


source("scripts/_functionsDeepTools.r")
source("scripts/_functionsHeatmap.r")

##############

str(arr1)
arr1 = arr1[,,labels]
str(arr1)
gc(); gc();

str(arr1)
arr1 = arr1[,selCol,]
str(arr1)
gc(); gc();

str(arr2)
arr2 = arr2[,selCol,]
str(arr2)
gc(); gc();

##############
if(grepl("\\.bed$", pathSelectBed)){
	dtSelectBed = fread(pathSelectBed, header = FALSE, sep = "\t")
	selRow = dtSelectBed[, V4]
}else{
	dtSelectBed = fread(pathSelectBed, header = TRUE, sep = "\t")
	selRow = dtSelectBed[, Id]
}

str(arr1[,,])
str(arr1[selRow,,])


######################
# 
# for(trimRatio in c(0)){
# 	for(maxv in c(50, 100, 150, 200)){
# 		metalabel = sprintf("%s_heatmap_%.2f_%.2f", "timeSeries", trimRatio, maxv)
# 		heatmapTrimFromArray(arr1[selRow,,], c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio, max_value = maxv)

# 		metalabel = sprintf("%s_heatmap_bg_%.2f_%.2f", "timeSeries", trimRatio, maxv)
# 		heatmapTrimFromArray(arr2, c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio, max_value = maxv)
# 	}
# }

# for(trimRatio in c(0)){
# 	for(maxv in c(50, 100, 150, 200)){
# 		metalabel = sprintf("%s_heatmap_%.2f_%.2f", "esCellCycle", trimRatio, maxv)
# 		heatmapTrimFromArray(arr1[selRow,,], c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio, max_value = maxv)

# 		metalabel = sprintf("%s_heatmap_bg_%.2f_%.2f", "esCellCycle", trimRatio, maxv)
# 		heatmapTrimFromArray(arr2, c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio, max_value = maxv)
# 	}
# }




######################
# time series

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_ratio_%.2f", "timeSeries", trimRatio)
# 	layeredRatioPlotTrim(arr1[selRow,,], arr2, c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio)
# }

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_fg_%.2f", "timeSeries", trimRatio)
# 	layeredFGPlotTrim(arr1[selRow,,], arr2, c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio)
# }

for(trimRatio in c(1, 0)){
	for(expTopRatio in c(0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0)){
		metalabel = sprintf("%s_fgNorm_%.2f_%.2f", "timeSeries", trimRatio, expTopRatio)
		layeredFGPlotTrimNormalize(arr1[selRow,,], arr2, c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio, expTopRatio=expTopRatio)
	}
}


# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_fracCount_%.2f", "timeSeries", trimRatio)
# 	layeredFracCountPlotTrim(arr1[selRow,,], c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio)

# 	# metalabel = sprintf("%s_fracCount_%.2f_neg", "timeSeries", trimRatio)
# 	# layeredFracCountPlotTrim(arr1[selRow,,], c("00hminus", "12hminus", "24hminus", "48hminus", "72hminus"), metalabel, top = trimRatio)
# }

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_count_%.2f", "timeSeries", trimRatio)
# 	layeredCountPlotTrim(arr1[selRow,,], c("00h", "12h", "24h", "48h", "72h"), metalabel, top = trimRatio)

# 	# metalabel = sprintf("%s_count_%.2f_neg", "timeSeries", trimRatio)
# 	# layeredCountPlotTrim(arr1[selRow,,], c("00hminus", "12hminus", "24hminus", "48hminus", "72hminus"), metalabel, top = trimRatio)
# }

# ######################
# # ES cell cycle

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_ratio_%.2f", "esCellCycle", trimRatio)
# 	layeredRatioPlotTrim(arr1[selRow,,], arr2, c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio)
# }
# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_fg_%.2f", "esCellCycle", trimRatio)
# 	layeredFGPlotTrim(arr1[selRow,,], arr2, c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio)
# }

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_fracCount_%.2f", "esCellCycle", trimRatio)
# 	layeredFracCountPlotTrim(arr1[selRow,,], c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio)

# 	# metalabel = sprintf("%s_fracCount_%.2f_neg", "esCellCycle", trimRatio)
# 	# layeredFracCountPlotTrim(arr1[selRow,,], c("G1minus", "G2Mminus", "Sminus", "Ntc"), metalabel, top = trimRatio)
# }

# for(trimRatio in c(1, 0)){
# 	metalabel = sprintf("%s_count_%.2f", "esCellCycle", trimRatio)
# 	layeredCountPlotTrim(arr1[selRow,,], c("G1", "S", "G2M", "AvgG1", "10pg"), metalabel, top = trimRatio)

# 	# metalabel = sprintf("%s_count_%.2f_neg", "esCellCycle", trimRatio)
# 	# layeredCountPlotTrim(arr1[selRow,,], c("G1minus", "G2Mminus", "Sminus", "Ntc"), metalabel, top = trimRatio)
# }


#####################

sessionInfo()

timestamp()

