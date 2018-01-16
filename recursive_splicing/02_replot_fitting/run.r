timestamp()

args=(commandArgs(TRUE))
idirBase = args[1]
label = args[2]
ofile = args[3]
strand = args[4]
probabilityFile = args[5]
 
# idirBase = "products/_06_seesawFitting_bin_5k"
# label = "chr6:94050766-94283009"
# ofile = ""

thPval = 1e-5
thProbability = 0.95

ifiles =  list.files(idirBase,
                      sprintf("%s.res.txt.gz", label), 
                      recursive=TRUE, full.names=TRUE)

pvalfiles = list.files(idirBase,
                       sprintf("%s.pval.txt", label), 
                       recursive=TRUE, full.names=TRUE)

path_rdata = sprintf("%s.rda", ofile)


source("scripts/r14_replot_fitting_bin_5k_coloringByPvalue/base_prepare.r", echo = TRUE)
source("scripts/r14_replot_fitting_bin_5k_coloringByPvalue/base.r", echo = TRUE)

sessionInfo()
timestamp()

