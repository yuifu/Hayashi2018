
timestamp()

args=(commandArgs(TRUE))
ifile = args[1]
oprefix = args[2]
id = args[3]
strand = args[4]

#
binsize = 5000

source("scripts/_06_seesawFitting_bin_5k/base_new_binsize.r", echo=TRUE)

sessionInfo()
timestamp()

