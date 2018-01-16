pathRegion = "../170825/products/_00_make_list/partial_novel.annotateWithReferenceIntron.noOverlapWithSmallerIntron.filterOutlierCells.sum.over5kbp.seesaw_p1.00e-05.raw.txt"


idirBase = "../170828/report/results__single_cell_sawtooth_fitting/"
idirBaseDownsample = "../170828/report/results_sensitivity_of_RS_detection/" # RsDetectionRate_Cadm1_00h_chr9/47634076-47787970.txt

selGene = c("Robo2", "Magi1", "Cadm1")
selGroup = c("72h", "00h", "00h")

thpval = 10^(-5)

odir = "products/r00_calcProbForRsDetection"

source("scripts/r00_calcProbForRsDetection/base.R", echo=TRUE)
