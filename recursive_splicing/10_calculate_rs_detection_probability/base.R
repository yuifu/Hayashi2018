
#########################
library(data.table); library(dplyr); library(magrittr); library(dtplyr)
library(ggplot2)
library(tidyr)
library(knitr)

# source("glm_and_ecdf.R")

if(!file.exists(odir)){
	  dir.create(odir, recursive = T)
}
############################

calcRsDetectionProb = function(oldx, oldy, newx, newy, label = "", n_replicate = 100000, outPrefix = "./plot"){
	
	# Prepare data for baseline logistic regression
	intron_reads = oldx
	y = oldy # binary (0, 1)
	
	# Learn logistic regression model
	model = glm(y~intron_reads, family=binomial)
	summary(model) %>% print
	# plot(intron_reads, y, main = label)
	# lines(sort(intron_reads), predict(model, data.frame(intron_reads = sort(intron_reads)), type="response"))

			pdf(sprintf("%s_logitRegression.pdf", outPrefix), width=8, height=4)
			plot(intron_reads, y, main = label)
			lines(sort(intron_reads), predict(model, data.frame(intron_reads = sort(intron_reads)), type="response"))
			dev.off()
	
	dttmp = data.table(
		x = sort(intron_reads),
		y = predict(model, data.frame(intron_reads = sort(intron_reads)), type="response")
		)
	ofile = sprintf("%s_logitRegression.txt", outPrefix)
	write.table(dttmp, ofile, sep = "\t", quote=FALSE, row.names=FALSE)	


	g = ggplot(data.table(x=oldx, y=oldy), aes(x=x, y=y))
	# g = g + geom_point()
	g = g + stat_smooth(method="glm", se=F, method.args = list(family="binomial"), colour="black", size=0.5)
	g = g + ggtitle(basename(outPrefix)) + xlab("log10(Fraction of intronic reads relative to aggregated timepoint data)") + ylab("RS detection probability")
	g = g + theme_bw()
	for(p in c(0.8, 0.85, 0.9, 0.95, 0.99) ){
		px = (log(p/(1-p)) -  model$coefficients[[1]])/ model$coefficients[[2]]
		g = g + geom_vline(xintercept = px, linetype="dotted", 
                color = "black", size=0.25)
	}
	pdf(sprintf("%s_logitRegression.ggplot.pdf", outPrefix), width=6, height=3)
	plot(g)
	dev.off()
		

	# Calculate probability of RS detection
	# given the observed intron reads and the learned logistic regression model
	newdata = data.frame(intron_reads = newx)
	newp = predict(model, newdata, type="response")
	return(newp)
	
	# print(str(newp))
	
	# # Calculate empirical cumulative distribution function for the number of samples with RS detection
	# print(sprintf("n_replicate: %d", n_replicate))
	# vec_replicate = numeric(n_replicate)
	# for(n in 1:n_replicate){
	# 	vec_replicate[n] = sum(newp >= runif(length(newp), 0, 1))
	# }
	# # print(summary(vec_replicate))
	# p_RS = ecdf(vec_replicate)
	# 
	# # Test
	# test_stat = sum(newy)
	# bl = min(min(vec_replicate), test_stat)
	# br = max(max(vec_replicate), test_stat)
	# subText = sprintf("p=%.3e", ifelse(p_RS(test_stat)<0.5, p_RS(test_stat), 1- p_RS(test_stat)))
	# 
	# hist(vec_replicate, main = label, breaks=seq(bl, br, length.out = 30), sub = subText, xlab = "Number of samples with RS detection")
	# abline(v=test_stat,col="black")
	# 
	# 		pdf(sprintf("%s_histEmpricalDistirbution.pdf", outPrefix), width=8, height=4)
	# 		hist(vec_replicate, main = label, breaks=seq(bl, br, length.out = 30), sub = subText, xlab = "Number of samples with RS detection")
	# 		abline(v=test_stat,col="black")
	# 		dev.off()
	# 
	# 	
	# result = c(test_stat = test_stat, lower = p_RS(test_stat), greater = 1-p_RS(test_stat))
	# 
	# return(result)

}



#########################
main00 = function(i){

	id = dtRegion[i, name]
	gn = dtRegion[i, gene_name]
	gid = dtRegion[i, gene_id]
	gr = selGroup[which(gn == selGene)]

	ifile = sprintf("%s/tidy_%s_%s_%s%s", idirBase, gn, gr, id, ".txt")
	
	ifileDownsample = sprintf("%s/RsDetectionRate_%s_%s_%s%s", 
														idirBaseDownsample, gn, gr, id, ".txt")
	

	dt1 = fread(ifile, header=TRUE)

	# dt1 = merge(dt, dtTime[, .(label, pseudotime)], by = "label")
	
	# ofile = sprintf("%s/tidy_pseudotime_%s_%s_%s%s", odir, gn, gr, id, ".txt")
	# write.table(dt1, ofile, sep = "\t", quote=FALSE, row.names=FALSE)


	titleText = sprintf("%s, %s, %s", gn, gid, id)
	print(titleText)
	outPrefix = sprintf("%s/%s_%s_%s", odir, gn, gid, id)
	
	dtDown = fread(ifileDownsample, header=TRUE)
	
	oldx = log10(dtDown[, rep(frac, n_trial)])
	oldy = lapply(1:nrow(dtDown), 
								function(nr){
									n_one = as.integer(dtDown[nr,detectionRate] * dtDown[nr, n_trial])
									n_zero = dtDown[nr, n_trial] - n_one
									rep(c(1, 0), c(n_one, n_zero))
									}
								) %>% unlist
	
	newx = log10(dt1[, fraction_intron])
	newy = as.numeric(dt1[, is_RS])
	
	rsDetectionProbs = calcRsDetectionProb(oldx, oldy, newx, newy, titleText, outPrefix=outPrefix)
	
	dt1 %<>% mutate(rsDetectionProb = rsDetectionProbs)
	dt1 %<>% mutate(rsDetectionReferenceGroup = gr)
	
	ofile = sprintf("%s/tidy_rsDetectionProb_%s_%s%s", odir, gn, id, ".txt")
	write.table(dt1, ofile, sep = "\t", quote=FALSE, row.names=FALSE)	
	
	# data.frame(test_stat = result[1], lower = result[2], greater = result[3]) %>%
	# 	kable %>% print

	# # Only one timepoint
	# titleText = sprintf("%s samples; %s, %s, %s", gr, gn, gid, id)
	# outPrefix = sprintf("%s/%s_%s_%s_%s", odir, gn, gid, id, gr)
	# newx = log10(dt1[group==gr, fraction_intron])
	# newy = as.numeric(dt1[group==gr, is_RS])
	# 
	# result = calcEmpiricalPval(oldx, oldy, newx, newy, titleText, outPrefix = outPrefix)
	# 
	# data.frame(test_stat = result[1], lower = result[2], greater = result[3]) %>%
	# 	kable %>% print

		
}

######################

dtRegion = fread(pathRegion, header=TRUE, sep="\t")

for(i in 1:nrow(dtRegion)){
	if(dtRegion[i, gene_name] %in% selGene){
		print(sprintf("%s, %s", dtRegion[i, gene_name], dtRegion[i, name] ))
		main00(i)
	}
}


#############################

sessionInfo()








































