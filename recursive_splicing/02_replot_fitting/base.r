
###############

library(data.table)
library(ggplot2)
library(dplyr)
library(dtplyr)
library(magrittr)

#############

if(! file.exists(dirname(ofile))) dir.create(dirname(ofile))


#################
# Load data
load(path_rdata)

########################

to_string = dtpval[, sprintf("%s\np=%.3e", group, pval.raw_new)]
names(to_string) = dtpval[, group]

# Color by p-value
dt = merge(dt, dtpval[, .(group, isRs)], by = "group")

## Reorder by p-value
# neworder = dtpval %>% arrange(pval.raw_new) %>% select(group) %>% unlist
# dt[, group := factor(group, levels = neworder)]

title = dtpval$id[1]

# labeller =  dtpval[, sprintf("%s (p=%.3e)", group, pval.raw)]
# names(labeller) = dtpval[,group]
# labeller = as_labeller(labeller)

xmin = dt[, min(bp)]

pdfname = ofile
pdf(pdfname, paper="a4", height = 11, width = 8)

ymax = dt[, max(depth)]

g = ggplot(dt, aes(x = bp, y = depth))
g = g + geom_histogram(stat="identity",  aes(fill = isRs), width=dt$bp[2]-dt$bp[1])
g = g + scale_fill_manual(values = c(`TRUE` = "red", `FALSE`="grey"))
g = g + geom_line(aes(x = bp, y = y.rss.pred), size = 0.2)
g = g + facet_wrap(~group, ncol = 6, strip.position = "left", labeller=as_labeller(to_string))
# g = g + geom_text(data = dtpval, aes(label = sprintf("p=%.3e", pval.raw)), x = xmin, y = ymax, hjust=0, vjust=1)
g = g + ggtitle(title)
g = g + scale_y_continuous(position = "right")
g = g + theme_bw()
# g = g + theme(aspect.ratio = 1/10)
g = g + theme(legend.position = "none")
g = g + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour="grey")
)
g = g + theme(strip.background = element_blank())
g = g + theme(strip.text = element_text(size = 2))
g = g + theme(strip.text.y = element_text(angle = 180))
g = g + theme(axis.text.x = element_text(size = 2))
g = g + theme(axis.text.y = element_text(size = 2))
g = g + theme(panel.spacing.y = unit(0, "lines"))
g = g + theme(axis.ticks = element_line(size = 0.1))
g = g + theme(axis.ticks.length = unit(0.1, "lines"))
g = g + theme(aspect.ratio = 1/2)

plot(g)


g = g + facet_wrap(~group, ncol = 6, strip.position = "left", labeller=as_labeller(to_string), scales = "free_y")
plot(g)
# ymax = dt[, max(depth.sm)]

# g = ggplot(dt, aes(x = bp, y = depth.sm))
# g = g + geom_histogram(stat="identity",  aes(fill = group), width=dt$bp[2]-dt$bp[1])
# g = g + scale_fill_manual(breaks = groups, values = color_labels)
# g = g + geom_line(aes(x = bp, y = y.sm.rss.pred))
# g = g + facet_wrap(~group, ncol = 10, strip.position = "left")
# g = g + geom_text(data = dtpval, aes(label = sprintf("p=%.3e", pval.sm)), x = xmin, y = ymax, hjust=0, vjust=1)
# g = g + ggtitle(title)
# g = g + scale_y_continuous(position = "right")
# g = g + theme_bw()
# # g = g + theme(aspect.ratio = 1/10)
# g = g + theme(legend.position = "none")
# g = g + theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank()
# )

# print(g)


dev.off()


