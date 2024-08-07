#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location
library("qqman")
results_log <- read.table("/NovaSeq_128/Kenny/PRS/hypertension/Output/BaseData/completeqc.assoc.logistic", head=TRUE)
jpeg("/NovaSeq_128/Kenny/PRS/hypertension/Output/BaseData/Logistic_manhattan.jpeg")
manhattan(results_log, chr="CHR",bp="BP",p="P",snp="SNP", ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), main = "Manhattan plot: logistic")
dev.off()
