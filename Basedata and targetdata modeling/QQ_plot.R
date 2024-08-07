#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location
library("qqman")
results_log <- read.table("/NovaSeq_128/Kenny/PRS/hypertension/Output/BaseData/completeqc.assoc.logistic", head=TRUE)
jpeg("/NovaSeq_128/Kenny/PRS/hypertension/Output/BaseData/QQ-Plot_logistic.jpeg")
qq(results_log$P, main = "Q-Q plot of GWAS p-values : log", xlim = c(0, 7), ylim = c(0, 12), pch = 18, col = "blue4", cex = 1.5, las = 1)
dev.off()
