hwe<-read.table (file="process/hwe_indQC.hwe", header=TRUE)
pdf("process/histhwe.pdf")
hist(hwe[,9], xlab = "P_value", ylab = "SNPs_count", breaks = 100, main="Histogram HWE")
dev.off()

hwe_zoom<-read.table (file="process/select_hwe_indQC.hwe", header=TRUE)
pdf("process/histhwe_below_theshold.pdf")
hist(hwe_zoom[,9], xlab = "P_value", ylab = "SNPs_count", breaks = 50, main="Histogram HWE: strongly deviating SNPs only")
dev.off()
