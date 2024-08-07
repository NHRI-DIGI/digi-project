maf_freq <- read.table("process/MAF_check.frq", header =TRUE, as.is=T)
pdf("process/MAF_distribution.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF", xlim = c(0,0.5), breaks = 50)
dev.off()


