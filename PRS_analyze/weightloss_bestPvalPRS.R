# 註解掉的「程式碼」都可無視
library(data.table)
library(magrittr)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(stats)
library(matrixStats)

setwd('/Users/kenny/Desktop/Laptop_input')

urPhen <- "weightloss" # Type phenotype here
phenotype <- read.table("/Users/kenny/Desktop/Laptop_input/weightloss_target.txt", header=T) %>%
  mutate(weightloss = weightloss - 1)

# Scree plot -> choose principle componets #################
tiff(filename = "/Users/kenny/Desktop/Laptop_output/screeplot.tiff", height = 6, width = 8, units = 'in', res = 600)
  eigenval <- read.table("/Users/kenny/Desktop/Laptop_input/pcaOut.eigenval")$V1
  plot(
    x = seq(1:length(eigenval)), y = eigenval,
    type = "o",
    xlab = "Principle Component", ylab = "Variance")
dev.off()
#####################
pcs <- fread("/Users/kenny/Desktop/Laptop_input/pcaOut.eigenvec", header = F) %>%
    setnames(., colnames(.), c("FID", "IID", paste0("PC",1:12))) # choose all PCs
pheno <- phenotype %>%
  merge(., pcs)
colnames(pheno)[3] <- 'mypheno'

#Choose 11 PC due to the scree plot#####
pheno <- pheno[,c(1:14)] # 11 PC

# linear model 改用 logistic regression model
# 人口組成（PCs）和pheno相關性
# pheno[,!colnames(pheno)%in%c("FID","IID")]
# null.r2 <- summary(lm(mypheno~., data=pheno[,!colnames(pheno)%in%c("FID","IID")]))$r.squared
prs.result <- NULL
p.threshold <- c(0.001,0.005,0.01,0.05,0.1,0.2,0.3,0.4,0.5)
for(i in p.threshold){
  pheno.prs <- paste0('/Users/kenny/Desktop/Laptop_input/FinalResults/',urPhen,"/pval_th.", i, ".profile") %>%
    fread(.) %>%
    .[,c("FID", "IID", "SCORE")] %>%
    merge(., pheno, by = c("FID", "IID"))
  pheno.prs <- data.frame(pheno.prs)
  # prs分數跟pheno相關性
  #model<-glm(mypheno~.,family=binomial(link="logit"), data=pheno.prs[,!colnames(pheno)%in%c("FID","IID")])%>% summary
  model <- glm(mypheno ~ ., family = binomial(link = "logit"), data = pheno.prs[,!colnames(pheno)%in%c("FID","IID")])
  nullmod <- glm(mypheno ~ 1, family = binomial(link = "logit"), data = pheno.prs[,!colnames(pheno)%in%c("FID","IID")])

  1 - logLik(model) / logLik(nullmod)
  model.r2 <- 1 - logLik(model) / logLik(nullmod)
  #prs.r2 <- model.r2-null.r2
  #prs.coef <- model$coeff["SCORE",]
  prs.coef <- summary(model)$coeff["SCORE",]
  prs.result %<>% rbind(.,
                        data.frame(Threshold = i, R2 = prs.r2,
                                   P = as.numeric(prs.coef[4]),
                                   BETA = as.numeric(prs.coef[1]),
                                   SE = as.numeric(prs.coef[2])))
}
print(prs.result[which.max(prs.result$R2),])
bestCutOff <- prs.result[which.max(prs.result[-1, 2]) + 1,]$Threshold # 不考慮0.001
prs.result

write.table(paste(urPhen,bestCutOff), file = paste0("/Users/kenny/Desktop/Laptop_output/",urPhen,'_BestPval.txt'), sep = "\t", quote = F,row.names= F, col.names = F)

# prs density plot (雙峰圖)
res <- read.table(paste0('/Users/kenny/Desktop/Laptop_input/FinalResults/',urPhen,"/pval_th.",bestCutOff,".profile"), header = T)
p <- phenotype
p2 <- select(p,FID,IID,urPhen)
colnames(p2)[3] <-'mypheno'
da <- merge(res, p2, by = 'FID')
da$Standardized_PRS_Score <- scale(da$SCORE)

da1 <- replace(da, da == 0,'Control')
da2 <- replace(da1, da1 == 1,'Case')
da2$urPhen <- factor(da2$mypheno)

gd <- ggdensity(da2, x = "Standardized_PRS_Score",
   add = "mean", rug = TRUE,
   color = "urPhen", fill = "urPhen",
   palette = c("#FC4E07", "#0073C2FF")) + theme(legend.title = element_blank()) + ggtitle(urPhen)
ggsave(plot = gd, filename = paste0("/Users/kenny/Desktop/Laptop_output/OutputImage/",urPhen,"_DensityPlot.tiff"), height = 4, width = 6)

# Export indQC plot
imiss <- read.table("/Users/kenny/Desktop/Laptop_input/base_miss_indQC.imiss", header = T, sep = "")
het <- read.table("/Users/kenny/Desktop/Laptop_input/base_het_indQC.het", header = T, sep = "")

# Visulize individals fail imiss and het
het$Het_propo <- (het$N.NM. - het$O.HOM.)/het$N.NM.
imiss_het <- merge(het,imiss,by = "FID")

# Make fail imiss and het plots
Addcolors <- densCols(imiss_het$F_MISS,imiss_het$Het_propo)

tiff(filename = '/Users/kenny/Desktop/Laptop_output/OutputImage/QCplot/indQCplot_base.tiff',height = 6, width = 8,units = 'in', res = 600)

indQCplot <- plot(imiss_het$F_MISS,imiss_het$Het_propo, col = Addcolors, pch = 20,
     xlim = c(0, 0.06),
     ylim = c(0.18, 0.212),
     xlab = "Missing rate", ylab = "Heterozygosity rate", axes = FALSE, cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)
axis(2, tick = T) #2, left. t, tickmarks
axis(1, tick = T) #1, below
legend("bottomright", legend = "Genotypes_missing_het", pch = 16, bty = "n", col = Addcolors)

abline(h = mean(imiss_het$Het_propo) -
         (3 * sd(imiss_het$Het_propo)), col = "RED",lty = 2)
abline(h = mean(imiss_het$Het_propo)
       +(3 * sd(imiss_het$Het_propo)), col = "RED",lty = 2)

# Missing Data Thresholds (Vertical Line)
abline(v = 0.03, col = "BLUE", lty = 2) # THRESHOLD=missing_rate
title("Distribution of missing and heterozygosity rate")

m1 <- mean(imiss_het$Het_propo) - (3 * sd(imiss_het$Het_propo))
m2 <- mean(imiss_het$Het_propo) + (3 * sd(imiss_het$Het_propo))

if (length(imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'FID']) != 0)
{
text(x = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'F_MISS' ], y = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'Het_propo'],labels = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'FID'],cex = 0.9, font = 2, pos = 4)
}
dev.off()

imiss <- read.table("/Users/kenny/Desktop/Laptop_input/target_miss_indQC.imiss", header = T, sep = "")
het <- read.table("/Users/kenny/Desktop/Laptop_input/target_het_indQC.het", header = T, sep ="")

# Visulize individals fail imiss and het
het$Het_propo = (het$N.NM. - het$O.HOM.) / het$N.NM.
imiss_het = merge(het,imiss,by = "FID")
# Make a plot
Addcolors = densCols(imiss_het$F_MISS,imiss_het$Het_propo)

tiff(filename = '/Users/kenny/Desktop/Laptop_output/OutputImage/QCplot/indQCplot_target.tiff', height = 6, width = 8, units = 'in', res = 600)

indQCplot <- plot(imiss_het$F_MISS,imiss_het$Het_propo, col = Addcolors, pch = 20,
     xlim = c(0,0.06),
     ylim = c(0.18,0.212),
     xlab = "Missing rate", ylab = "Heterozygosity rate", axes = FALSE, cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.3)
axis(2,tick = T) #2, left. t, tickmarks
axis(1,tick = T) #1, below
legend("bottomright", legend = "Genotypes_missing_het", pch = 16, bty = "n", col = Addcolors)
abline(h = mean(imiss_het$Het_propo) -
         (3*sd(imiss_het$Het_propo)), col = "RED", lty = 2)
abline(h = mean(imiss_het$Het_propo)
       + (3*sd(imiss_het$Het_propo)), col = "RED", lty = 2)
##Missing Data Thresholds (Vertical Line)
abline(v = 0.03, col = "BLUE", lty = 2) # THRESHOLD=missing_rate
title("Distribution of missing and heterozygosity rate")

m1 <- mean(imiss_het$Het_propo) - (3*sd(imiss_het$Het_propo))
m2 <- mean(imiss_het$Het_propo) + (3*sd(imiss_het$Het_propo))

if (length(imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'FID']) != 0)
{
text(x = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1,'F_MISS'], y = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1, 'Het_propo'],labels = imiss_het[imiss_het$F_MISS > 0.03 | imiss_het$Het_propo > m2 | imiss_het$Het_propo < m1, 'FID'], cex = 0.9, font = 2, pos = 4)
}
dev.off()
