n<-read.table("process/snpQC_rowN.txt", header = F, sep = " ")

# preQC: as.numeric(n[1])
# nodup: as.numeric(n[2])
# sQC: as.numeric(n[3])

Remain_SNP_number<-c(as.numeric(n[2])-1,as.numeric(n[3]))
Removed_SNP_number<-c(as.numeric(n[1])-as.numeric(n[2]), as.numeric(n[2])-as.numeric(n[3]))
Steps<- c('Remove_duplicated_SNPs',paste0('Standard_GWAS_QC: ',"\n",' -- Minor allele frequency (MAF) >0.01'," \n",' -- HWE deviation: P > 0.00001'))



reportb<-data.frame(Steps,Remain_SNP_number,Removed_SNP_number)


write.table(reportb,'../../Input/Laptop_input/Base_snpQCreport.txt',quote=F,row.names=F, col.names=T,sep="\t")
