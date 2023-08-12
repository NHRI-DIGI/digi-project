imiss <- read.table("process/miss_indQC.imiss", header = T, sep = "")
het <- read.table("process/het_indQC.het", header = T, sep ="")
het$Het_propo = (het$N.NM. - het$O.HOM.)/het$N.NM.
het_std<-sd(het$Het_propo)

# set missing call rate > 1%
fail_imiss<- imiss[imiss$F_MISS>0.01,c('FID','IID')]
fail_het<- het[het$Het_propo<mean(het$Het_propo)-3*het_std | het$Het_propo>mean(het$Het_propo)+3*het_std ,c('FID','IID')]

write.table(fail_imiss,'process/fail_imiss_ind.txt',quote=F,row.names=F, col.names=F)
write.table(fail_het,'process/fail_het_ind.txt',quote=F,row.names=F,col.names=F)
