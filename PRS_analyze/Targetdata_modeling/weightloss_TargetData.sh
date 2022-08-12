# Build output route, and the pheno name -> phenoarray
# Merge VCF in merge_dir
# tools in ${pathway}Input/Script/bin/
### Important!!!! Make sure the chr = 1 or chr1
pathway='/staging2/reserve/flagship/u3121714/kenny/'
phenoarray=weightloss_6
merge_dir='/staging2/reserve/flagship/u3121714/kenny/Input/TargetData/hg38/weightloss/VCF_merge/'
#===========================================================
cd ${pathway}Output
mkdir TargetData

cd ${pathway}Output/TargetData
mkdir process

# Remove Duplicate SNPs
cat ${merge_dir}"$phenoarray".merge.vcf | egrep "^##" > process/"$phenoarray"Headers.txt
sed '/^##/d' ${merge_dir}"$phenoarray".merge.vcf  > process/"$phenoarray".vcf
head -n1 process/"$phenoarray".vcf > process/"$phenoarray"Header_2.txt

# CHR: X->23 Y->24 M->25 (avoid Invalid variant bp coordinate ERROR afterward)
awk -F'\t' -vOFS='\t' '/^#/ {next}; {gsub("X", "23", $1) ; gsub("Y", "24", $1); gsub("M", "25", $1); print}' process/"$phenoarray".vcf > process/noh_xym_"$phenoarray".vcf
cat process/"$phenoarray"Header_2.txt  process/noh_xym_"$phenoarray".vcf > process/xym_"$phenoarray".vcf

# Replace ID to chr_POS, Remove Duplicate SNPs
awk 'BEGIN{OFS="_";}  NR==1{print $2} NR >1{print $1,$2}' process/xym_"$phenoarray".vcf > process/chr_pos
# If original vcf file chr column = chr1... must do the code below.
# If not, don't do it.
sed -e 's/^chr//' -i process/chr_pos
#######################################
awk 'FNR==NR{a[NR]=$0;next} {$3=a[FNR]}1' OFS="\t" process/chr_pos  process/xym_"$phenoarray".vcf > process/noh_chr_xym_"$phenoarray".vcf
( sed 1q process/"$phenoarray"Header_2.txt; sed  1d process/noh_chr_xym_"$phenoarray".vcf ) > process/chr_xym_"$phenoarray".vcf
awk '{seen[$3]++; if(seen[$3]==1){ print}}' process/chr_xym_"$phenoarray".vcf > process/nodup_chr_xym_"$phenoarray".vcf

# Indiviual QC
${pathway}Input/Script/bin/plink --vcf process/nodup_chr_xym_"$phenoarray".vcf --make-bed --out process/indQC
${pathway}Input/Script/bin/plink --bfile process/indQC --missing --out process/miss_indQC
${pathway}Input/Script/bin/plink --bfile process/indQC --het --out process/het_indQC
Rscript ${pathway}Input/Script/indQC.r
cp process/miss_indQC.imiss ../../Input/Laptop_input/
mv ../../Input/Laptop_input/miss_indQC.imiss ../../Input/Laptop_input/target_miss_indQC.imiss
cp process/het_indQC.het ../../Input/Laptop_input/
mv ../../Input/Laptop_input/het_indQC.het ../../Input/Laptop_input/target_het_indQC.het
cp process/fail_imiss_ind.txt ../../Input/Laptop_input/
mv ../../Input/Laptop_input/fail_imiss_ind.txt ../../Input/Laptop_input/target_fail_imiss_ind.txt
cp process/fail_het_ind.txt ../../Input/Laptop_input/
mv ../../Input/Laptop_input/fail_het_ind.txt ../../Input/Laptop_input/target_fail_het_ind.txt
cat process/fail_het_ind.txt  process/fail_imiss_ind.txt | sort -k1 | uniq > fail_indQC.txt
${pathway}Input/Script/bin/plink --vcf process/nodup_chr_xym_"$phenoarray".vcf --remove fail_indQC.txt --make-bed  --out process/indQC_nodup_chr_xym_"$phenoarray"
${pathway}Input/Script/bin/plink --bfile process/indQC_nodup_chr_xym_"$phenoarray" --recode vcf-iid --out process/indQC_nodup_chr_xym_"$phenoarray"
sed -i '/^##/d' process/indQC_nodup_chr_xym_"$phenoarray".vcf

# Convert vcf file to bfile
${pathway}Input/Script/bin/plink --vcf process/indQC_nodup_chr_xym_"$phenoarray".vcf --make-bed --allow-no-sex --out process/indQC_nodup_chr_xym_"$phenoarray"

# SNP QC
${pathway}Input/Script/bin/plink \
    --bfile process/indQC_nodup_chr_xym_"$phenoarray" \
    --allow-no-sex \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --write-snplist \
    --make-just-fam \
    --out process/snpQC_target

# Export SNP QC report
preQC=$(< "process/chr_xym_"$phenoarray".vcf" wc -l )
nodup=$(< "process/nodup_chr_xym_"$phenoarray".vcf" wc -l )
sQC=$(< "process/snpQC_target.snplist" wc -l )
echo $preQC $nodup $sQC > process/snpQC_rowN.txt
Rscript ${pathway}Input/Script/snpQCreport_target.r

# Extract lowly correlated SNPs for PCA
${pathway}Input/Script/bin/plink \
    --bfile process/indQC_nodup_chr_xym_"$phenoarray" \
    --allow-no-sex \
    --extract process/snpQC_target.snplist \
    --indep-pairwise 200 50 0.2 \
    --out LDr2 # output: LDr2.prune.in

# Generate final target data file
${pathway}Input/Script/bin/plink \
    --bfile process/indQC_nodup_chr_xym_"$phenoarray" \
    --make-bed \
    --out finalTarget \
    --extract process/snpQC_target.snplist
