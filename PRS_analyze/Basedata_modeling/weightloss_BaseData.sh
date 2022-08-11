# Build output route, and the pheno name -> phenoarray
# Merge VCF in merge_dir
# tool in ${pathway}Input/Script/bin/
# case and control txt file in ${pathway}Input/Laptop_input/
### Important!!!! Make sure the chr = 1 or chr1
pathway='/staging2/reserve/flagship/u3121714/kenny/'
phenoarray=weightloss_12
merge_dir='/staging2/reserve/flagship/u3121714/kenny/Input/BaseData/hg38/weightloss/VCF_merge/'

#===========================================================
# Input case and contol txt
data=${pathway}Input/Laptop_input/weightloss_base.txt
# Classify control and case, 2 -> case，1 -> control; Due to manage the effect allele, so choose the oposite to set up
control=($(awk '$3 == 2 {printf $1" "}' $data))
case=($(awk '$3 == 1 {printf $1" "}' $data))

cd ${pathway}Output
mkdir Laptop_output
cd ${pathway}Output/Laptop_output
mkdir OutputImage  OutputImage/QCplot
cd ${pathway}Output
mkdir Basedata

cd ${pathway}Output/Basedata
mkdir $phenoarray process
mkdir $phenoarray/ef_al_process

# Make Heaader
cat ${merge_dir}"$phenoarray".merge.vcf | egrep "^##" > process/"$phenoarray"_Headers.txt
sed '/^##/d' ${merge_dir}"$phenoarray".merge.vcf > process/"$phenoarray".vcf
head -n1 process/"$phenoarray".vcf > process/"$phenoarray"_Header_2.txt

# CHR: X->23 Y->24 M->25 (avoid Invalid variant bp coordinate ERROR afterward)
awk -F'\t' -vOFS='\t' '/^#/ {next}; {gsub("X", "23", $1) ; gsub("Y", "24", $1); gsub("M", "25", $1); print}' process/"$phenoarray".vcf > process/noh_xym_"$phenoarray".vcf
cat process/"$phenoarray"_Header_2.txt process/noh_xym_"$phenoarray".vcf > process/xym_"$phenoarray".vcf

# Replace ID to chr_POS, Remove Duplicate SNPs
awk 'BEGIN{OFS="_";}  NR==1{print $2} NR >1{print $1,$2}' process/xym_"$phenoarray".vcf > process/chr_pos
# If original vcf file chr column = chr1... must do the code below.
# If not, don't do it.
sed -e 's/^chr//' -i process/chr_pos
#######################################
awk 'FNR==NR{a[NR]=$0;next} {$3=a[FNR]}1' OFS="\t" process/chr_pos  process/xym_"$phenoarray".vcf > process/noh_chr_xym_"$phenoarray".vcf
( sed 1q process/"$phenoarray"_Header_2.txt; sed 1d process/noh_chr_xym_"$phenoarray".vcf ) > process/chr_xym_"$phenoarray".vcf
awk '{seen[$3]++; if(seen[$3]==1){ print}}' chr_xym_"$phenoarray".vcf > process/nodup_chr_xym_"$phenoarray".vcf

# Indiviual QC
${pathway}Input/Script/bin/plink --vcf process/nodup_chr_xym_"$phenoarray".vcf --make-bed --out process/indQC
${pathway}Input/Script/bin/plink --bfile process/indQC --missing --out process/miss_indQC
${pathway}Input/Script/bin/plink --bfile process/indQC --het --out process/het_indQC
Rscript ${pathway}Input/Script/indQC.r
cp process/miss_indQC.imiss ../../Input/Laptop_input/
mv ../../Input/Laptop_input/miss_indQC.imiss ../../Input/Laptop_input/base_miss_indQC.imiss
cp process/het_indQC.het ../../Input/Laptop_input/
mv ../../Input/Laptop_input/het_indQC.het ../../Input/Laptop_input/base_het_indQC.het
cp process/fail_imiss_ind.txt ../../Input/Laptop_input/
mv ../../Input/Laptop_input/fail_imiss_ind.txt ../../Input/Laptop_input/base_fail_imiss_ind.txt
cp process/fail_het_ind.txt ../../Input/Laptop_input/
mv ../../Input/Laptop_input/fail_het_ind.txt ../../Input/Laptop_input/base_fail_het_ind.txt
cat process/fail_het_ind.txt  process/fail_imiss_ind.txt | sort -k1 | uniq > fail_indQC.txt
${pathway}Input/Script/bin/plink --vcf process/nodup_chr_xym_"$phenoarray".vcf --remove fail_indQC.txt --make-bed  --out process/indQC_nodup_chr_xym_"$phenoarray"
${pathway}Input/Script/bin/plink --bfile process/indQC_nodup_chr_xym_"$phenoarray" --recode vcf-iid --out process/indQC_nodup_chr_xym_"$phenoarray"
sed -i '/^##/d' process/indQC_nodup_chr_xym_"$phenoarray".vcf

# Add phenotype info (Every pheno output call "$phenoarray"_p)
${pathway}Input/Script/bin/plink --vcf process/indQC_nodup_chr_xym_"$phenoarray".vcf --pheno ${pathway}Input/Laptop_input/weightloss_base.txt --pheno-name $phenoarray --make-bed --out $phenoarray/"$phenoarray"_p

# SNP QC (Do only one pheno)
${pathway}Input/Script/bin/plink \
    --bfile $phenoarray/"$phenoarray"_p \
    --allow-no-sex \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --write-snplist \
    --make-just-fam \
    --out snpQC_base

# Export SNP QC report
preQC=$(< "process/chr_xym_"$phenoarray".vcf" wc -l )
nodup=$(< "process/nodup_chr_xym_"$phenoarray".vcf" wc -l )
sQC=$(< "snpQC_base.snplist" wc -l )
echo $preQC $nodup $sQC > process/snpQC_rowN.txt
Rscript ${pathway}Input/Script/snpQCreport_base.r

#Association Test
${pathway}Input/Script/bin/plink --bfile $phenoarray/"$phenoarray"_p --assoc --allow-no-sex --out $phenoarray/chi

# Remove Ambiguous SNPs
awk '!( ($4=="A" && $7=="T") || \
        ($4=="T" && $7=="A") || \
        ($ 4=="G" && $7=="C") || \
        ($4=="C" && $7=="G")) {print}' $phenoarray/chi.assoc > $phenoarray/ambRm_chi.assoc

# Effect allele Adjustment
cd $phenoarray
${pathway}Input/Script/bin/plink --bfile "$phenoarray"_p --recode vcf-iid --out ef_al_process/"$phenoarray"_p
cd ef_al_process

# Grep first row
sed -i '/^##/d' "$phenoarray"_p.vcf
head -n 1 "$phenoarray"_p.vcf > "$phenoarray"_p.txt

# Make control.vcf
declare -a myArray
myArray=(`cat "$phenoarray"_p.txt`)
declare -a num
num=()

for j in "${!control[@]}"
do
  for i in "${!myArray[@]}"
  do
    if [ "${control[$j]}" == "${myArray[$i]}" ]
    then
      num+=(`expr $i + 1`)
    fi
  done
done

num_complete=$(echo ${num[@]} | tr ' ' ',')

cut --fields="$num_complete" --complement ${pathway}Output/Basedata/$phenoarray/ef_al_process/"$phenoarray"_p.vcf > ${pathway}Output/Basedata/$phenoarray/ef_al_process/cont.vcf

# Make case.vcf
declare -a num_1
num_1=()

for j in "${!case[@]}"
do
  for i in "${!myArray[@]}"
  do
    if [ "${case[$j]}" == "${myArray[$i]}" ]
    then
      num_1+=(`expr $i + 1`)
    fi
  done
done

num_complete_1=$(echo ${num_1[@]} | tr ' ' ',')

cut --fields="$num_complete_1" --complement ${pathway}Output/Basedata/$phenoarray/ef_al_process/"$phenoarray"_p.vcf > ${pathway}Output/Basedata/$phenoarray/ef_al_process/case.vcf

# Make freq_case and freq_cont
${pathway}Input/Script/bin/plink --vcf cont.vcf --make-bed --out cont
${pathway}Input/Script/bin/plink --vcf case.vcf --make-bed --out case
${pathway}Input/Script/bin/plink --bfile cont --freq --out freq_cont
${pathway}Input/Script/bin/plink --bfile case --freq --out freq_case

# Extract 'POS' 'MAF'(of control&case), change these two header'MAF' into 'CONTROL_MAF'&'CASE_MAF'
awk '{print$1"\t"$2}' "$phenoarray"_p.vcf > pos_"$phenoarray"_p.txt
awk '{print$5}' freq_cont.frq > cont.txt
sed -i 's/MAF/CONTROL_MAF/g' cont.txt
awk '{print$5}' freq_case.frq > case.txt
sed -i 's/MAF/CASE_MAF/g' case.txt

# Merge 'POS','CASE_MAF','CONTROL_MAF'
paste pos_"$phenoarray"_p.txt case.txt cont.txt > ef.txt

# Find out which POS is control>case
awk '$4>$3{print$2}' ef.txt > efp_ct.txt
perl ${pathway}Input/Script/ef_ad.pl #otput: ef_ambRm_chi.assoc
cd ../../

# Add log(OR) column (BaseData)
awk '{s=(NR==1)?"logOR":log($10);$0=$0 OFS s}1' $phenoarray/ef_ambRm_chi.assoc > $phenoarray/finalBase.txt
cp -r $phenoarray/ef_ambRm_chi.assoc ../finalBase_OR.txt
sed -i 's|-inf|NA|g' $phenoarray/finalBase.txt
