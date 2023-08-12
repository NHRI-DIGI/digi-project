#設定output的資料路徑,phenoarray名稱
#tool位置 /NovaSeq_128/Kenny/tool/bin/
#餵進來case/control 的txt file 在 ${pathway}Input/Laptop_input/
pathway='/NovaSeq_128/Kenny/CKD/test/'
phenoarray=CKD
tool_dir='/NovaSeq_128/Kenny/tool/bin/'
Outputbasedata='/NovaSeq_128/Kenny/CKD/test/Output/Basedata/'
OutputPRS='/NovaSeq_128/Kenny/CKD/test/Output/PRS/'
#============================以下不要動===============================
cd  ${pathway}Output/
mkdir PRS
cd  ${pathway}Output/PRS
mkdir $phenoarray

echo "0.001 0 0.001" > range_list
echo "0.005 0 0.005" >> range_list
echo "0.01 0 0.01" >> range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

# Clumping
${tool_dir}plink \
    --bfile /NovaSeq_128/Kenny/PRS/Input/Reference_genome/nodup_chr_merge_CHBS_1000g \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump ${Outputbasedata}"$phenoarray"_24/finalBase.txt \
    --clump-snp-field SNP \
    --clump-field P \
    --out ${OutputPRS}$phenoarray/clump #output: blahblah.clump
awk 'NR!=1{print $3}' ${OutputPRS}$phenoarray/clump.clumped >  ${OutputPRS}$phenoarray/clumped.snp

  # Prs & p-val thresholding
awk '{print $2,$9}' ${Outputbasedata}"$phenoarray"_24/finalBase.txt > ${OutputPRS}$phenoarray/SNP.pvalue
${tool_dir}plink \
    --bfile ${pathway}Output/Targetdata/finalTarget \
    --score ${Outputbasedata}"$phenoarray"_24/finalBase.txt 2 4 11 header \
    --q-score-range range_list ${OutputPRS}$phenoarray/SNP.pvalue \
    --extract ${OutputPRS}$phenoarray/clumped.snp \
    --out ${OutputPRS}$phenoarray/pval_th  #output: blahblah.profile

# Accounting for Population Stratification (PCA)
${tool_dir}plink \
    --bfile ${pathway}Output/Targetdata/finalTarget \
    --extract ${pathway}Output/Targetdata/LDr2.prune.in \
    --pca 20 \
    --out pcaOut  # output:pcaOut.eigenvec; Use elbow plot to see PC numbers
cp pcaOut.eigenv* ${pathway}Input/Laptop_input/


# copy final result to 'FinalResults' folder
cd ${pathway}Output/
mkdir FinalResults
cd FinalResults

rsync -a ${OutputPRS}$phenoarray .
cd ${pathway}Output/FinalResults/$phenoarray
rm clump* pval_th.log pval_th.nopred pval_th.nosex
cp -R ${pathway}Output/FinalResults/  ${pathway}Input/Laptop_input/

###
###Run bestPvalPRS.R to see model efficiency
###
