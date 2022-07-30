mkdir /staging/reserve/flagship/u6302720/kenny/Output/1000g
cd /staging/reserve/flagship/u6302720/kenny/Output/1000g
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr"
suffix=".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" ;
for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130501_sample_info/20130501_g1k.ped
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

for chr in {1..22}; do
    ../../Input/Script/bin/bcftools norm -m-any ALL.chr"${chr}".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | \
      ../../Input/Script/bin/bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        ../../Input/Script/bin/bcftools norm -Ob --rm-dup both \
          > ALL.chr"${chr}".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.bcf ;

    ../../Input/Script/bin/bcftools index ALL.chr"${chr}".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.bcf ;
done

for chr in {1..22}; do
    ../../Input/Script/bin/plink --noweb \
      --bcf ALL.chr"${chr}".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased ;
done



for i in {1..22}
do
../../Input/Script/bin/bcftools view -S ../../Input/CHBCHS.txt ALL.chr"$i".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz --force-samples > CHBS_ALL.chr"$i".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf
sed -i '/^##/d' CHBS_ALL.chr"$i".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf
done


grep '^#' CHBS_ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf > merge_CHBS_1000g.vcf
grep -v '^#' -h CHBS_ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr2.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr3.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr4.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr5.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr6.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr7.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr8.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr9.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr10.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr11.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr12.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr13.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr15.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr18.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr19.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr20.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
CHBS_ALL.chr22.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf \
   >> merge_CHBS_1000g.vcf


# Replace ID to chr_POS, Remove Duplicate SNPs
head -n1 merge_CHBS_1000g.vcf > CHBShead2
awk 'BEGIN{OFS="_";}  NR==1{print $2} NR >1{print $1,$2}' merge_CHBS_1000g.vcf > chr_pos_1000Gmerge
awk 'FNR==NR{a[NR]=$0;next} {$3=a[FNR]}1' OFS="\t" chr_pos_1000Gmerge  merge_CHBS_1000g.vcf > noh_merge_CHBS_1000g.vcf
( sed 1q CHBShead2; sed  1d noh_merge_CHBS_1000g.vcf ) > chr_merge_CHBS_1000g.vcf
(head -n 1 chr_merge_CHBS_1000g.vcf && tail -n +2 chr_merge_CHBS_1000g.vcf| sort -u  -k3,3 ) > ../../Input/nodup_chr_merge_CHBS_1000g.vcf
../../Input/Script/bin/plink --vcf ../../Input/nodup_chr_merge_CHBS_1000g.vcf --make-bed --out ../../Input/nodup_chr_merge_CHBS_1000g
