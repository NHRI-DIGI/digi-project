# Imputation courses using Beagle and GLIMPSE
---
## 專題實作講義(PPT file)

* [01. Imputation](https://drive.google.com/file/d/1jfO8SPqhp3oyhaYlAbzl6cAkb2Sj22zR/view?usp=sharing)
* Beagle code :
java -jar /NovaSeq_128/Digit_2024/tool/bin/beagle.27Feb25.75f.jar \
  gt=/NovaSeq_128/Digit_2024/[Name]/BGE_illumina_chr22.vcf.gz \
  ref=/NovaSeq_128/Digit_2024/[Name]/TWReference_2500_MAF00002_chr22_rmmulti_phasing.vcf.gz \
  map=/NovaSeq_128/Digit_2024/[Name]/chr22.GRCh38.map \
  gp=true \
  out=/NovaSeq_128/Digit_2024/[Name]/[Output_name] ex: Test_chr22
* Environment variables configuration
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/chunk/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/phase/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/ligate/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/bcftools-1.17/bin:$PATH’ >> ~/.bashrc


* Chunking a chromosome
INDIR=/NovaSeq_128/Digit_2024/[Name]
OUTDIR=/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE

for chr in 22; do
 /home/lyf/GLIMPSE-1.1.1/chunk/bin/GLIMPSE_chunk \
    --input ${INDIR}/TWReference_2500_MAF00002_chr${chr}_rmmulti_phasing.vcf.gz \
    --region chr${chr} \
    --window-size 2000000 \
    --buffer-size 200000 \
    --output ${OUTDIR}/chunks.chr${chr}.Glimpse.txt \
    --log ${OUTDIR}/chunk.chr${chr}.log
done

* Compute genotype likelihoods (GLs) for individual samples at specified variant positions
INDIR=/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/lcWGS_demo_bam
OUTDIR=/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/GT_by_Ref
VCFDIR=/NovaSeq_128/Digit_2024/[Name]
REFDIR=/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/reference_genome

for BAM in ${INDIR}/*.bam; do
    SAMPLE=$(basename ${BAM} .bam)
    for i in 22; do
        mkdir -p ${OUTDIR}/chr${i}
        VCF=${VCFDIR}/TWReference_2500_MAF00002_chr${i}_rmmulti_phasing.sites.vcf.gz
        TSV=${VCFDIR}/TWReference_2500_MAF00002_chr${i}_rmmulti_phasing.sites.tsv.gz
        REF=${REFDIR}/hs38DH.chr${i}.fa
        OUT=${OUTDIR}/chr${i}/${SAMPLE}.chr${i}.vcf.gz
        bcftools mpileup -f ${REF} -I -E -a FORMAT/DP -T ${VCF} -r chr${i} ${BAM} -Ou | \
        bcftools call -Aim -C alleles -T ${TSV} -Oz -o ${OUT}
        bcftools index -f ${OUT}
    done
done

* [02. For beginner of Linux](https://drive.google.com/file/d/1Rv-wirTVwoVo0o7aL1jU59XEmkkgQLfz/view?usp=sharing)
---

