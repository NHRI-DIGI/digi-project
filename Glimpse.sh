# Environment variables configuration
echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/chunk/bin:$PATH' >> ~/.bashrc
echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/phase/bin:$PATH' >> ~/.bashrc
echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/ligate/bin:$PATH' >> ~/.bashrc
echo 'export PATH=/usr/local/genome/bcftools-1.17/bin:$PATH' >> ~/.bashrc


# Chunking a chromosome

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

# Compute genotype likelihoods (GLs) for individual samples at specified variant positions

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

# Merge GLs across multiple individuals
BASE=/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/GT_by_Ref

for i in 22; do
    DIR=${BASE}/chr${i}
    OUT=${DIR}/Merge.chr${i}.vcf.gz
    bcftools merge \
    -Oz \
    -o ${OUT} \
    ${DIR}/*.vcf.gz

    bcftools index -f ${OUT}
done

# Impute and phase genotypes across the entire chromosome
BASE_VCF_DIR="/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/GT_by_Ref"
REF_DIR="/NovaSeq_128/Digit_2024/[Name]"
MAP_DIR="/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/GLIMPSE_map"
CHUNK_DIR="/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE“

for i in 22; do
    CHR="chr${i}"; OUTPUT_DIR="${BASE_VCF_DIR}/${CHR}/GLIMPSE"; mkdir -p "${OUTPUT_DIR}"
    while IFS=$'\t' read -r ID CHR2 IRG ORG WINDOW_CM WINDOW_BP COUNT1 COUNT2 || [[ -n "$ID" ]]; do
        [[ -z "$ID" ]] && continue
        OUT="${OUTPUT_DIR}/${ID}_${IRG//:/_}.vcf.gz"
        /home/lyf/GLIMPSE-1.1.1/phase/bin/GLIMPSE_phase \
            --input "${BASE_VCF_DIR}/${CHR}/Merge.${CHR}.vcf.gz" --reference "${REF_DIR}/TWReference_2500_MAF00002_${CHR}_rmmulti_phasing.vcf.gz" \
            --map "${MAP_DIR}/${CHR}.b38.gmap.fixed.gz" --input-region "${IRG}" --output-region "${ORG}" \
            --output "${OUT}" --seed 15052011 --burnin 2 --main 2 --thread 25
        bcftools index -f "${OUT}"
    done < "${CHUNK_DIR}/chunks.${CHR}.Glimpse.txt"
done
# Ligate imputed chromosomal chunks into a single file

BASE="/NovaSeq_128/Digit_2024/[Name]/For_GLIMPSE/GT_by_Ref"
for i in 22; do
    CHR="chr${i}"
    DIR="${BASE}/${CHR}/GLIMPSE"

    ls ${DIR}/*.vcf.gz > ${DIR}/list.${CHR}.txt

    /home/lyf/GLIMPSE-1.1.1/ligate/bin/GLIMPSE_ligate \
        --input ${DIR}/list.${CHR}.txt \
        --output ${DIR}/Merge.${CHR}.Impute.vcf.gz

    bcftools index -f ${DIR}/Merge.${CHR}.Impute.vcf.gz
done
