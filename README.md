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
  out=/NovaSeq_128/Digit_2024/[Name]/[Output_name] ex: Test_chr22* Environment variables configuration
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/chunk/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/phase/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/GLIMPSE-1.1.1/ligate/bin:$PATH' >> ~/.bashrc
* echo 'export PATH=/usr/local/genome/bcftools-1.17/bin:$PATH’ >> ~/.bashrc

* [02. For beginner of Linux](https://drive.google.com/file/d/1Rv-wirTVwoVo0o7aL1jU59XEmkkgQLfz/view?usp=sharing)
---

