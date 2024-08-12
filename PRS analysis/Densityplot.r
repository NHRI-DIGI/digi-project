library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ggpubr)

setwd('/Users/kenny/Desktop')

phenotype <- fread("/Users/kenny/Desktop/HBOC_target_phenotype.txt") %>%
  select("V2", "V3") %>%
  rename(ID = V2, pheno = V3)

score <- fread("/Users/kenny/Desktop/PRSice.all.score", header = TRUE) %>%
  select(2, 3) %>%
  rename(ID = 1, score = 2)

score_pheno <- phenotype %>%
  left_join(score, by = "ID") %>%
  select("score", "pheno")

score_pheno$pheno <- str_replace(score_pheno$pheno, "1", "CONTROL")
score_pheno$pheno <- str_replace(score_pheno$pheno, "2", "CASE")

densityplot <- ggdensity(score_pheno, x = "score", color = "pheno", fill = "pheno", add = "mean", rug = TRUE, palette = c("#FC4E07", "#0073C2FF"))+theme(legend.title=element_blank())
ggsave(plot=densityplot, filename=paste0("/Users/kenny/Desktop/HBOC_densityplot.tiff"), height = 4, width = 6)
