############################################

## Pre-processing files for Erin and Clare

## Erin diet data and Na:K outlier individuals removed

## SBP/DBP

## William Reay (2022)

############################################

library(dplyr)
library(data.table)

setwd("~/Documents/Na_K_BP_Erin_Clare/")

## Pheno data with cleaned Na:K

Pheno_gen <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")

Diet <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Diet_physical_activity_UKBB_July_2021_Freeze.txt")

Diet_non_miss_baseline_quest <- Diet %>% filter(!is.na(`20080-0.0`))

## Merge

Pheno_full_diet <- merge(Pheno_gen, Diet_non_miss_baseline_quest, by = "eid")

table(Pheno_full_diet$`100020-0.0`)

## 37357 usual diet, 8282 not usual

table(Pheno_full_diet$`20085-0.0`)

## Of the 8282 not following their usual diet yesterday

## 3.1% were ill, 0.07% were fasting, 42% were away, and 54% for an unspecified other reason

write.table(Pheno_full_diet, file="Erin_diet_questionnaire_subcohort_UKBB.txt",
            sep = "\t", row.names = F, quote = F)
