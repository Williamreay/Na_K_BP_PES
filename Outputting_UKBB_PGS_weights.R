#################################################

## Outputting PGS weights for the UKBB - full sample

## William Reay (2022)

#################################################

library(dplyr)
library(data.table)

setwd("~/Documents/Na_K_BP_Erin_Clare/")

## Select the best performing average 2 from two HCS cohorts

SBP_r2 <- fread("HCS_cohort_tuning_of_PRS_PES/SBP_HCS_tuning_results.txt")
SBP_r2$Mean <- (SBP_r2$SBP_score_med+SBP_r2$SBP_score_unmed)/2

## P < 0.01 best for SBP

DBP_r2 <- fread("HCS_cohort_tuning_of_PRS_PES/DBP_HCS_tuning_results.txt")
DBP_r2$Mean <- (DBP_r2$DBP_score_med+DBP_r2$DBP_score_unmed)/2

# P < 0.005 best for SBP

## Output weights with dbSNP SNPs for UKBB

SBP_sumstats <- fread("Genetic_data/Cleaned_GERA_sumstats_SBP.txt.gz")
DBP_sumstats <- fread("Genetic_data/Cleaned_GERA_sumstats_DBP.txt.gz")

## SBP formatting

SBP_weights <- fread("Genetic_data/HCS_PGS_or_PES_weights_CHR_BP/0.01_exMHC_SBP_HCS_weights.txt", header = T)
Merged_SBP <- merge(SBP_sumstats, SBP_weights, by ="SNP")
Merged_SBP$SBP_PGS_Pt_0_01 <- Merged_SBP$`0.01`

Merged_SBP <- Merged_SBP %>% select(dbSNP, EA.y, SBP_PGS_Pt_0_01)

write.table(Merged_SBP, file="UKBB_PGS/SBP_PGS_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## DBP formatting

DBP_weights <- fread("Genetic_data/HCS_PGS_or_PES_weights_CHR_BP/0.005_exMHC_DBP_HCS_weights.txt", header = T)
Merged_DBP <- merge(DBP_sumstats, DBP_weights, by ="SNP")
Merged_DBP$DBP_PGS_Pt_0_005 <- Merged_DBP$`0.005`

Merged_DBP <- Merged_DBP %>% select(dbSNP, EA.y, DBP_PGS_Pt_0_005)

write.table(Merged_DBP, file="UKBB_PGS/DBP_PGS_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

