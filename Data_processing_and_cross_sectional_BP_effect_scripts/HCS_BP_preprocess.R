##########################

## Processing the blood pressure data in the genotyped subset of the Hunter Community Study

## William Reay (2022)

########################

setwd("~/Documents/Na_K_BP_Erin_Clare/HCS_cohort_tuning_of_PRS_PES/")

library(dplyr)
library(data.table)
library(easyGgplot2)
library(readxl)
library(table1)

## Load HCS data 

HCS_pheno <- read.csv("HCS_phenotype_data.csv", header = T)

## Load covariates

HCS_cov <- fread("~/Documents/Hunter_cohort/Eigenvectors/HCS_PC1_PC5.tab.txt", header = F)

HCS_cov <- rename(HCS_cov, "ID"="V1", "PC1"="V3", "PC2"="V4", "PC3"="V5",
                  "PC4"="V6", "PC5"="V7")

Electoral_HCSID_conversion <- read_excel("~/cloudstor/Pneumonia_cytokine_lung_function/HCS_PES_profiles/Electoral_HCSID_conversion.xlsx")

## Derive EIDs from Hxxxx IDs

HCS_cov <- merge(HCS_cov, Electoral_HCSID_conversion, by = "ID")

HCS_cov <- rename(HCS_cov, "electoralId"="IID")

## Merge with pheno

Gen_merged_pheno <- merge(HCS_pheno, HCS_cov, by = "electoralId")

## Identify individuals with non-missing sex, age, SBP and DBP

Filt_gen_merged_pheno <- Gen_merged_pheno %>% filter(mbs != 'NA') %>% 
                                              filter(mbd != 'NA') %>%
                                              filter(bage != 'NA') %>%
                                              filter(sex.x != 'NA')

## Define medication indicator for antihypertensives and seperate cohorts

## ATC-Antihypertensives-wave 1, ACE AngiotensinII inhibitors-wave 1, BetaBlockerw1,
## ATC-Diuretics-wave 1, ATC-Beta blocking agents-wave 1, ATC-Ca channel blockers-wave 1
## ATC-Agents acting on the renin-angiotensin system-wave 1, Ca Channel Blockers-wave 1
## Diuretics-wave 1, 

Med_Filt_gen_pheno <- Filt_gen_merged_pheno %>% filter(C02w1 != 'NA') %>%
                                                filter(ACEangiotensinw1 != 'NA') %>%
                                                filter(BetaBlockerw1 != 'NA') %>%
                                                filter(C03w1 != 'NA') %>% filter(C07w1 != 'NA') %>%
                                                filter(C08w1 != 'NA') %>% filter(C09w1 != 'NA') %>%
                                                filter(CCBw1 != 'NA') %>% filter(Diureticw1 != 'NA')

## Define binary indicator of antihypertensive use at baseline from above df

Med_Filt_gen_pheno$BP_med <- ifelse((Med_Filt_gen_pheno$C02w1 == 1 | Med_Filt_gen_pheno$ACEangiotensinw == 1 | Med_Filt_gen_pheno$BetaBlockerw1 == 1 | Med_Filt_gen_pheno$C03w1 == 1 |
                                       Med_Filt_gen_pheno$C07w1 == 1 | Med_Filt_gen_pheno$C08w1 == 1 | Med_Filt_gen_pheno$C09w1 == 1 | Med_Filt_gen_pheno$CCBw1 == 1 | Med_Filt_gen_pheno$Diureticw1 == 1), "YES", "NO")

## Table 

table1( ~ bage + mbs + mbd + factor(sex.x) | factor(factor(BP_med)),
        data = Med_Filt_gen_pheno)

## Output cohort

write.table(Med_Filt_gen_pheno, file="HCS_BP_cohort.txt",
            sep = "\t", row.names = F, quote = F)

