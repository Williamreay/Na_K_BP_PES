#################################################

## Na/K and hypertension: interplay of diet and genetics

## eGFRcys benchmarking of UNa+

## William Reay (2022)

#################################################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggridges)
library(table1)
library(RColorBrewer)
library(effects)
library(nephro)

## Read in data generated from previous script

Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")
Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)

## Read in Cystatin C and filter for missingness

Biochem_df <- fread("../Biochem_immunology_data_UKBB_July_2021_Freeze.txt")
Cystatin_C_df <- Biochem_df %>% select(eid, `30720-0.0`)
Cystatin_C_df <- Cystatin_C_df %>% filter(!is.na(`30720-0.0`))
## Convert to mg/dl for formula calc
Cystatin_C_df$Cynstatin_C <- Cystatin_C_df$`30720-0.0`
write.table(Cystatin_C_df, file="Non_missing_Cystatin_C_df.txt", sep="\t", row.names = F, quote = F)

Subcohort_1_dat <- merge(Cystatin_C_df, Subcohort_1_dat, by="eid")
Subcohort_1_dat$Sex_coded <- ifelse(Subcohort_1_dat$Sex == "Male", 1,0)
## Using the nephro package v1.3
## Model two that accounts for age
## GFR is estimated with the CKD-EPI equation for cystatin C proposed by Inker et al., N Engl J Med 2012
## Also use the Stevens et al 2008 equation
Subcohort_1_dat$eGFRcys_CKD <- CKDEpi.cys(Subcohort_1_dat$Cynstatin_C, 
                                      Subcohort_1_dat$Sex_coded,
                                      Subcohort_1_dat$Age)

Subcohort_1_dat$Ethinicity <- 0

Subcohort_1_dat$eGFRcys_Stevens <- Stevens.cys2(Subcohort_1_dat$Cynstatin_C, 
                                          Subcohort_1_dat$Sex_coded,
                                          Subcohort_1_dat$Age,
                                          Subcohort_1_dat$Ethinicity)

eGFR_df <- Subcohort_1_dat

table1( ~ factor(Sex) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL + Cynstatin_C + eGFRcys_CKD + eGFRcys_Stevens | factor(`BP Medication`),
        data = eGFR_df)

## Winsorize at 15 ml/min/1.73m2 or 200 ml/min/1.73m2 = PMID: 34272381

eGFR_df_filt <- eGFR_df %>% filter(eGFRcys_CKD > 15 & eGFR_df$eGFRcys_CKD < 200)

## 41 oarticipants lost

Unmed_cohort <- eGFR_df_filt %>% filter(BP_meds == 0)

## Write new df

write.table(eGFR_df_filt, file="Winsorized_eGFR_adjusted_subcohort_1.txt",
            sep = "\t", row.names = F, quote = F)

## Test two eGFR methods

Compare <- Unmed_cohort %>% 
  ggplot(aes(eGFRcys_CKD,eGFRcys_Stevens,colour=Sex))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "firebrick") +
  theme_bw() +
  xlab("eGFRcys-CKD-Epi (mL/min/1.73m2)") +
  ylab("eGFRcys-Stevens et al. (mL/min/1.73m2)")

## Split into over and above 90, and stratify by sex

Unmed_cohort$Binary_CKD <- ifelse(Unmed_cohort$eGFRcys_CKD >= 90, "> 90", "< 90")
Unmed_cohort$Binary_Stevens <- ifelse(Unmed_cohort$eGFRcys_Stevens >= 90, "> 90", "< 90")

write.table(Unmed_cohort, file="eGFRcys_analyses/Unmed_binary_eGFR.txt",
            sep = "\t", row.names = F, quote = F)
