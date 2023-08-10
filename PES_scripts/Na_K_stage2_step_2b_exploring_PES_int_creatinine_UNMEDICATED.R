#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2b - exploring transport PES x sodium interaction - UNMEDICATED ONLY

## William Reay (2022)

#################################################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggridges)
library(interactions)
library(sjPlot)
library(sjmisc)
library(viridis)
library(car)


Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")

Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)

Subcohort_1_dat <- Subcohort_1_dat %>% filter(BP_meds == "0")

## Read in PES

Na_K_PES <- fread("UKBB_PES/PES_Na_K_SBP_DBP.txt")

Colnames_PES <- make.unique(names(Na_K_PES))

colnames(Na_K_PES) <- Colnames_PES

Na_K_PES <- Na_K_PES %>% select(IID, contains("PES"))

Merged_PES <- merge(Na_K_PES, Subcohort_1_dat, by = "IID")

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")
Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_AVG))

## Test effect of scaled creatinine as well in model
Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_urine_millimolL))
Merged_PGS_PES$Scaled_CRET <- as.numeric(scale(Merged_PGS_PES$Creatinine_urine_micromolL))

CREAT <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                 Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                 PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Scaled_CRET, data = Merged_PGS_PES)

TE_CREAT <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
           Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Scaled_CRET + Transport_SBP_PES_AVG*Scaled_Na, data = Merged_PGS_PES)

Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                              Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                              PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                              PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                              PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                              Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                              Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                              PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                              PC10*Scaled_Na + PC11*Scaled_Na +
                              PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG + Transport_SBP_PES_AVG*Scaled_CRET + Scaled_Na*Scaled_CRET,
                            data = Merged_PGS_PES)

## BMI incl

BMI_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                              Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                              PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                              PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                              PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                              Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                              Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                              PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                              PC10*Scaled_Na + PC11*Scaled_Na +
                              PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG + Transport_SBP_PES_AVG*Scaled_CRET + Scaled_Na*Scaled_CRET +
                                Transport_SBP_PES_AVG*scale(BMI) + Scaled_Na*scale(BMI),
                            data = Merged_PGS_PES)




