#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2b - exploring transport PES x sodium interaction - UNMEDICATED ONLY

## Stratify by eGFRcys

## William Reay (2023)

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
library(ggsci)

Subcohort_1_dat <- fread("eGFRcys_analyses/Unmed_binary_eGFR.txt")

Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)
## Read in PES

Na_K_PES <- fread("UKBB_PES/PES_Na_K_SBP_DBP.txt")

Colnames_PES <- make.unique(names(Na_K_PES))

colnames(Na_K_PES) <- Colnames_PES

Na_K_PES <- Na_K_PES %>% select(IID, contains("PES"))

Merged_PES <- merge(Na_K_PES, Subcohort_1_dat, by = "IID")

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")

## Test interaction with eGFRCys

Merged_PGS_PES$Scaled_eGFRcys_1 <- as.numeric(scale(Merged_PGS_PES$eGFRcys_CKD))
Merged_PGS_PES$Scaled_eGFRcys_2 <- as.numeric(scale(Merged_PGS_PES$eGFRcys_Stevens))

eGFR1_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*scale(Transport_SBP_PES_AVG) + Age*scale(Transport_SBP_PES_AVG) + Age2*scale(Transport_SBP_PES_AVG) + Assesment_centre_month*scale(Transport_SBP_PES_AVG) + 
                                       Assesment_centre*scale(Transport_SBP_PES_AVG) + PC1*scale(Transport_SBP_PES_AVG) + PC2*scale(Transport_SBP_PES_AVG) + PC3*scale(Transport_SBP_PES_AVG) + PC4*scale(Transport_SBP_PES_AVG) + PC5*scale(Transport_SBP_PES_AVG) + PC6*scale(Transport_SBP_PES_AVG) +
                                       PC7*scale(Transport_SBP_PES_AVG) + PC8*scale(Transport_SBP_PES_AVG) + PC9*scale(Transport_SBP_PES_AVG) +
                                       PC10*scale(Transport_SBP_PES_AVG) + PC11*scale(Transport_SBP_PES_AVG) +
                                       PC12*scale(Transport_SBP_PES_AVG) + PC13*scale(Transport_SBP_PES_AVG) + PC14*scale(Transport_SBP_PES_AVG) + PC15*scale(Transport_SBP_PES_AVG) + PC16*scale(Transport_SBP_PES_AVG) + PC17*scale(Transport_SBP_PES_AVG) + PC18*scale(Transport_SBP_PES_AVG) + PC19*scale(Transport_SBP_PES_AVG) + PC20*scale(Transport_SBP_PES_AVG) + 
                                       Sex*Scaled_eGFRcys_1 + Age*Scaled_eGFRcys_1 + Age2*Scaled_eGFRcys_1 + Assesment_centre_month*Scaled_eGFRcys_1 + 
                                       Assesment_centre*Scaled_eGFRcys_1 + PC1*Scaled_eGFRcys_1 + PC2*Scaled_eGFRcys_1 + PC3*Scaled_eGFRcys_1 + PC4*Scaled_eGFRcys_1 + PC5*Scaled_eGFRcys_1 + PC6*Scaled_eGFRcys_1 +
                                       PC7*Scaled_eGFRcys_1 + PC8*Scaled_eGFRcys_1 + PC9*Scaled_eGFRcys_1 +
                                       PC10*Scaled_eGFRcys_1 + PC11*Scaled_eGFRcys_1 +
                                       PC12*Scaled_eGFRcys_1 + PC13*Scaled_eGFRcys_1 + PC14*Scaled_eGFRcys_1 + PC15*Scaled_eGFRcys_1 + PC16*Scaled_eGFRcys_1 + PC17*Scaled_eGFRcys_1 + PC18*Scaled_eGFRcys_1 + PC19*Scaled_eGFRcys_1 + PC20*Scaled_eGFRcys_1 +Scaled_eGFRcys_1*scale(Transport_SBP_PES_AVG),
                                     data = Merged_PGS_PES)

eGFR2_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*scale(Transport_SBP_PES_AVG) + Age*scale(Transport_SBP_PES_AVG) + Age2*scale(Transport_SBP_PES_AVG) + Assesment_centre_month*scale(Transport_SBP_PES_AVG) + 
                                    Assesment_centre*scale(Transport_SBP_PES_AVG) + PC1*scale(Transport_SBP_PES_AVG) + PC2*scale(Transport_SBP_PES_AVG) + PC3*scale(Transport_SBP_PES_AVG) + PC4*scale(Transport_SBP_PES_AVG) + PC5*scale(Transport_SBP_PES_AVG) + PC6*scale(Transport_SBP_PES_AVG) +
                                    PC7*scale(Transport_SBP_PES_AVG) + PC8*scale(Transport_SBP_PES_AVG) + PC9*scale(Transport_SBP_PES_AVG) +
                                    PC10*scale(Transport_SBP_PES_AVG) + PC11*scale(Transport_SBP_PES_AVG) +
                                    PC12*scale(Transport_SBP_PES_AVG) + PC13*scale(Transport_SBP_PES_AVG) + PC14*scale(Transport_SBP_PES_AVG) + PC15*scale(Transport_SBP_PES_AVG) + PC16*scale(Transport_SBP_PES_AVG) + PC17*scale(Transport_SBP_PES_AVG) + PC18*scale(Transport_SBP_PES_AVG) + PC19*scale(Transport_SBP_PES_AVG) + PC20*scale(Transport_SBP_PES_AVG) + 
                                    Sex*Scaled_eGFRcys_2 + Age*Scaled_eGFRcys_2 + Age2*Scaled_eGFRcys_2 + Assesment_centre_month*Scaled_eGFRcys_2 + 
                                    Assesment_centre*Scaled_eGFRcys_2 + PC1*Scaled_eGFRcys_2 + PC2*Scaled_eGFRcys_2 + PC3*Scaled_eGFRcys_2 + PC4*Scaled_eGFRcys_2 + PC5*Scaled_eGFRcys_2 + PC6*Scaled_eGFRcys_2 +
                                    PC7*Scaled_eGFRcys_2 + PC8*Scaled_eGFRcys_2 + PC9*Scaled_eGFRcys_2 +
                                    PC10*Scaled_eGFRcys_2 + PC11*Scaled_eGFRcys_2 +
                                    PC12*Scaled_eGFRcys_2 + PC13*Scaled_eGFRcys_2 + PC14*Scaled_eGFRcys_2 + PC15*Scaled_eGFRcys_2 + PC16*Scaled_eGFRcys_2 + PC17*Scaled_eGFRcys_2 + PC18*Scaled_eGFRcys_2 + PC19*Scaled_eGFRcys_2 + PC20*Scaled_eGFRcys_2 +Scaled_eGFRcys_2*scale(Transport_SBP_PES_AVG),
                                  data = Merged_PGS_PES)

table1( ~ factor(Sex) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL + Cynstatin_C | factor(Binary_CKD),
        data = Merged_PGS_PES)

table1( ~ factor(Sex) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL + Cynstatin_C | factor(Binary_Stevens),
        data = Merged_PGS_PES)

## Test with continous adjustment for eGFR

Na_eGFR1_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*scale(Transport_SBP_PES_AVG) + Age*scale(Transport_SBP_PES_AVG) + Age2*scale(Transport_SBP_PES_AVG) + Assesment_centre_month*scale(Transport_SBP_PES_AVG) + 
                                    Assesment_centre*scale(Transport_SBP_PES_AVG) + PC1*scale(Transport_SBP_PES_AVG) + PC2*scale(Transport_SBP_PES_AVG) + PC3*scale(Transport_SBP_PES_AVG) + PC4*scale(Transport_SBP_PES_AVG) + PC5*scale(Transport_SBP_PES_AVG) + PC6*scale(Transport_SBP_PES_AVG) +
                                    PC7*scale(Transport_SBP_PES_AVG) + PC8*scale(Transport_SBP_PES_AVG) + PC9*scale(Transport_SBP_PES_AVG) +
                                    PC10*scale(Transport_SBP_PES_AVG) + PC11*scale(Transport_SBP_PES_AVG) +
                                    PC12*scale(Transport_SBP_PES_AVG) + PC13*scale(Transport_SBP_PES_AVG) + PC14*scale(Transport_SBP_PES_AVG) + PC15*scale(Transport_SBP_PES_AVG) + PC16*scale(Transport_SBP_PES_AVG) + PC17*scale(Transport_SBP_PES_AVG) + PC18*scale(Transport_SBP_PES_AVG) + PC19*scale(Transport_SBP_PES_AVG) + PC20*scale(Transport_SBP_PES_AVG) + 
                                    Sex*scale(Na_urine_millimolL) + Age*scale(Na_urine_millimolL) + Age2*scale(Na_urine_millimolL) + Assesment_centre_month*scale(Na_urine_millimolL) + 
                                    Assesment_centre*scale(Na_urine_millimolL) + PC1*scale(Na_urine_millimolL) + PC2*scale(Na_urine_millimolL) + PC3*scale(Na_urine_millimolL) + PC4*scale(Na_urine_millimolL) + PC5*scale(Na_urine_millimolL) + PC6*scale(Na_urine_millimolL) +
                                    PC7*scale(Na_urine_millimolL) + PC8*scale(Na_urine_millimolL) + PC9*scale(Na_urine_millimolL) +
                                    PC10*scale(Na_urine_millimolL) + PC11*scale(Na_urine_millimolL) +
                                    PC12*scale(Na_urine_millimolL) + PC13*scale(Na_urine_millimolL) + PC14*scale(Na_urine_millimolL) + PC15*scale(Na_urine_millimolL) + PC16*scale(Na_urine_millimolL) + PC17*scale(Na_urine_millimolL) + PC18*scale(Na_urine_millimolL) + PC19*scale(Na_urine_millimolL) + PC20*scale(Na_urine_millimolL) +scale(Na_urine_millimolL)*scale(Transport_SBP_PES_AVG) +
                                      Scaled_eGFRcys_1*scale(Na_urine_millimolL) + Scaled_eGFRcys_1*scale(Transport_SBP_PES_AVG),
                                  data = Merged_PGS_PES)

Na_eGFR2_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*scale(Transport_SBP_PES_AVG) + Age*scale(Transport_SBP_PES_AVG) + Age2*scale(Transport_SBP_PES_AVG) + Assesment_centre_month*scale(Transport_SBP_PES_AVG) + 
                                       Assesment_centre*scale(Transport_SBP_PES_AVG) + PC1*scale(Transport_SBP_PES_AVG) + PC2*scale(Transport_SBP_PES_AVG) + PC3*scale(Transport_SBP_PES_AVG) + PC4*scale(Transport_SBP_PES_AVG) + PC5*scale(Transport_SBP_PES_AVG) + PC6*scale(Transport_SBP_PES_AVG) +
                                       PC7*scale(Transport_SBP_PES_AVG) + PC8*scale(Transport_SBP_PES_AVG) + PC9*scale(Transport_SBP_PES_AVG) +
                                       PC10*scale(Transport_SBP_PES_AVG) + PC11*scale(Transport_SBP_PES_AVG) +
                                       PC12*scale(Transport_SBP_PES_AVG) + PC13*scale(Transport_SBP_PES_AVG) + PC14*scale(Transport_SBP_PES_AVG) + PC15*scale(Transport_SBP_PES_AVG) + PC16*scale(Transport_SBP_PES_AVG) + PC17*scale(Transport_SBP_PES_AVG) + PC18*scale(Transport_SBP_PES_AVG) + PC19*scale(Transport_SBP_PES_AVG) + PC20*scale(Transport_SBP_PES_AVG) + 
                                       Sex*scale(Na_urine_millimolL) + Age*scale(Na_urine_millimolL) + Age2*scale(Na_urine_millimolL) + Assesment_centre_month*scale(Na_urine_millimolL) + 
                                       Assesment_centre*scale(Na_urine_millimolL) + PC1*scale(Na_urine_millimolL) + PC2*scale(Na_urine_millimolL) + PC3*scale(Na_urine_millimolL) + PC4*scale(Na_urine_millimolL) + PC5*scale(Na_urine_millimolL) + PC6*scale(Na_urine_millimolL) +
                                       PC7*scale(Na_urine_millimolL) + PC8*scale(Na_urine_millimolL) + PC9*scale(Na_urine_millimolL) +
                                       PC10*scale(Na_urine_millimolL) + PC11*scale(Na_urine_millimolL) +
                                       PC12*scale(Na_urine_millimolL) + PC13*scale(Na_urine_millimolL) + PC14*scale(Na_urine_millimolL) + PC15*scale(Na_urine_millimolL) + PC16*scale(Na_urine_millimolL) + PC17*scale(Na_urine_millimolL) + PC18*scale(Na_urine_millimolL) + PC19*scale(Na_urine_millimolL) + PC20*scale(Na_urine_millimolL) +scale(Na_urine_millimolL)*scale(Transport_SBP_PES_AVG) +
                                       Scaled_eGFRcys_2*scale(Na_urine_millimolL) + Scaled_eGFRcys_2*scale(Transport_SBP_PES_AVG),
                                     data = Merged_PGS_PES)

## Split up into CKDEpi eGFR > 90 and eGFR < 90
## Over 90
Over_90_CKD_Merged_PGS_PES <- Merged_PGS_PES %>% filter(Binary_CKD == "> 90")
Over_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG))
Over_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))

## Test effect of scaled UNa
Over_90_CKD_Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$Na_urine_millimolL))
## UNa

UNa_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                  Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                                  PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                                  PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                                  PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                                  Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                  Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                                  PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                                  PC10*Scaled_Na + PC11*Scaled_Na +
                                  PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG,
                                data = Over_90_CKD_Merged_PGS_PES)

## Under 90

Under_90_CKD_Merged_PGS_PES <- Merged_PGS_PES %>% filter(Binary_CKD == "< 90")
Under_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG))
Under_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))
Under_90_CKD_Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$Na_urine_millimolL))


Under_90_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                  Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                                  PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                                  PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                                  PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                                  Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                  Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                                  PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                                  PC10*Scaled_Na + PC11*Scaled_Na +
                                  PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG,
                                data = Under_90_CKD_Merged_PGS_PES)

## Test if GxE estimate is significantly different between the two
## CKD eGFR > 90 = 5.745e-02  5.012e-02
## CKD eGFR < 90 = 9.479e-02  5.683e-02  

(5.745e-02-9.479e-02)/(sqrt((5.012e-02^2)+(5.683e-02^2)))
pnorm(abs(-0.4927827), lower.tail = F)*2

## No sig difference - 0.622

## Repeat above for the Stevens et al. eGFR
rm(Under_90_CKD_Merged_PGS_PES)
rm(Over_90_CKD_Merged_PGS_PES)

Over_90_CKD_Merged_PGS_PES <- Merged_PGS_PES %>% filter(Binary_Stevens == "> 90")
Over_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG))
Over_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))

## Test effect of scaled UNa
Over_90_CKD_Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Over_90_CKD_Merged_PGS_PES$Na_urine_millimolL))
## UNa

UNa_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                  Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                                  PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                                  PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                                  PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                                  Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                  Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                                  PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                                  PC10*Scaled_Na + PC11*Scaled_Na +
                                  PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG,
                                data = Over_90_CKD_Merged_PGS_PES)

## Under 90

Under_90_CKD_Merged_PGS_PES <- Merged_PGS_PES %>% filter(Binary_Stevens == "< 90")
Under_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG))
Under_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))
Under_90_CKD_Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Under_90_CKD_Merged_PGS_PES$Na_urine_millimolL))


Under_90_Transport_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                       Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG + PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC6*Transport_SBP_PES_AVG +
                                       PC7*Transport_SBP_PES_AVG + PC8*Transport_SBP_PES_AVG + PC9*Transport_SBP_PES_AVG +
                                       PC10*Transport_SBP_PES_AVG + PC11*Transport_SBP_PES_AVG +
                                       PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                                       Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                       Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na + PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC6*Scaled_Na +
                                       PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na +
                                       PC10*Scaled_Na + PC11*Scaled_Na +
                                       PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na +Scaled_Na*Transport_SBP_PES_AVG,
                                     data = Under_90_CKD_Merged_PGS_PES)

## Test if GxE estimate is significantly different between the two
## CKD eGFR > 90 = 7.058e-02  6.052e-02 
## CKD eGFR < 90 =  6.853e-02  4.799e-02

(7.058e-02-6.853e-02)/(sqrt((6.052e-02^2)+(4.799e-02^2)))
pnorm(abs(0.02654134), lower.tail = F)*2

## No significant difference between 0.9788256

## Make plots below showing the difference

Under_90_CKD_Merged_PGS_PES$Decile_PES <- as.numeric(ntile(Under_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG, 10))
Over_90_CKD_Merged_PGS_PES$Decile_PES <- as.numeric(ntile(Over_90_CKD_Merged_PGS_PES$Transport_SBP_PES_AVG, 10))


Na_assoc_per_decile <- function(decile, urinary, score, df) {
  new_df <- df[df[[score]] == decile,]
  new_df$scaled_urinary <- as.numeric(scale(new_df[[urinary]]))
  mod <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
              Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
              PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + scaled_urinary, data = new_df)
  return(summary(mod))
}

Vec_decile <- c(1:10)

## PES under

PES_under <- mapply(Na_assoc_per_decile, decile =  Vec_decile, 
              MoreArgs = list("Na_urine_millimolL","Decile_PES", Under_90_CKD_Merged_PGS_PES), SIMPLIFY = T)

## PES over

PES_over <- mapply(Na_assoc_per_decile, decile =  Vec_decile, 
                    MoreArgs = list("Na_urine_millimolL","Decile_PES", Over_90_CKD_Merged_PGS_PES), SIMPLIFY = T)

PES_under_extract <- apply(PES_under, 2, function(x) return(as.data.frame(x$coefficients)[56, 1:4]))
PES_over_extract <- apply(PES_over, 2, function(x) return(as.data.frame(x$coefficients)[56, 1:4]))

## Df for each

PES_over_plot <- data.frame()

for (i in 1:length(PES_over_extract)) {
  PES_over_plot <- rbind(PES_over_plot, PES_over_extract[[i]])
}

Decile <- rep(c("1st-10th percentile", "10th-20th percentile", "20th-30th percentile", "30th-40th percentile", "40th-50th percentile",
                "50th-60th percentile", "60th-70th percentile", "70th-80th percentile", "80th-90th percentile", "Top decile"),each=1)

PES_over_plot$Decile <- Decile

PES_under_plot <- data.frame()

for (i in 1:length(PES_under_extract)) {
  PES_under_plot <- rbind(PES_under_plot, PES_under_extract[[i]])
}

Decile <- rep(c("1st-10th percentile", "10th-20th percentile", "20th-30th percentile", "30th-40th percentile", "40th-50th percentile",
                "50th-60th percentile", "60th-70th percentile", "70th-80th percentile", "80th-90th percentile", "Top decile"),each=1)

PES_under_plot$Decile <- Decile

## Statistical comparison of deciles between the two different eGFR conditions

Stat_compare_eGFR <- function(df_1, df_2, decile) {
  new_df_1 <- df_1
  new_df_2 <- df_2
  Z <- (new_df_1$Estimate[[decile]]-new_df_2$Estimate[[decile]])/(sqrt((new_df_1$`Std. Error`[[decile]]^2)+(new_df_2$`Std. Error`[[decile]]^2)))
  P <- pnorm(abs(Z), lower.tail = F)*2
  Z_out <- as.data.frame(Z)
  P_out <- as.data.frame(P)
  Out <- as.data.frame(cbind(Z_out, P_out))
  return(Out)
}

Val <- c(1,2,3,4,5,6,7,8,9,10)
sapply(Val, Stat_compare_eGFR, df_1 = PES_over_plot, df_2 = PES_under_plot)

Stat_compare <- function(df, decile) {
  new_df <- df
  Z <- (new_df$Estimate[[10]]-new_df$Estimate[[decile]])/(sqrt((new_df$`Std. Error`[[10]]^2)+(new_df$`Std. Error`[[decile]]^2)))
  P <- pnorm(abs(Z), lower.tail = F)*2
  Z_out <- as.data.frame(Z)
  P_out <- as.data.frame(P)
  Out <- as.data.frame(cbind(Z_out, P_out))
  return(Out)
}

sapply(Val, Stat_compare, df=PES_over_plot)
sapply(Val, Stat_compare, df=PES_under_plot)

## Make plots
PES_over_plot$L_CI <- PES_over_plot$Estimate - (1.96*PES_over_plot$`Std. Error`)
PES_over_plot$U_CI <- PES_over_plot$Estimate + (1.96*PES_over_plot$`Std. Error`)
PES_over_plot$Val <- c(1,2,3,4,5,6,7,8,9,10)
mean(PES_over_plot$Estimate)
mean(PES_over_plot$Estimate)-sd(PES_over_plot$Estimate)
mean(PES_over_plot$Estimate)+sd(PES_over_plot$Estimate)  

P1 <- ggplot(PES_over_plot, aes(Val, Estimate)) +
  geom_errorbar(aes(ymin = L_CI, ymax = U_CI, color = as.factor(Val))) +
  geom_point(aes(color = as.factor(Val)), position = position_dodge(0.3)) +
  #stat_smooth(method = "loess", se = F, size = 0.5, colour = "black") +
  stat_smooth(method = "lm", se = F, size = 0.5, colour = "blue") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD") +
  xlab("Transport PES (deciles)") +
  geom_hline(yintercept = 1.02, lty = "dashed") +
  ggtitle("UNa+ (eGFRcys > 90 ml/min/1.73m2)")

PES_under_plot$L_CI <- PES_under_plot$Estimate - (1.96*PES_under_plot$`Std. Error`)
PES_under_plot$U_CI <- PES_under_plot$Estimate + (1.96*PES_under_plot$`Std. Error`)
PES_under_plot$Val <- c(1,2,3,4,5,6,7,8,9,10)
mean(PES_under_plot$Estimate)
mean(PES_under_plot$Estimate)-sd(PES_under_plot$Estimate)
mean(PES_under_plot$Estimate)+sd(PES_under_plot$Estimate)  

P2 <- ggplot(PES_under_plot, aes(Val, Estimate)) +
  geom_errorbar(aes(ymin = L_CI, ymax = U_CI, color = as.factor(Val))) +
  geom_point(aes(color = as.factor(Val)), position = position_dodge(0.3)) +
  #stat_smooth(method = "loess", se = F, size = 0.5, colour = "black") +
  stat_smooth(method = "lm", se = F, size = 0.5, colour = "blue") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD") +
  xlab("Transport PES (deciles)") +
  geom_hline(yintercept = 1.244, lty = "dashed") +
  ggtitle("UNa+ (eGFRcys < 90 ml/min/1.73m2)")

Na_final_1 <- P1 + scale_x_continuous(breaks = 1:10,
                                    labels = paste(c("1st", "2nd", "3rd", "4th", 
                                                     "5th", "6th", "7th", "8th", "9th", "10th"), " Decile")) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_aaas()

Na_final_2 <- P2 + scale_x_continuous(breaks = 1:10,
                                      labels = paste(c("1st", "2nd", "3rd", "4th", 
                                                       "5th", "6th", "7th", "8th", "9th", "10th"), " Decile")) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_aaas()

ggarrange(Na_final_1, Na_final_2)


