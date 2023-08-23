#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2c - exploring transport PES x sodium interaction at diff thresholds - UNMEDICATED ONLY

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

Na_K_PES <- fread("UKBB_PES/PES_BP_assoc_and_int/Transport_PES_sensitivity_scores/Sensitivity_diff_thresholds_transport_PES.txt")

Colnames_PES <- make.unique(names(Na_K_PES))

colnames(Na_K_PES) <- Colnames_PES

Na_K_PES <- Na_K_PES %>% select(IID, contains("PES"))

Merged_PES <- merge(Na_K_PES, Subcohort_1_dat, by = "IID")

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")
Merged_PGS_PES$Transport_SBP_PES_1_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_1_AVG))
Merged_PGS_PES$Transport_SBP_PES_0_5_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_0_5_AVG))
Merged_PGS_PES$Transport_SBP_PES_0_05_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_0_05_AVG))

Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_urine_millimolL))


## Pt < 1

Ratio_Transport_int_G_C_C_1 <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_1_AVG + Age*Transport_SBP_PES_1_AVG + Age2*Transport_SBP_PES_1_AVG + Assesment_centre_month*Transport_SBP_PES_1_AVG + 
                                    Assesment_centre*Transport_SBP_PES_1_AVG + PC1*Transport_SBP_PES_1_AVG + PC2*Transport_SBP_PES_1_AVG+ 
                                    PC3*Transport_SBP_PES_1_AVG + PC4*Transport_SBP_PES_1_AVG + PC5*Transport_SBP_PES_1_AVG + 
                                    PC6*Transport_SBP_PES_1_AVG + PC7*Transport_SBP_PES_1_AVG + PC8*Transport_SBP_PES_1_AVG + PC9*Transport_SBP_PES_1_AVG + 
                                    PC10*Transport_SBP_PES_1_AVG+ PC11*Transport_SBP_PES_1_AVG +
                                    PC12*Transport_SBP_PES_1_AVG + PC13*Transport_SBP_PES_1_AVG + PC14*Transport_SBP_PES_1_AVG + PC15*Transport_SBP_PES_1_AVG + 
                                    PC16*Transport_SBP_PES_1_AVG + PC17*Transport_SBP_PES_1_AVG + PC18*Transport_SBP_PES_1_AVG + PC19*Transport_SBP_PES_1_AVG + PC20*Transport_SBP_PES_1_AVG + 
                                    Scaled_Na*Transport_SBP_PES_1_AVG +
                                    Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                    Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + 
                                    PC6*Scaled_Na + PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                                    PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                                  data = Merged_PGS_PES)

summary(Ratio_Transport_int_G_C_C_1)
## P < 1
#Transport_SBP_PES_1_AVG:Scaled_Na                -1.114e-02  3.652e-02  -0.305 0.760352  

## P < 0.5

Ratio_Transport_int_G_C_C_0_5 <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_0_5_AVG + Age*Transport_SBP_PES_0_5_AVG + Age2*Transport_SBP_PES_0_5_AVG + Assesment_centre_month*Transport_SBP_PES_0_5_AVG + 
                                    Assesment_centre*Transport_SBP_PES_0_5_AVG + PC1*Transport_SBP_PES_0_5_AVG + PC2*Transport_SBP_PES_0_5_AVG+ 
                                    PC3*Transport_SBP_PES_0_5_AVG + PC4*Transport_SBP_PES_0_5_AVG + PC5*Transport_SBP_PES_0_5_AVG + 
                                    PC6*Transport_SBP_PES_0_5_AVG + PC7*Transport_SBP_PES_0_5_AVG + PC8*Transport_SBP_PES_0_5_AVG + PC9*Transport_SBP_PES_0_5_AVG + 
                                    PC10*Transport_SBP_PES_0_5_AVG+ PC11*Transport_SBP_PES_0_5_AVG +
                                    PC12*Transport_SBP_PES_0_5_AVG + PC13*Transport_SBP_PES_0_5_AVG + PC14*Transport_SBP_PES_0_5_AVG + PC15*Transport_SBP_PES_0_5_AVG + 
                                    PC16*Transport_SBP_PES_0_5_AVG + PC17*Transport_SBP_PES_0_5_AVG + PC18*Transport_SBP_PES_0_5_AVG + PC19*Transport_SBP_PES_0_5_AVG + PC20*Transport_SBP_PES_0_5_AVG + 
                                    Scaled_Na*Transport_SBP_PES_0_5_AVG +
                                    Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                    Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + 
                                    PC6*Scaled_Na + PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                                    PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                                  data = Merged_PGS_PES)

summary(Ratio_Transport_int_G_C_C_0_5)

## Transport_SBP_PES_0_5_AVG:Scaled_Na                -8.113e-03  3.648e-02  -0.222 0.824020  

## P < 0.05

Ratio_Transport_int_G_C_C_0_05 <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_0_05_AVG + Age*Transport_SBP_PES_0_05_AVG + Age2*Transport_SBP_PES_0_05_AVG + Assesment_centre_month*Transport_SBP_PES_0_05_AVG + 
                                      Assesment_centre*Transport_SBP_PES_0_05_AVG + PC1*Transport_SBP_PES_0_05_AVG + PC2*Transport_SBP_PES_0_05_AVG+ 
                                      PC3*Transport_SBP_PES_0_05_AVG + PC4*Transport_SBP_PES_0_05_AVG + PC5*Transport_SBP_PES_0_05_AVG + 
                                      PC6*Transport_SBP_PES_0_05_AVG + PC7*Transport_SBP_PES_0_05_AVG + PC8*Transport_SBP_PES_0_05_AVG + PC9*Transport_SBP_PES_0_05_AVG + 
                                      PC10*Transport_SBP_PES_0_05_AVG+ PC11*Transport_SBP_PES_0_05_AVG +
                                      PC12*Transport_SBP_PES_0_05_AVG + PC13*Transport_SBP_PES_0_05_AVG + PC14*Transport_SBP_PES_0_05_AVG + PC15*Transport_SBP_PES_0_05_AVG + 
                                      PC16*Transport_SBP_PES_0_05_AVG + PC17*Transport_SBP_PES_0_05_AVG + PC18*Transport_SBP_PES_0_05_AVG + PC19*Transport_SBP_PES_0_05_AVG + PC20*Transport_SBP_PES_0_05_AVG + 
                                      Scaled_Na*Transport_SBP_PES_0_05_AVG +
                                      Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                      Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + 
                                      PC6*Scaled_Na + PC7*Scaled_Na + PC8*Scaled_Na + PC9*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                                      PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                                    data = Merged_PGS_PES)

summary(Ratio_Transport_int_G_C_C_0_05)

## Transport_SBP_PES_0_05_AVG:Scaled_Na                 8.955e-04  3.656e-02   0.024  0.98046 


