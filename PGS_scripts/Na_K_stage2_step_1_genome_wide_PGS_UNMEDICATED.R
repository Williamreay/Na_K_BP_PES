#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 1 - PGS impact on UKBB, any interplay with Na, K or the ratio

## UNMEDICATED

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

## Read in files, set relevant variables as categorical

Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")

Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)
Unmed <- Subcohort_1_dat %>% filter(BP_meds == 0)

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS <- merge(PGS_BP, Unmed, by = "IID")

## Test effect of SBP PGS on SBP and DBP PGS on DBP
## Compare effect size of urinary measures per SD to PGS per SD

Merged_PGS$Scaled_SBP_PGS_full <- as.numeric(scale(Merged_PGS$SBP_PGS_Pt_0_01_AVG))
Merged_PGS$Scaled_DBP_PGS_full <- as.numeric(scale(Merged_PGS$DBP_PGS_Pt_0_005_AVG))


## MEDICATED AND UNMEDICATED

SBP_PGS_baseline <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                         Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Scaled_SBP_PGS_full, data = Merged_PGS)

SBP_sodium <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +  PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                     Na_urine_millimolL, data = Merged_PGS)

SBP_PGS_baseline_with_sodium <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +  PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                    Na_urine_millimolL + Scaled_SBP_PGS_full, data = Merged_PGS)

anova(SBP_PGS_baseline_with_sodium, SBP_sodium, test="F")

Merged_PGS$Quintile <- ntile(Merged_PGS$Scaled_SBP_PGS_full, 5)

Figure_plot <- ggplot(data = Merged_PGS, aes(x = as.factor(Quintile), y = Mean_SBP,
                                             fill = as.factor(Quintile))) + 
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("SBP PGS (quintiles)") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 138.228, lty="dashed", colour = "red")

Figure_plot + scale_fill_brewer(palette = "Paired")

anova(SBP_PGS_baseline, SBP_PGS_baseline_with_sodium, test="F")


DBP_PGS_baseline <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                         PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_DBP_PGS_full, data = Merged_PGS)

SBP_DBP_PGS_corr <- cor(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Scaled_DBP_PGS_full)


## Is there a non-linear interaction between the PGS and the urinary measures: code the variables as continous as as high + low PES

## Scale Urinary measures to SD = 1 so more interpretable interaction

Merged_PGS$Scaled_K_full <- as.numeric(scale(Merged_PGS$K_urine_millimolL))
Merged_PGS$Scaled_Na_full <- as.numeric(scale(Merged_PGS$Na_urine_millimolL))

## SBP ~ Na
SBP_PGS_baseline_Na_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_SBP_PGS_full*Scaled_Na_full, data = Merged_PGS)


## Add GxC and ExC

SBP_PGS_baseline_Na_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                        Assesment_centre*Scaled_SBP_PGS_full + PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full + PC3*Scaled_SBP_PGS_full + PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC6*Scaled_SBP_PGS_full +
                                        PC7*Scaled_SBP_PGS_full + PC8*Scaled_SBP_PGS_full + PC9*Scaled_SBP_PGS_full +
                                        PC10*Scaled_SBP_PGS_full + PC11*Scaled_SBP_PGS_full +
                                        PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full + 
                                        Sex*Scaled_Na_full + Age*Scaled_Na_full + Age2*Scaled_Na_full + Assesment_centre_month*Scaled_Na_full + 
                                        Assesment_centre*Scaled_Na_full + PC1*Scaled_Na_full + PC2*Scaled_Na_full + PC3*Scaled_Na_full + PC4*Scaled_Na_full + PC5*Scaled_Na_full + PC6*Scaled_Na_full +
                                        PC7*Scaled_Na_full + PC8*Scaled_Na_full + PC9*Scaled_Na_full +
                                        PC10*Scaled_Na_full + PC11*Scaled_Na_full +
                                        PC12*Scaled_Na_full + PC13*Scaled_Na_full + PC14*Scaled_Na_full + PC15*Scaled_Na_full + PC16*Scaled_Na_full + PC17*Scaled_Na_full + PC18*Scaled_Na_full + PC19*Scaled_Na_full + PC20*Scaled_Na_full +Scaled_Na_full*Scaled_SBP_PGS_full,
                                      data = Merged_PGS)


## DBP ~ Na

DBP_PGS_baseline_Na_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_DBP_PGS_full*Scaled_Na_full, data = Merged_PGS)


DBP_PGS_baseline_Na_int_G_C_E_C <- lm(Mean_DBP ~ Sex*Scaled_DBP_PGS_full + Age*Scaled_DBP_PGS_full + Age2*Scaled_DBP_PGS_full + Assesment_centre_month*Scaled_DBP_PGS_full + 
                                        Assesment_centre*Scaled_DBP_PGS_full + PC1*Scaled_DBP_PGS_full + PC2*Scaled_DBP_PGS_full + PC3*Scaled_DBP_PGS_full + PC4*Scaled_DBP_PGS_full + PC5*Scaled_DBP_PGS_full + PC6*Scaled_DBP_PGS_full +
                                        PC7*Scaled_DBP_PGS_full + PC8*Scaled_DBP_PGS_full + PC9*Scaled_DBP_PGS_full +
                                        PC10*Scaled_DBP_PGS_full + PC11*Scaled_DBP_PGS_full +
                                        PC12*Scaled_DBP_PGS_full + PC13*Scaled_DBP_PGS_full + PC14*Scaled_DBP_PGS_full + PC15*Scaled_DBP_PGS_full + PC16*Scaled_DBP_PGS_full + PC17*Scaled_DBP_PGS_full + PC18*Scaled_DBP_PGS_full + PC19*Scaled_DBP_PGS_full + PC20*Scaled_DBP_PGS_full + 
                                        Sex*Scaled_Na_full + Age*Scaled_Na_full + Age2*Scaled_Na_full + Assesment_centre_month*Scaled_Na_full + 
                                        Assesment_centre*Scaled_Na_full + PC1*Scaled_Na_full + PC2*Scaled_Na_full + PC3*Scaled_Na_full + PC4*Scaled_Na_full + PC5*Scaled_Na_full + PC6*Scaled_Na_full +
                                        PC7*Scaled_Na_full + PC8*Scaled_Na_full + PC9*Scaled_Na_full +
                                        PC10*Scaled_Na_full + PC11*Scaled_Na_full +
                                        PC12*Scaled_Na_full + PC13*Scaled_Na_full + PC14*Scaled_Na_full + PC15*Scaled_Na_full + PC16*Scaled_Na_full + PC17*Scaled_Na_full + PC18*Scaled_Na_full + PC19*Scaled_Na_full + PC20*Scaled_Na_full +Scaled_Na_full*Scaled_DBP_PGS_full,
                                      data = Merged_PGS)

## Repeat above for K+

SBP_PGS_baseline_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_SBP_PGS_full*Scaled_K_full, data = Merged_PGS)


## Add GxC and ExC

SBP_PGS_baseline_K_int_G_C_E_C <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                       Assesment_centre*Scaled_SBP_PGS_full + PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full + PC3*Scaled_SBP_PGS_full + PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC6*Scaled_SBP_PGS_full +
                                       PC7*Scaled_SBP_PGS_full + PC8*Scaled_SBP_PGS_full + PC9*Scaled_SBP_PGS_full +
                                       PC10*Scaled_SBP_PGS_full + PC11*Scaled_SBP_PGS_full +
                                       PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full + 
                                       Sex*Scaled_K_full + Age*Scaled_K_full + Age2*Scaled_K_full + Assesment_centre_month*Scaled_K_full + 
                                       Assesment_centre*Scaled_K_full + PC1*Scaled_K_full + PC2*Scaled_K_full + PC3*Scaled_K_full + PC4*Scaled_K_full + PC5*Scaled_K_full + PC6*Scaled_K_full +
                                       PC7*Scaled_K_full + PC8*Scaled_K_full + PC9*Scaled_K_full +
                                       PC10*Scaled_K_full + PC11*Scaled_K_full +
                                       PC12*Scaled_K_full + PC13*Scaled_K_full + PC14*Scaled_K_full + PC15*Scaled_K_full + PC16*Scaled_K_full + PC17*Scaled_K_full + PC18*Scaled_K_full + PC19*Scaled_K_full + PC20*Scaled_K_full +Scaled_K_full*Scaled_SBP_PGS_full,
                                     data = Merged_PGS)


## DBP ~ K

DBP_PGS_baseline_K_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_DBP_PGS_full*Scaled_K_full, data = Merged_PGS)


## Add GxC and ExC

DBP_PGS_baseline_K_int_G_C_E_C <- lm(Mean_DBP ~ Sex*Scaled_DBP_PGS_full + Age*Scaled_DBP_PGS_full + Age2*Scaled_DBP_PGS_full + Assesment_centre_month*Scaled_DBP_PGS_full + 
                                       Assesment_centre*Scaled_DBP_PGS_full + PC1*Scaled_DBP_PGS_full + PC2*Scaled_DBP_PGS_full + PC3*Scaled_DBP_PGS_full + PC4*Scaled_DBP_PGS_full + PC5*Scaled_DBP_PGS_full + PC6*Scaled_DBP_PGS_full +
                                       PC7*Scaled_DBP_PGS_full + PC8*Scaled_DBP_PGS_full + PC9*Scaled_DBP_PGS_full +
                                       PC10*Scaled_DBP_PGS_full + PC11*Scaled_DBP_PGS_full +
                                       PC12*Scaled_DBP_PGS_full + PC13*Scaled_DBP_PGS_full + PC14*Scaled_DBP_PGS_full + PC15*Scaled_DBP_PGS_full + PC16*Scaled_DBP_PGS_full + PC17*Scaled_DBP_PGS_full + PC18*Scaled_DBP_PGS_full + PC19*Scaled_DBP_PGS_full + PC20*Scaled_DBP_PGS_full + 
                                       Sex*Scaled_K_full + Age*Scaled_K_full + Age2*Scaled_K_full + Assesment_centre_month*Scaled_K_full + 
                                       Assesment_centre*Scaled_K_full + PC1*Scaled_K_full + PC2*Scaled_K_full + PC3*Scaled_K_full + PC4*Scaled_K_full + PC5*Scaled_K_full + PC6*Scaled_K_full +
                                       PC7*Scaled_K_full + PC8*Scaled_K_full + PC9*Scaled_K_full +
                                       PC10*Scaled_K_full + PC11*Scaled_K_full +
                                       PC12*Scaled_K_full + PC13*Scaled_K_full + PC14*Scaled_K_full + PC15*Scaled_K_full + PC16*Scaled_K_full + PC17*Scaled_K_full + PC18*Scaled_K_full + PC19*Scaled_K_full + PC20*Scaled_K_full +Scaled_K_full*Scaled_DBP_PGS_full,
                                     data = Merged_PGS)

