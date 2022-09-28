#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 1 - PGS impact on UKBB, any interplay with Na, K or the ratio

## UNMEDICATED COHORT

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

Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")

Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS <- merge(PGS_BP, Subcohort_1_dat, by = "IID")

Merged_PGS <- Merged_PGS %>% filter(BP_meds == 0)

Merged_PGS$Scaled_SBP_PGS_full <- as.numeric(scale(Merged_PGS$SBP_PGS_Pt_0_01_AVG))
Merged_PGS$Scaled_DBP_PGS_full <- as.numeric(scale(Merged_PGS$DBP_PGS_Pt_0_005_AVG))

## REPEAT ANALYSES IN UNMEDICATED

SBP_PGS_baseline <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                         Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Scaled_SBP_PGS_full, data = Merged_PGS)

SBP_PGS_baseline_with_sodium <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                     Scaled_SBP_PGS_full + Na_urine_millimolL, data = Merged_PGS)

anova(SBP_PGS_baseline, SBP_PGS_baseline_with_sodium, test="F")


DBP_PGS_baseline <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                         PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_DBP_PGS_full, data = Merged_PGS)

## INTERACTIONS - scale Na and K to 1 so more interpretable

Merged_PGS$Scaled_K_full <- as.numeric(scale(Merged_PGS$K_urine_millimolL))
Merged_PGS$Scaled_Na_full <- as.numeric(scale(Merged_PGS$Na_urine_millimolL))

## SBP ~ Na
SBP_PGS_baseline_Na_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_SBP_PGS_full*Scaled_Na_full, data = Merged_PGS)
## DBP ~ Na

DBP_PGS_baseline_Na_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_DBP_PGS_full*Scaled_Na_full, data = Merged_PGS)

## SBP ~ K

SBP_PGS_baseline_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_SBP_PGS_full*Scaled_K_full, data = Merged_PGS)



SBP_PGS_baseline_K_int_all_cov_int <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                           PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full + PC3*Scaled_SBP_PGS_full + 
                                           PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC10*Scaled_SBP_PGS_full + PC11*Scaled_SBP_PGS_full +
                                           PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full +
                                           Assesment_centre*Scaled_SBP_PGS_full + Scaled_SBP_PGS_full*Scaled_K_full, data = Merged_PGS)

## DBP ~ K

DBP_PGS_baseline_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_DBP_PGS_full*Scaled_K_full, data = Merged_PGS)



## SBP ~ Na:K

SBP_PGS_baseline_Na_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                             PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                             PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                             Assesment_centre + Scaled_SBP_PGS_full*Na_K_ratio, data = Merged_PGS)

DBP_PGS_baseline_Na_K_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                             PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                             PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                             Assesment_centre + Scaled_DBP_PGS_full*Na_K_ratio, data = Merged_PGS)

## Quartile K/SBP plot

## Define high PGS

quantile(Merged_PGS$Scaled_SBP_PGS_full, probs = seq(.1, .9, by = .1))

Merged_PGS$High_PGS <- ifelse(Merged_PGS$Scaled_SBP_PGS_full > 1.279041346, "High PGS", "Comparator PGS (< 90th percentile)")

Merged_PGS$K_quartile <- ntile(Merged_PGS$K_urine_millimolL, 4)

ggplot(data = Merged_PGS, aes(x=as.factor(K_quartile), y = Mean_SBP, fill = as.factor(High_PGS))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "firebrick") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "plum")) +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4")) +
  ylab("Mean SBP (mmHg)") +
  xlab("Urinary potassium (millimol/L)") +
  ggtitle("Genome-wide PGS - SBP (unmedicated cohort)") +
  labs(fill = " ")

cor.test(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Scaled_K_full)
