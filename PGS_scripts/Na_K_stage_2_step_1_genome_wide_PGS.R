#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 1 - PGS impact on UKBB, any interplay with Na, K or the ratio

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

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS <- merge(PGS_BP, Subcohort_1_dat, by = "IID")

## Test effect of SBP PGS on SBP and DBP PGS on DBP
## Compare effect size of urinary measures per SD to PGS per SD

Merged_PGS$Scaled_SBP_PGS_full <- as.numeric(scale(Merged_PGS$SBP_PGS_Pt_0_01_AVG))
Merged_PGS$Scaled_DBP_PGS_full <- as.numeric(scale(Merged_PGS$DBP_PGS_Pt_0_005_AVG))


## MEDICATED AND UNMEDICATED

SBP_PGS_baseline <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                         Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Scaled_SBP_PGS_full, data = Merged_PGS)

SBP_PGS_baseline_with_sodium <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +  PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                     Scaled_SBP_PGS_full + Na_urine_millimolL, data = Merged_PGS)

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


DBP_PGS_baseline <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                         PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_DBP_PGS_full, data = Merged_PGS)

SBP_DBP_PGS_corr <- cor(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Scaled_DBP_PGS_full)

## Test whether adding SBP and DBP PGS improves models

SBP_PGS_baseline_with_DBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_SBP_PGS_full + Scaled_DBP_PGS_full, data = Merged_PGS)

anova(SBP_PGS_baseline, SBP_PGS_baseline_with_DBP, test = "F")

## Are the PGS correlated with any of the urinary measures

cor.test(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Na_urine_millimolL)
cor.test(Merged_PGS$Scaled_DBP_PGS_full, Merged_PGS$Na_urine_millimolL)
cor.test(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$K_urine_millimolL)
cor.test(Merged_PGS$Scaled_DBP_PGS_full, Merged_PGS$K_urine_millimolL)
