#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 1: Step 3 - cross-sectional relationship between Na, K, and Na:K with BP

## William Reay (2022)

#################################################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggridges)

## Read in data generated from previous script

Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")

table1( ~ factor(Sex) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL | factor(`BP Medication`),
        data = Subcohort_1_dat)

Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)

## Test assoc between raw and ln transformed Na and K with SBP and DBP separately.

## Adjust for age, age2, sex, assesesment centre, month of visit, blood pressure meds

## SBP ##

Raw_Na_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + Na_urine_millimolL, data = Subcohort_1_dat)

Raw_K_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + K_urine_millimolL, data = Subcohort_1_dat)

Ln_Na_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + Log_Na, data = Subcohort_1_dat)

Ln_K_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + Log_K, data = Subcohort_1_dat)

## DBP ##

Raw_Na_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + Na_urine_millimolL, data = Subcohort_1_dat)

Raw_K_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + K_urine_millimolL, data = Subcohort_1_dat)

Ln_Na_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + Log_Na, data = Subcohort_1_dat)

Ln_K_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                     Assesment_centre + Log_K, data = Subcohort_1_dat)

## Plot regression line of effect - first with no covariates

Na_plot <- ggplot(data = Subcohort_1_dat, aes(x = Na_urine_millimolL, y = Mean_SBP, colour = BP_medication)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Urinary Sodium (millimol/L)") +
  theme(legend.position = "none") +
  labs(colour = "BP Medication at Baseline")

K_plot <- ggplot(data = Subcohort_1_dat, aes(x = K_urine_millimolL, y = Mean_SBP, colour = BP_medication)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab(" ") +
  theme_bw() +
  xlab("Spot Urinary Potassium (millimol/L)") +
  theme(legend.position = "none") +
  labs(colour = "BP Medication at Baseline")

ggarrange(Na_plot, K_plot, common.legend = T)

## Per month of assessment centre visit

Na_plot_2 <- ggplot(data = Subcohort_1_dat, aes(x = Na_urine_millimolL, y = Mean_SBP, colour = Assesment_centre_month)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Urinary Sodium (millimol/L)") +
  theme(legend.position = "none") +
  labs(colour = "Month visited assessment centre")

K_plot_2 <- ggplot(data = Subcohort_1_dat, aes(x = K_urine_millimolL, y = Mean_SBP, colour = Assesment_centre_month)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab(" ") +
  theme_bw() +
  xlab("Spot Urinary Potassium (millimol/L)") +
  theme(legend.position = "none") +
  labs(colour = "Month visited assessment centre")

ggarrange(Na_plot_2, K_plot_2, common.legend = T)

## Sex strat

Na_plot_sex <- ggplot(data = Subcohort_1_dat, aes(x = Na_urine_millimolL, y = Mean_SBP, colour = Sex)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Urinary Sodium (millimol/L)") +
  theme(legend.position = "none")

K_plot_sex <- ggplot(data = Subcohort_1_dat, aes(x = K_urine_millimolL, y = Mean_SBP, colour = Sex)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab(" ") +
  theme_bw() +
  xlab("Spot Urinary Potassium (millimol/L)") +
  theme(legend.position = "none") 

ggarrange(Na_plot_sex, K_plot_sex, common.legend = T)


## Test Na:K ratio > 1 vs Na:K ratio < 1

Subcohort_1_dat$Na_K_ratio_binary <- ifelse(Subcohort_1_dat$Na_K_ratio > 1, "Over_1", "Below_1")

Ratio_SBP_binary <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + Na_K_ratio, data = Subcohort_1_dat)

Ratio_DBP_binary <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                  Assesment_centre + Na_K_ratio, data = Subcohort_1_dat)



## Plot

Ratio_1 <- ggplot(data = Subcohort_1_dat, aes(x=Na_K_ratio, y = Mean_SBP, fill = Na_K_ratio)) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "#56B4E9") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("> 1","<= 1")) +
  ylab("Mean SBP (mmHg)") +
  xlab("Spot Urinary Na/K ratio")

Ratio_2 <- ggplot(data = Subcohort_1_dat, aes(x=Na_K_ratio, y = Mean_DBP, fill = Na_K_ratio)) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 82.25, linetype = "dashed", colour = "darkred") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("> 1","<= 1")) +
  ylab("Mean DBP (mmHg)") +
  xlab("Spot Urinary Na/K ratio")

ggarrange(Ratio_1, Ratio_2)


## Adjust for BMI

BMI_adj_Na_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                    Assesment_centre + BMI + Na_urine_millimolL, data = Subcohort_1_dat)

BMI_adj_Na_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                    Assesment_centre + BMI + Na_urine_millimolL, data = Subcohort_1_dat)

BMI_adj_K_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + K_urine_millimolL, data = Subcohort_1_dat)

BMI_adj_K_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + K_urine_millimolL, data = Subcohort_1_dat)

## Adjust for BMI and creatinine

Creat_BMI_adj_Na_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + Creatinine_urine_micromolL + Na_urine_millimolL, data = Subcohort_1_dat)

Creat_BMI_adj_Na_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + Creatinine_urine_micromolL + Na_urine_millimolL, data = Subcohort_1_dat)

Creat_BMI_adj_K_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + BMI + Creatinine_urine_micromolL + K_urine_millimolL, data = Subcohort_1_dat)

Creat_BMI_adj_K_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + BMI + Creatinine_urine_micromolL + K_urine_millimolL, data = Subcohort_1_dat)

## Better ratio plot

breaks <- c(0, 1, 2, 3, 4, 5, 16)

tags <- c("(< 1)", "(1-2)", "(2-3)", "(3-4)", "(3-4)", "(> 5)")

Subcohort_1_dat$Cut_ratio <- cut(Subcohort_1_dat$Na_K_ratio, 
                                 breaks=breaks, 
                                 include.lowest=TRUE, 
                                 right=T, 
                                 labels=tags)

Better_ratio_1 <- ggplot(data = Subcohort_1_dat, aes(x=as.factor(Cut_ratio), y = Mean_SBP, fill = as.factor(Cut_ratio))) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "darkred") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Mean SBP (mmHg)") +
  xlab("Spot Urinary Na/K ratio")

Better_ratio_2 <- ggplot(data = Subcohort_1_dat, aes(x=as.factor(Cut_ratio), y = Mean_DBP, fill = as.factor(Cut_ratio))) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 82.25, linetype = "dashed", colour = "darkred") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Mean DBP (mmHg)") +
  xlab("Spot Urinary Na/K ratio")

ggarrange(Better_ratio_1, Better_ratio_2)

## Stratify by meds

## Remove medicated individuals

Unmed_cohort <- Subcohort_1_dat %>% filter(BP_meds == 0)

Raw_Na_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                       Assesment_centre + Na_urine_millimolL, data = Unmed_cohort)

Raw_K_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                      Assesment_centre + K_urine_millimolL, data = Unmed_cohort)

## DBP ##

Raw_Na_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                       Assesment_centre + Na_urine_millimolL, data = Unmed_cohort)

Raw_K_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                      Assesment_centre + K_urine_millimolL, data = Unmed_cohort)

## Ridge plot strat by meds

R1<- ggplot(Subcohort_1_dat, aes(x = Mean_SBP, y = as.factor(Cut_ratio), fill = BP_medication)) +
  geom_density_ridges(alpha = 0.3) +
  theme_ridges() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 128.228, linetype = "dashed") + 
  xlab("Mean SBP (mmHg)") +
  ylab("Spot Urinary Na/K ratio") +
  theme(legend.position = "none") +
  labs(fill = "BP Medication at Baseline")
  
R2<- ggplot(Subcohort_1_dat, aes(x = Mean_DBP, y = as.factor(Cut_ratio), fill = BP_medication)) +
  geom_density_ridges(alpha = 0.3) +
  theme_ridges() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 82.25, linetype = "dashed") + 
  xlab("Mean DBP (mmHg)") +
  ylab(" ") +
  theme(legend.position = "none") +
  labs(fill = "BP Medication at Baseline")

ggarrange(R1, R2, common.legend = T, legend = "bottom")
