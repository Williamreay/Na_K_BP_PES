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
                         Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Scaled_SBP_PGS_full, data = Merged_PGS)

SBP_PGS_baseline_with_sodium <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
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
                         PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_DBP_PGS_full, data = Merged_PGS)

SBP_DBP_PGS_corr <- cor(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Scaled_DBP_PGS_full)

## Test whether adding SBP and DBP PGS improves models

SBP_PGS_baseline_with_DBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                  PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_SBP_PGS_full + Scaled_DBP_PGS_full, data = Merged_PGS)

anova(SBP_PGS_baseline, SBP_PGS_baseline_with_DBP, test = "F")

## Are the PGS correlated with any of the urinary measures

cor.test(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$Na_urine_millimolL)
cor.test(Merged_PGS$Scaled_DBP_PGS_full, Merged_PGS$Na_urine_millimolL)
cor.test(Merged_PGS$Scaled_SBP_PGS_full, Merged_PGS$K_urine_millimolL)
cor.test(Merged_PGS$Scaled_DBP_PGS_full, Merged_PGS$K_urine_millimolL)

## Is there a non-linear interaction between the PGS and the urinary measures: code the variables as continous as as high + low PES

## Scale Urinary measures to SD = 1 so more interpretable interaction

Merged_PGS$Scaled_K_full <- as.numeric(scale(Merged_PGS$K_urine_millimolL))
Merged_PGS$Scaled_Na_full <- as.numeric(scale(Merged_PGS$Na_urine_millimolL))

## SBP ~ Na
SBP_PGS_baseline_Na_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                         Assesment_centre + Scaled_SBP_PGS_full*Scaled_Na_full, data = Merged_PGS)

SBP_PGS_baseline_Na_int_all_cov_int <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + BP_meds*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                            PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full + PC3*Scaled_SBP_PGS_full + 
                                            PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC10*Scaled_SBP_PGS_full + PC11*Scaled_SBP_PGS_full +
                                            PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full +
                                           Assesment_centre*Scaled_SBP_PGS_full + Scaled_SBP_PGS_full*Scaled_Na_full, data = Merged_PGS)

## DBP ~ Na

DBP_PGS_baseline_Na_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                           PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_DBP_PGS_full*Scaled_Na_full, data = Merged_PGS)

## SBP ~ K

SBP_PGS_baseline_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                           Assesment_centre + Scaled_SBP_PGS_full*Scaled_K_full, data = Merged_PGS)


SBP_PGS_baseline_K_int_all_cov_int <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + BP_meds*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                            PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full + PC3*Scaled_SBP_PGS_full + 
                                            PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC10*Scaled_SBP_PGS_full + PC11*Scaled_SBP_PGS_full +
                                            PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full +
                                            Assesment_centre*Scaled_SBP_PGS_full + Scaled_SBP_PGS_full*Scaled_K_full, data = Merged_PGS)


## DBP ~ K

DBP_PGS_baseline_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                          PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                          PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_DBP_PGS_full*Scaled_K_full, data = Merged_PGS)



## SBP ~ Na:K

SBP_PGS_baseline_Na_K_int_simple_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                             PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                             PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                          Assesment_centre + Scaled_SBP_PGS_full*Na_K_ratio, data = Merged_PGS)

DBP_PGS_baseline_Na_K_int_simple_int <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                             PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                             PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                             Assesment_centre + Scaled_DBP_PGS_full*Na_K_ratio, data = Merged_PGS)

## Test if there is a effect of having high PGS and high Na:K ratio

Merged_PGS$High_PGS <- ifelse(Merged_PGS$Scaled_SBP_PGS_full > 1.278867933, "High PGS", "Comparator PGS (< 90th percentile)")

quantile(Merged_PGS$Na_K_ratio, probs = seq(.1, .9, by = .1))

Merged_PGS$High_Na_K <- ifelse(Merged_PGS$Na_K_ratio > 2.4700657, "High_Ratio", "Reamining_Ratio")

High_Na_K_cohort <- Merged_PGS %>% filter(High_Na_K == "High_Ratio")


High_in_both <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                     Assesment_centre +
                     PC1 + PC2 + PC3 + PC4 + PC5 + High_PGS, data = High_Na_K_cohort)

ggplot(data = High_Na_K_cohort, aes(x=High_PGS, y = Mean_SBP, fill = High_PGS)) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 142.1815, linetype = "dashed", colour = "darkred") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("> 90th percentile PGS","< 90th percentile PGS")) +
  ylab("Mean SBP (mmHg)") +
  xlab("SBP PGS") +
  ggtitle("Genome-wide PGS - SBP")

## Make quartiles of Na, K, and Na/K

Merged_PGS$Na_quartile <- ntile(Merged_PGS$Na_urine_millimolL, 4)
Merged_PGS$K_quartile <- ntile(Merged_PGS$K_urine_millimolL, 4)
Merged_PGS$Na_K_quartile <- ntile(Merged_PGS$Na_K_ratio, 4)


ggplot(data = Merged_PGS, aes(x=as.factor(Na_quartile), y = Mean_SBP, fill = as.factor(High_PGS))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "firebrick") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "dodgerblue")) +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4")) +
  ylab("Mean SBP (mmHg)") +
  xlab("Urinary sodium (millimol/L)") +
  ggtitle("Genome-wide PGS - SBP") +
  labs(fill = " ")

ggplot(data = Merged_PGS, aes(x=as.factor(K_quartile), y = Mean_SBP, fill = as.factor(High_PGS))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "firebrick") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "plum")) +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4")) +
  ylab("Mean SBP (mmHg)") +
  xlab("Urinary potassium (millimol/L)") +
  ggtitle("Genome-wide PGS - SBP") +
  labs(fill = " ")

ggplot(data = Merged_PGS, aes(x=as.factor(Na_K_quartile), y = Mean_SBP, fill = as.factor(High_PGS))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "firebrick") +
  theme_bw() +
  scale_fill_manual(values=c("#999999", "lightgoldenrod")) +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=c("Quartile 1","Quartile 2", "Quartile 3", "Quartile 4")) +
  ylab("Mean SBP (mmHg)") +
  xlab("Na:K") +
  ggtitle("Genome-wide PGS - SBP") +
  labs(fill = " ")

## Test slope of Na and K on SBP in high vs low PGS (< 10%)

quantile(Merged_PGS$Scaled_SBP_PGS_full, probs = seq(.1, .9, by = .1))

High_PGS_dat <- Merged_PGS %>% filter(High_PGS == "High PGS")
Low_PGS_dat <- Merged_PGS %>% filter(Scaled_SBP_PGS_full < -1.283519547)

## Na ~ SBP (high PGS)

SBP_Na_high_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                     Na_urine_millimolL, data = High_PGS_dat)

SBP_Na_low_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                        Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                        PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                        Na_urine_millimolL, data = Low_PGS_dat)

## Compare slopes of sodium using Clogg et al. 1995 method
## Z = (b1-b2)/(sqrt(se1^2+se2^2))
## High 0.023559   0.002465 
## Low 0.028354   0.002328

Z_Na_SBP_Clog <- (0.023559-0.028354)/sqrt(((0.002465)^2)+((0.002328)^2))

## Trend towards larger slope in low Z = -1.414226

pnorm(abs(Z_Na_SBP_Clog), lower.tail=FALSE) * 2

Na_plot_high_PGS <- ggplot(data = High_PGS_dat, aes(x = Na_urine_millimolL, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Urinary Sodium (millimol/L)") +
  theme(legend.position = "none") +
  ggtitle("> 90th percentile SBP PGS")

Na_plot_low_PGS <- ggplot(data = Low_PGS_dat, aes(x = Na_urine_millimolL, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE,colour = "red") +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Urinary Sodium (millimol/L)") +
  theme(legend.position = "none") +
  ggtitle("< 10th percentile SBP PGS")

ggarrange(Na_plot_high_PGS, Na_plot_low_PGS, nrow = 2)


## Try K 

SBP_K_high_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                        Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                        PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                        K_urine_millimolL, data = High_PGS_dat)

SBP_K_low_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                       Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                       PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                       K_urine_millimolL, data = Low_PGS_dat)

## High -0.038084   0.003058
## Low -0.030711   0.002849

Z_K_SBP_Clog <- (-0.038084--0.030711)/sqrt((0.003058^2)+(0.002849^2))

## Trend towards Nominal sig in same direction -1.764089

pnorm(abs(Z_K_SBP_Clog), lower.tail=FALSE) * 2
