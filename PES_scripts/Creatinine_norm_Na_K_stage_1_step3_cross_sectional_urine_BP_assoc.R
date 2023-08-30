#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 1: Step 3 - cross-sectional relationship between Na, K, and Na:K with BP

## Normalised to urinary creatinine

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

## Read in data generated from previous script

Subcohort_1_dat <- fread("Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt")



Subcohort_1_dat$BP_meds <- as.factor(Subcohort_1_dat$BP_meds)
Subcohort_1_dat$Assesment_centre_month <- as.factor(Subcohort_1_dat$Assesment_centre_month)
Subcohort_1_dat$Assesment_centre <- as.factor(Subcohort_1_dat$Assesment_centre)

## Make Na:Cr and K:Cr

Subcohort_1_dat$Na_Cr <- Subcohort_1_dat$Na_urine_millimolL/Subcohort_1_dat$Creatinine_urine_micromolL
Subcohort_1_dat$K_Cr <- Subcohort_1_dat$K_urine_millimolL/Subcohort_1_dat$Creatinine_urine_micromolL

table1( ~ factor(Sex) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL + Na_Cr + K_Cr | factor(`BP Medication`),
        data = Subcohort_1_dat)

Unmed_cohort <- Subcohort_1_dat %>% filter(BP_meds == 0)

## Write new df

write.table(Subcohort_1_dat, file="Creatinine_adjusted_subcohort_1.txt",
            sep = "\t", row.names = F, quote = F)

## Plots comparing Na:Cr and K:Cr to raw values - full and unmedicated

Subcohort_1_dat$Decile_Na_Cr <- as.factor(ntile(Subcohort_1_dat$Na_Cr, 10))
Subcohort_1_dat$Decile_K_Cr <- as.factor(ntile(Subcohort_1_dat$K_Cr, 10))

P1 <- ggplot(Subcohort_1_dat, aes(x=Decile_Na_Cr, y=Na_urine_millimolL, fill=Sex)) +
  geom_boxplot() +
  theme_bw() +
  xlab("UNa:UCr (decile)") +
  ylab("Raw UNa (mmol/L)") +
  ylim(0, 382)

P2 <- ggplot(Subcohort_1_dat, aes(x=Decile_K_Cr, y=K_urine_millimolL, fill=Sex)) +
  geom_boxplot() +
  theme_bw() +
  xlab("UK:UCr (decile)") +
  ylab("Raw UK (mmol/L)") +
  ylim(0, 220)

Full_compare <- ggarrange(P1, P2, common.legend = T, legend = "none")

Unmed_cohort$Decile_Na_Cr <- as.factor(ntile(Unmed_cohort$Na_Cr, 10))
Unmed_cohort$Decile_K_Cr <- as.factor(ntile(Unmed_cohort$K_Cr, 10))

P3 <- ggplot(Unmed_cohort, aes(x=Decile_Na_Cr, y=Na_urine_millimolL, fill=Sex)) +
  geom_boxplot() +
  theme_bw() +
  xlab("UNa:UCr (decile)") +
  ylab("Raw UNa (mmol/L)") +
  ylim(0, 382)

P4 <- ggplot(Unmed_cohort, aes(x=Decile_K_Cr, y=K_urine_millimolL, fill=Sex)) +
  geom_boxplot() +
  theme_bw() +
  xlab("UK:UCr (decile)") +
  ylab("Raw UK (mmol/L)") +
  ylim(0, 220)

Unmed_compare <- ggarrange(P3, P4, legend = "bottom", common.legend = T)

Raw_compare <- ggarrange(Full_compare, Unmed_compare, nrow = 2, common.legend = T, legend = "bottom")

annotate_figure(Raw_compare,
                top="Full cohort", "Unmedicated")

## Test assoc between raw and ln transformed UNa:UCr and UK:UCr with SBP and DBP separately.

## Adjust for age, age2, sex, assesesment centre, month of visit, blood pressure meds

## SBP ##

Scaled_Na_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + scale(Na_Cr), data = Subcohort_1_dat)

Scaled_K_SBP_ALL <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + scale(K_Cr), data = Subcohort_1_dat)

## DBP ##

Scaled_Na_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                          Assesment_centre + scale(Na_Cr), data = Subcohort_1_dat)

Scaled_K_DBP_ALL <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                         Assesment_centre + scale(K_Cr), data = Subcohort_1_dat)

## Plot SBP per quintile of UNa:UCr

Subcohort_1_dat$Quint_Na_Cr <- as.factor(ntile(Subcohort_1_dat$Na_Cr, 5))
Subcohort_1_dat$Quint_K_Cr <- as.factor(ntile(Subcohort_1_dat$K_Cr, 5))


Ratio_SBP_Na <- ggplot(data = Subcohort_1_dat, aes(x=Quint_Na_Cr, y = Mean_SBP, fill = Quint_Na_Cr)) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Mean SBP (mmHg)") +
  xlab("UNa:UCr (quintile)") +
  scale_fill_brewer(palette = "Pastel1")

Ratio_SBP_K <- ggplot(data = Subcohort_1_dat, aes(x=Quint_K_Cr, y = Mean_SBP, fill = Quint_K_Cr)) +
  geom_violin() +
  geom_boxplot(width=0.4, fill="white") +
  geom_hline(yintercept = 138.228, linetype = "dashed", colour = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Mean SBP (mmHg)") +
  xlab("UNa:UK (quintile)") +
  scale_fill_brewer(palette = "Pastel1")

ggarrange(Ratio_SBP_Na, Ratio_SBP_K)

## Plot direct relationship




## Adjust for BMI

BMI_adj_Na_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                    Assesment_centre + BMI + scale(Na_Cr), data = Subcohort_1_dat)

BMI_adj_Na_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                    Assesment_centre + BMI + scale(Na_Cr), data = Subcohort_1_dat)

BMI_adj_K_SBP <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + scale(K_Cr), data = Subcohort_1_dat)

BMI_adj_K_SBP_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                      Assesment_centre + scale(BMI)*scale(K_Cr), data = Subcohort_1_dat)

BMI_adj_Na_SBP_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                          Assesment_centre + scale(BMI)*scale(Na_Cr), data = Subcohort_1_dat)


BMI_adj_K_DBP <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month +
                       Assesment_centre + BMI + scale(K_Cr), data = Subcohort_1_dat)


## Corr between BMI and UNa:UCr and UK:UCr

BMI_Na_Cr <- lm(BMI ~ Sex + Age + Age2 + BP_meds + scale(Na_Cr), data = Subcohort_1_dat)

BMI_K_Cr <- lm(BMI ~ Sex + Age + Age2 + BP_meds + scale(K_Cr), data = Subcohort_1_dat)

## Stratify by meds

## Remove medicated individuals

Raw_Na_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                       Assesment_centre + scale(Na_Cr), data = Unmed_cohort)

Raw_K_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                      Assesment_centre + scale(K_Cr), data = Unmed_cohort)

## DBP ##

Raw_Na_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                       Assesment_centre + scale(Na_Cr), data = Unmed_cohort)

Raw_K_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                      Assesment_centre + scale(K_Cr), data = Unmed_cohort)

## BMI adjustment

BMI_Raw_Na_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                         Assesment_centre + scale(Na_Cr) + BMI, data = Unmed_cohort)

BMI_Raw_K_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                        Assesment_centre + scale(K_Cr) + BMI, data = Unmed_cohort)

Int_BMI_Raw_K_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                            Assesment_centre + scale(K_Cr)*scale(BMI), data = Unmed_cohort)

Int_BMI_urinary_K_SBP_UNMED <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month +
                                Assesment_centre + scale(K_urine_millimolL)*scale(BMI), data = Unmed_cohort)

## DBP ##

BMI_Raw_Na_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                         Assesment_centre + scale(Na_Cr) + BMI, data = Unmed_cohort)

BMI_Raw_K_DBP_UNMED <- lm(Mean_DBP ~ Sex + Age + Age2 + Assesment_centre_month +
                        Assesment_centre + scale(K_Cr) + BMI, data = Unmed_cohort)


## Test straight linear
## Quadratic
## Spline

## Creatinine to sodium

library(mgcv)

fit.1 <- lm(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
                      Assesment_centre + Creatinine_urine_micromolL, data = Unmed_cohort)

fit.2 <- lm(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
              Assesment_centre + Creatinine_urine_micromolL + I(Creatinine_urine_micromolL^2), data = Unmed_cohort)

fit.3 <- gam(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
              Assesment_centre + s(Creatinine_urine_micromolL), method = "REML", data = Unmed_cohort)




p1 <- Unmed_cohort %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "steelblue"
  )+
  labs(title = "GLM") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")

p2 <- Unmed_cohort %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    color = "steelblue"
  )+
  labs(title = "Quadratic effects") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")

p3 <- Unmed_cohort %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    color = "steelblue"
  )+
  labs(title = "Spline Model") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")

ggarrange(p1,p2,p3, nrow = 1)

mean(Unmed_cohort$Creatinine_urine_micromolL) - (5*sd(Unmed_cohort$Creatinine_urine_micromolL))
mean(Unmed_cohort$Creatinine_urine_micromolL) + (5*sd(Unmed_cohort$Creatinine_urine_micromolL))

Winsorized <- Unmed_cohort %>% filter(Creatinine_urine_micromolL < 36869.75)

fit.4 <- lm(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
              Assesment_centre + Creatinine_urine_micromolL, data = Winsorized)

fit.5 <- lm(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
              Assesment_centre + Creatinine_urine_micromolL + I(Creatinine_urine_micromolL^2), data = Winsorized)

fit.6 <- gam(Na_urine_millimolL ~ Sex + Age + Age2 + Assesment_centre_month +
               Assesment_centre + s(Creatinine_urine_micromolL), method = "REML", data = Winsorized)

p4 <- Winsorized %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "steelblue"
  )+
  labs(title = "GLM") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  geom_vline(xintercept = 14306.84, lty="dashed", colour = "purple") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")

p5 <- Winsorized %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    color = "steelblue"
  )+
  labs(title = "Quadratic effects") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  geom_vline(xintercept = 14306.84, lty="dashed", colour = "purple") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")

p6 <- Winsorized %>% 
  ggplot(aes(Creatinine_urine_micromolL,Na_urine_millimolL))+
  geom_point(alpha = 0.6)+
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    color = "steelblue"
  )+
  labs(title = "Spline Model") +
  geom_vline(xintercept = 8635.34, lty="dashed", colour = "red") +
  geom_vline(xintercept = 14306.84, lty="dashed", colour = "purple") +
  theme_bw() +
  xlab("Urinary creatinine (micromol/L)") +
  ylab("Urinary sodium (millimol/L)")


ggarrange(p4, p5, p6, nrow = 1)
