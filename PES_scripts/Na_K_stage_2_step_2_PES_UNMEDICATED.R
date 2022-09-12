#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2 - PES impact on UKBB, any interplay with Na, K or the ratio - UNMEDICATED ONLY

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

## Read in files, set relevant variables as categorical

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

## Scale PES - unmedicated

Merged_PES[, c(2:7)] <-  lapply(Merged_PES[, c(2:7)], function(x) c(scale(x)))

List_of_SBP_PES <- as.list(colnames(Merged_PES %>% select(contains("SBP")) %>% select(contains("PES"))))
List_of_DBP_PES <- as.list(colnames(Merged_PES %>% select(contains("DBP")) %>% select(contains("PES"))))

List_PES <- c(List_of_DBP_PES, List_of_SBP_PES)

List_traits <- as.list(c("Mean_DBP", "Mean_DBP", "Mean_DBP", "Mean_SBP", "Mean_SBP", "Mean_SBP"))

## Test without sodium or potassium adjustment

PES_BP_assoc <- function(score, BP, df) {
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}

## Run list of scores

Run_PES_BP <- mapply(PES_BP_assoc, score = List_PES, BP =  List_traits,
                     MoreArgs = list(Merged_PES),
                     SIMPLIFY = T)

## Extract relevant data

Run_BP_extract <- apply(Run_PES_BP, 2, function(x) return(as.data.frame(x$coefficients)[53, 1:4]))

Run_BP_results <- data.frame()

for (i in 1:length(Run_BP_extract)) {
  Run_BP_results <- rbind(Run_BP_results, Run_BP_extract[[i]])
}

Run_BP_results$PES <- unlist(List_PES)
Run_BP_results$BP_trait <- unlist(List_traits)
Run_BP_results$PGS <- "PGS unadjusted"

## TEST impact of PGS adjustment

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")

## Test effect of SBP PGS on SBP and DBP PGS on DBP
## Compare effect size of urinary measures per SD to PGS per SD

Merged_PGS_PES$Scaled_SBP_PGS_full <- as.numeric(scale(Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))
Merged_PGS_PES$Scaled_DBP_PGS_full <- as.numeric(scale(Merged_PGS_PES$DBP_PGS_Pt_0_005_AVG))

List_PGS <- as.list(c("Scaled_DBP_PGS_full", "Scaled_DBP_PGS_full", "Scaled_DBP_PGS_full", "Scaled_SBP_PGS_full",
                      "Scaled_SBP_PGS_full", "Scaled_SBP_PGS_full"))

PES_BP_assoc_PGS_adjust <- function(score, BP, PGS, df) {
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", PGS, "+ ", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}

## Run list of scores

Run_PES_BP_PGS_adjust <- mapply(PES_BP_assoc_PGS_adjust, score = List_PES, BP =  List_traits, PGS = List_PGS,
                                MoreArgs = list(Merged_PGS_PES),
                                SIMPLIFY = T)

## Extract relevant data

Run_BP_extract_PGS_adjust <- apply(Run_PES_BP_PGS_adjust, 2, function(x) return(as.data.frame(x$coefficients)[54, 1:4]))

Run_BP_results_PGS_adjust <- data.frame()

for (i in 1:length(Run_BP_extract_PGS_adjust)) {
  Run_BP_results_PGS_adjust <- rbind(Run_BP_results_PGS_adjust, Run_BP_extract_PGS_adjust[[i]])
}

Run_BP_results_PGS_adjust$PES <- unlist(List_PES)
Run_BP_results_PGS_adjust$BP_trait <- unlist(List_traits)
Run_BP_results_PGS_adjust$PGS <- "PGS adjusted"

PES_BP_results_unmed_cohort <- rbind(Run_BP_results, Run_BP_results_PGS_adjust)

write.table(PES_BP_results_unmed_cohort, file = "UKBB_PES/PES_BP_assoc_and_int/PES_assoc_UNMEDICATED_cohort_with_and_without_PGS_adjustment.txt",
            sep = "\t", row.names = F, quote = F)


Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_urine_millimolL))
Merged_PGS_PES$Scaled_K <- as.numeric(scale(Merged_PGS_PES$K_urine_millimolL))
Merged_PGS_PES$Binary_ratio <- ifelse(Merged_PGS_PES$Na_K_ratio > 1.33, "High", "Low")

## Interaction analyses

PES_Na_K_int <- function(score, Urine, BP, df) {
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", Urine, "*", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}


Run_int <- mapply(PES_Na_K_int, score = List_PES, BP =  List_traits, Urine = "Scaled_Na",
                  MoreArgs = list(Merged_PGS_PES),
                  SIMPLIFY = T)

Run_int_K <- mapply(PES_Na_K_int, score = List_PES, BP =  List_traits, Urine = "Scaled_K",
                    MoreArgs = list(Merged_PGS_PES),
                    SIMPLIFY = T)
Extract_int_Na <- apply(Run_int, 2, function(x) return(as.data.frame(x$coefficients)[55, 1:4]))

Extract_int_K <- apply(Run_int_K, 2, function(x) return(as.data.frame(x$coefficients)[55, 1:4]))

## Test with GxC terms and CxE terms

Renal_int <- lm(Mean_DBP ~ Sex*Renal_DBP_PES_AVG + Age*Renal_DBP_PES_AVG + Age2*Renal_DBP_PES_AVG + Assesment_centre_month*Renal_DBP_PES_AVG + 
                                     Assesment_centre*Renal_DBP_PES_AVG + PC1*Renal_DBP_PES_AVG + PC2*Renal_DBP_PES_AVG + PC3*Renal_DBP_PES_AVG + PC4*Renal_DBP_PES_AVG + PC5*Renal_DBP_PES_AVG + PC10*Renal_DBP_PES_AVG + PC11*Renal_DBP_PES_AVG +
                                     PC12*Renal_DBP_PES_AVG + PC13*Renal_DBP_PES_AVG + PC14*Renal_DBP_PES_AVG + PC15*Renal_DBP_PES_AVG + PC16*Renal_DBP_PES_AVG + PC17*Renal_DBP_PES_AVG + PC18*Renal_DBP_PES_AVG + PC19*Renal_DBP_PES_AVG + PC20*Renal_DBP_PES_AVG + Scaled_K*Renal_DBP_PES_AVG,
                data = Merged_PGS_PES)

Renal_int_G_C_C_E <- lm(Mean_DBP ~ Sex*Renal_DBP_PES_AVG + Age*Renal_DBP_PES_AVG + Age2*Renal_DBP_PES_AVG + Assesment_centre_month*Renal_DBP_PES_AVG + 
                                     Assesment_centre*Renal_DBP_PES_AVG + PC1*Renal_DBP_PES_AVG + PC2*Renal_DBP_PES_AVG + PC3*Renal_DBP_PES_AVG + PC4*Renal_DBP_PES_AVG + PC5*Renal_DBP_PES_AVG + PC10*Renal_DBP_PES_AVG + PC11*Renal_DBP_PES_AVG +
                                     PC12*Renal_DBP_PES_AVG + PC13*Renal_DBP_PES_AVG + PC14*Renal_DBP_PES_AVG + PC15*Renal_DBP_PES_AVG + PC16*Renal_DBP_PES_AVG + PC17*Renal_DBP_PES_AVG + PC18*Renal_DBP_PES_AVG + PC19*Renal_DBP_PES_AVG + PC20*Renal_DBP_PES_AVG + Scaled_K*Renal_DBP_PES_AVG +
                          Sex*Scaled_K + Age*Scaled_K + Age2*Scaled_K + Assesment_centre_month*Scaled_K + 
                          Assesment_centre*Scaled_K + PC1*Scaled_K + PC2*Scaled_K + PC3*Scaled_K + PC4*Scaled_K + PC5*Scaled_K + PC10*Scaled_K + PC11*Scaled_K +
                          PC12*Scaled_K + PC13*Scaled_K + PC14*Scaled_K + PC15*Scaled_K + PC16*Scaled_K + PC17*Scaled_K + PC18*Scaled_K + PC19*Scaled_K + PC20*Scaled_K,
                          
                data = Merged_PGS_PES)

Transport_int <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                     Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG+ PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC10*Transport_SBP_PES_AVG+ PC11*Transport_SBP_PES_AVG +
                                     PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + Scaled_Na*Transport_SBP_PES_AVG +
                      ,
                data = Merged_PGS_PES)

Transport_int_G_C_C_E <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                                     Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG+ 
                              PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC10*Transport_SBP_PES_AVG+ PC11*Transport_SBP_PES_AVG +
                                     PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + 
                              PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                              Scaled_Na*Transport_SBP_PES_AVG +
                              Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                              Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                              PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                    data = Merged_PGS_PES)

Transport_int_G_C_C_E_PGS <- lm(Mean_SBP ~ Sex*Scaled_SBP_PGS_full + Age*Scaled_SBP_PGS_full + Age2*Scaled_SBP_PGS_full + Assesment_centre_month*Scaled_SBP_PGS_full + 
                                  Assesment_centre*Scaled_SBP_PGS_full + PC1*Scaled_SBP_PGS_full + PC2*Scaled_SBP_PGS_full+ 
                                  PC3*Scaled_SBP_PGS_full + PC4*Scaled_SBP_PGS_full + PC5*Scaled_SBP_PGS_full + PC10*Scaled_SBP_PGS_full+ PC11*Scaled_SBP_PGS_full +
                                  PC12*Scaled_SBP_PGS_full + PC13*Scaled_SBP_PGS_full + PC14*Scaled_SBP_PGS_full + PC15*Scaled_SBP_PGS_full + 
                                  PC16*Scaled_SBP_PGS_full + PC17*Scaled_SBP_PGS_full + PC18*Scaled_SBP_PGS_full + PC19*Scaled_SBP_PGS_full + PC20*Scaled_SBP_PGS_full + 
                                  Scaled_Na*Scaled_SBP_PGS_full +
                                  Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                                  Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                                  PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                                data = Merged_PGS_PES)

Transport_int_null <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                           Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                           PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Scaled_Na + Transport_SBP_PES_AVG,
                         data = Merged_PGS_PES)

## Transport interaction plot

interact_plot(Transport_int, pred = "Scaled_Na", modx = "Transport_SBP_PES_AVG", modx.values = c(-1.283521369, 0, 1.279041346),
          interval = T)

## Try another plot

x <- Merged_PGS_PES$Transport_SBP_PES_AVG
y <- Merged_PGS_PES$Scaled_SBP_PGS_full

Merged_PGS_PES$PES <-
  case_when(x > mean(x)+sd(x) ~ "High",
            x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "Average",
            x < mean(x)-sd(x) ~ "Low")

Merged_PGS_PES$PGS <-
  case_when(y > mean(y)+sd(y) ~ "High",
            y < mean(y)+sd(y) & y > mean(y)-sd(y) ~ "Average",
            y < mean(y)-sd(y) ~ "Low")

Plot_2 <- Merged_PGS_PES %>% 
  ggplot() +
  aes(x = Scaled_Na, y = Mean_SBP, group = PES, color = PES) +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Spot Urinary Na (SD)") +
  ylab("SBP (mmHg)") +
  ggtitle("Main additive effects")

## Add in full multiple regression effects

Merged_PGS_PES$Trans_fit <- predict(Transport_int_G_C_C_E)
Merged_PGS_PES$Main_fit <- predict(Transport_int_null)

Plot1 <- Merged_PGS_PES %>% 
  ggplot() +
  aes(x = Mean_SBP, y = Trans_fit, group = PES, color = PES) +
  geom_smooth(method = "lm", se = TRUE, aes(y=Trans_fit)) +
  theme_bw() +
  xlab("SBP (mmHg)") +
  ylab("Pr") +
  ggtitle("Predicted SBP (G, E, C, GxE, GxC, ExC)")

ggarrange(Plot_2, Plot1, ncol = 2,
          legend = "bottom", common.legend = T)

Merged_PGS_PES %>% 
  ggplot() +
  aes(x = Scaled_Na, y = Mean_SBP, group = PGS, color = PGS) +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Spot Urinary Na (SD)") +
  ylab("SBP (mmHg)") +
  ggtitle("SBP ~ Urinary Na (stratified by PGS)")


sqrt(mean(Transport_int_null$residuals^2))
sqrt(mean(Transport_int$residuals^2))

## Make plots for extreme quantiles

quantile(Merged_PGS_PES$Transport_SBP_PES_AVG, probs = seq(.1, .9, by = .1))

High_PES_dat <- Merged_PGS_PES %>% filter(Transport_SBP_PES_AVG > 1.276391e+00)
High_PES_dat$Scaled_Na_new <- as.numeric(scale(High_PES_dat$Na_urine_millimolL))
Low_PES_dat <- Merged_PGS_PES %>% filter(Transport_SBP_PES_AVG < -1.279876e+00)
Low_PES_dat$Scaled_Na_new <- as.numeric(scale(Low_PES_dat$Na_urine_millimolL))

SBP_Na_high_PES <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                       Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                       PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                       Scaled_Na_new, data = High_PES_dat)

SBP_Na_low_PES <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                      Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                      Scaled_Na_new, data = Low_PES_dat)

Low_PES_dat$fitted <- predict(SBP_Na_low_PES)
High_PES_dat$fitted <- predict(SBP_Na_high_PES)

Na_plot_high_PES <- ggplot(data = High_PES_dat, aes(x = Scaled_Na_new, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE, aes(y=fitted))
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Sodium (SD)") +
  theme(legend.position = "none") +
  ggtitle("> 90th percentile Transport PES") 

Na_plot_low_PES <- ggplot(data = Low_PES_dat, aes(x = Scaled_Na_new, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE,colour = "red", aes(y=fitted)) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Sodium (SD)") +
  theme(legend.position = "none") +
  ggtitle("< 10th percentile Transport PES")

ggarrange(Na_plot_high_PES, Na_plot_low_PES, nrow = 2)

quantile(Merged_PGS_PES$Scaled_SBP_PGS_full, probs = seq(.1, .9, by = .1))

High_PGS_dat <- Merged_PGS_PES %>% filter(Scaled_SBP_PGS_full > 1.279041346)
High_PGS_dat$Scaled_Na_new <- as.numeric(scale(High_PGS_dat$Na_urine_millimolL))
Low_PGS_dat <- Merged_PGS_PES %>% filter(Scaled_SBP_PGS_full < -1.283521369)
Low_PGS_dat$Scaled_Na_new <- as.numeric(scale(Low_PGS_dat$Na_urine_millimolL))

SBP_Na_high_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                       Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                       PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                       Scaled_Na_new, data = High_PGS_dat)

SBP_Na_low_PGS <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                      Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                      Scaled_Na_new, data = Low_PGS_dat)

Low_PGS_dat$fitted <- predict(SBP_Na_low_PGS)
High_PGS_dat$fitted <- predict(SBP_Na_high_PGS)


Na_plot_high_PGS <- ggplot(data = High_PGS_dat, aes(x = Scaled_Na, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE, aes(y=fitted)) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Sodium Potassium (SD)") +
  theme(legend.position = "none") +
  ggtitle("> 90th percentile genome-wide PGS")

Na_plot_low_PGS <- ggplot(data = Low_PGS_dat, aes(x = Scaled_Na, y = Mean_SBP)) + 
  geom_smooth(method = "lm", se = TRUE,colour = "red", aes(y=fitted)) +
  ylab("Mean SBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Sodium (SD)") +
  theme(legend.position = "none") +
  ggtitle("< 10th percentile genome-wide PGS")

ggarrange(Na_plot_high_PGS, Na_plot_low_PGS, nrow = 2)

## Plot effect on SBP per Na decile in low and high


