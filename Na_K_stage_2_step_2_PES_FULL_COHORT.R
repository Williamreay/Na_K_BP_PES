#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2 - PES impact on UKBB, any interplay with Na, K or the ratio

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

## Read in PES

Na_K_PES <- fread("UKBB_PES/PES_Na_K_SBP_DBP.txt")

Colnames_PES <- make.unique(names(Na_K_PES))

colnames(Na_K_PES) <- Colnames_PES

Na_K_PES <- Na_K_PES %>% select(IID, contains("PES"))

Merged_PES <- merge(Na_K_PES, Subcohort_1_dat, by = "IID")

## Scale PES - medicated

Merged_PES[, c(2:7)] <-  lapply(Merged_PES[, c(2:7)], function(x) c(scale(x)))

List_of_SBP_PES <- as.list(colnames(Merged_PES %>% select(contains("SBP")) %>% select(contains("PES"))))
List_of_DBP_PES <- as.list(colnames(Merged_PES %>% select(contains("DBP")) %>% select(contains("PES"))))

List_PES <- c(List_of_DBP_PES, List_of_SBP_PES)

List_traits <- as.list(c("Mean_DBP", "Mean_DBP", "Mean_DBP", "Mean_SBP", "Mean_SBP", "Mean_SBP"))

## Test without sodium or potassium adjustment

PES_BP_assoc <- function(score, BP, df) {
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
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

Run_BP_extract <- apply(Run_PES_BP, 2, function(x) return(as.data.frame(x$coefficients)[54, 1:4]))

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
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
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

Run_BP_extract_PGS_adjust <- apply(Run_PES_BP_PGS_adjust, 2, function(x) return(as.data.frame(x$coefficients)[55, 1:4]))

Run_BP_results_PGS_adjust <- data.frame()

for (i in 1:length(Run_BP_extract_PGS_adjust)) {
  Run_BP_results_PGS_adjust <- rbind(Run_BP_results_PGS_adjust, Run_BP_extract_PGS_adjust[[i]])
}

Run_BP_results_PGS_adjust$PES <- unlist(List_PES)
Run_BP_results_PGS_adjust$BP_trait <- unlist(List_traits)
Run_BP_results_PGS_adjust$PGS <- "PGS adjusted"

PES_BP_results_full_cohort <- rbind(Run_BP_results, Run_BP_results_PGS_adjust)

write.table(PES_BP_results_full_cohort, file = "UKBB_PES/PES_BP_assoc_and_int/PES_assoc_FULL_cohort_with_and_without_PGS_adjustment.txt",
            sep = "\t", row.names = F, quote = F)

PES_BP_results_full_cohort$L_UI <- PES_BP_results_full_cohort$Estimate - (1.96*PES_BP_results_full_cohort$`Std. Error`)
PES_BP_results_full_cohort$U_UI <- PES_BP_results_full_cohort$Estimate + (1.96*PES_BP_results_full_cohort$`Std. Error`)

PES_BP_results_full_cohort$PES_name <- c("Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport")

SBP_results <- PES_BP_results_full_cohort %>% filter(BP_trait == "Mean_SBP")
DBP_results <- PES_BP_results_full_cohort %>% filter(BP_trait == "Mean_DBP")

SBP_full <- ggplot(SBP_results, aes(x=PES_name, y=Estimate, ymin=L_UI, ymax=U_UI,colour = PES_name)) +
  geom_pointrange() +
  facet_wrap(~PGS, strip.position="top",nrow=2, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position = "bottom") +
  ylab("Effect per SD increase in PES on SBP (mmHg)") +
  theme(legend.position = "none") +
  xlab("PES score") +
  ylim(0, 0.6) +
  ggtitle("SBP")

DBP_full <- ggplot(DBP_results, aes(x=PES_name, y=Estimate, ymin=L_UI, ymax=U_UI,colour = PES_name)) +
  geom_pointrange() +
  facet_wrap(~PGS, strip.position="top",nrow=2, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  theme(legend.position = "bottom") +
  ylab("Effect per SD increase in PES on DBP (mmHg)") +
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y =  element_blank()) +
  xlab(" ") +
  ylim(0, 0.6) +
  ggtitle("DBP")

ggarrange(SBP_full, DBP_full)

## Assoc of PES with measured Na and K

Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_urine_millimolL))
Merged_PGS_PES$Scaled_K <- as.numeric(scale(Merged_PGS_PES$K_urine_millimolL))
Merged_PGS_PES$Binary_ratio <- ifelse(Merged_PGS_PES$Na_K_ratio > 1.33, "High", "Low")

cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Absorption_DBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Absorption_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Renal_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Renal_DBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Transport_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_Na, Merged_PGS_PES$Transport_DBP_PES_AVG)

cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Absorption_DBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Absorption_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Renal_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Renal_DBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Transport_SBP_PES_AVG)
cor.test(Merged_PGS_PES$Scaled_K, Merged_PGS_PES$Transport_DBP_PES_AVG)

## Nom positive correlation with the absorption PES without any covariate adjustment 

## Test with covariate adjustment

PES_Na_K_corr <- function(score, Urine, df) {
  fmla <- as.formula(paste(Urine, "~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}

Na_corr <- sapply(List_PES, PES_Na_K_corr, Urine = "Scaled_Na", df = Merged_PGS_PES)

Na_results_corr <- apply(Na_corr, 2, function(x) return(as.data.frame(x$coefficients)[54, 1:4]))

K_corr <- sapply(List_PES, PES_Na_K_corr, Urine = "Scaled_K", df = Merged_PGS_PES)

K_results_corr <- apply(K_corr, 2, function(x) return(as.data.frame(x$coefficients)[54, 1:4]))

## Interaction models - don't test without GxC

PES_Na_K_int <- function(score, Urine, BP, df) {
  fmla <- as.formula(paste(BP, "~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
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

Run_int_ratio <- mapply(PES_Na_K_int, score = List_PES, BP =  List_traits, Urine = "Binary_ratio",
                    MoreArgs = list(Merged_PGS_PES),
                    SIMPLIFY = T)

## Extract relevant data

Extract_int_Na <- apply(Run_int, 2, function(x) return(as.data.frame(x$coefficients)[56, 1:4]))

Extract_int_K <- apply(Run_int_K, 2, function(x) return(as.data.frame(x$coefficients)[56, 1:4]))

Extract_int_Ratio <- apply(Run_int_ratio, 2, function(x) return(as.data.frame(x$coefficients)[56, 1:4]))

## Nominal interaction between Renal DBP PES and scaled K on DBP

Renal_int_G_by_C <- lm(Mean_DBP ~ Sex*Renal_DBP_PES_AVG + Age*Renal_DBP_PES_AVG + Age2*Renal_DBP_PES_AVG + BP_meds*Renal_DBP_PES_AVG + Assesment_centre_month*Renal_DBP_PES_AVG + 
                                     Assesment_centre*Renal_DBP_PES_AVG + PC1*Renal_DBP_PES_AVG + PC2*Renal_DBP_PES_AVG + PC3*Renal_DBP_PES_AVG + PC4*Renal_DBP_PES_AVG + PC5*Renal_DBP_PES_AVG + PC10*Renal_DBP_PES_AVG + PC11*Renal_DBP_PES_AVG +
                                     PC12*Renal_DBP_PES_AVG + PC13*Renal_DBP_PES_AVG + PC14*Renal_DBP_PES_AVG + PC15*Renal_DBP_PES_AVG + PC16*Renal_DBP_PES_AVG + PC17*Renal_DBP_PES_AVG + PC18*Renal_DBP_PES_AVG + PC19*Renal_DBP_PES_AVG + PC20*Renal_DBP_PES_AVG + Scaled_K*Renal_DBP_PES_AVG,
                data = Merged_PGS_PES)


Renal_int_G_by_C_and_c_by_E <- lm(Mean_DBP ~ Sex*Renal_DBP_PES_AVG + Age*Renal_DBP_PES_AVG + Age2*Renal_DBP_PES_AVG + BP_meds*Renal_DBP_PES_AVG + Assesment_centre_month*Renal_DBP_PES_AVG + 
                         Assesment_centre*Renal_DBP_PES_AVG + PC1*Renal_DBP_PES_AVG + PC2*Renal_DBP_PES_AVG + PC3*Renal_DBP_PES_AVG + PC4*Renal_DBP_PES_AVG + PC5*Renal_DBP_PES_AVG + PC10*Renal_DBP_PES_AVG + PC11*Renal_DBP_PES_AVG +
                         PC12*Renal_DBP_PES_AVG + PC13*Renal_DBP_PES_AVG + PC14*Renal_DBP_PES_AVG + PC15*Renal_DBP_PES_AVG + PC16*Renal_DBP_PES_AVG + PC17*Renal_DBP_PES_AVG + PC18*Renal_DBP_PES_AVG + PC19*Renal_DBP_PES_AVG + PC20*Renal_DBP_PES_AVG + Scaled_K*Renal_DBP_PES_AVG +
                           Sex*Scaled_K + Age*Scaled_K + Age2*Renal_DBP_PES_AVG + BP_meds*Scaled_K + Assesment_centre_month*Scaled_K + 
                           Assesment_centre*Scaled_K + PC1*Scaled_K + PC2*Scaled_K + PC3*Scaled_K + PC4*Scaled_K + PC5*Scaled_K + PC10*Scaled_K + PC11*Scaled_K +
                           PC12*Scaled_K+ PC13*Scaled_K + PC14*Scaled_K + PC15*Scaled_K + PC16*Scaled_K + PC17*Scaled_K + PC18*Scaled_K + PC19*Scaled_K + PC20*Scaled_K,
                       data = Merged_PGS_PES)

plot_model(Renal_int, type = "pred", terms = c("Scaled_K", "Renal_DBP_PES_AVG"))

quantile(Merged_PGS_PES$Renal_DBP_PES_AVG, probs = seq(.1, .9, by = .1))

High_PES_dat <- Merged_PGS_PES %>% filter(Renal_DBP_PES_AVG > 1.2786882298)
High_PES_dat$Scaled_K_new <- as.numeric(scale(High_PES_dat$K_urine_millimolL))
Low_PES_dat <- Merged_PGS_PES %>% filter(Renal_DBP_PES_AVG < -1.2793258672)
Low_PES_dat$Scaled_K_new <- as.numeric(scale(Low_PES_dat$K_urine_millimolL))

DBP_K_high_PES <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                       Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                       PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                       Scaled_K_new, data = High_PES_dat)

DBP_K_low_PES <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                      Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                      Scaled_K_new, data = Low_PES_dat)

Low_PES_dat$fitted <- predict(DBP_K_low_PES)
High_PES_dat$fitted <- predict(DBP_K_high_PES)


K_plot_high_PES <- ggplot(data = High_PES_dat, aes(x = Scaled_K_new, y = Mean_DBP)) + 
  geom_smooth(method = "lm", se = TRUE, aes(y=fitted)) +
  ylab("Mean DBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Potassium (SD)") +
  theme(legend.position = "none") +
  ggtitle("> 90th percentile Renal PES")

K_plot_low_PES <- ggplot(data = Low_PES_dat, aes(x = Scaled_K_new, y = Mean_DBP)) + 
  geom_smooth(method = "lm", se = TRUE,colour = "red", aes(y=fitted)) +
  ylab("Mean DBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Potassium (SD)") +
  theme(legend.position = "none") +
  ggtitle("< 10th percentile Renal PES")

ggarrange(K_plot_high_PES, K_plot_low_PES, nrow = 2)

quantile(Merged_PGS_PES$Scaled_DBP_PGS_full, probs = seq(.1, .9, by = .1))

High_PGS_dat <- Merged_PGS_PES %>% filter(Scaled_DBP_PGS_full > 1.280854725)
High_PGS_dat$Scaled_K_new <- as.numeric(scale(High_PGS_dat$K_urine_millimolL))
Low_PGS_dat <- Merged_PGS_PES %>% filter(Scaled_DBP_PGS_full < -1.285330676)
Low_PGS_dat$Scaled_K_new <- as.numeric(scale(Low_PGS_dat$K_urine_millimolL))

DBP_K_high_PGS <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                       Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                       PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                       Scaled_K_new, data = High_PGS_dat)

DBP_K_low_PGS <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                      Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                      PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                      Scaled_K_new, data = Low_PGS_dat)

Low_PGS_dat$fitted <- predict(DBP_K_low_PGS)
High_PGS_dat$fitted <- predict(DBP_K_high_PGS)


K_plot_high_PGS <- ggplot(data = High_PGS_dat, aes(x = Scaled_K_new, y = Mean_DBP)) + 
  geom_smooth(method = "lm", se = TRUE, aes(y=fitted)) +
  ylab("Mean DBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Potassium (SD)") +
  theme(legend.position = "none") +
  ggtitle("> 90th percentile genome-wide PGS")

K_plot_low_PGS <- ggplot(data = Low_PGS_dat, aes(x = Scaled_K_new, y = Mean_DBP)) + 
  geom_smooth(method = "lm", se = TRUE,colour = "red", aes(y=fitted)) +
  ylab("Mean DBP (mmHg)") +
  theme_bw() +
  xlab("Spot Scaled Urinary Potassium (SD)") +
  theme(legend.position = "none") +
  ggtitle("< 10th percentile genome-wide PGS")

ggarrange(K_plot_high_PGS, K_plot_low_PGS, nrow = 2)
