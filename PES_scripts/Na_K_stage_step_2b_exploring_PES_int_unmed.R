#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2b - exploring transport PES x sodium interaction - UNMEDICATED ONLY

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

Na_K_PES <- fread("UKBB_PES/PES_Na_K_SBP_DBP.txt")

Colnames_PES <- make.unique(names(Na_K_PES))

colnames(Na_K_PES) <- Colnames_PES

Na_K_PES <- Na_K_PES %>% select(IID, contains("PES"))

Merged_PES <- merge(Na_K_PES, Subcohort_1_dat, by = "IID")

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")
Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_AVG))

## Test effect of scaled Na:K ratio
Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_K_ratio))

TE <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
  Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
  PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Transport_SBP_PES_AVG*Scaled_Na, data = Merged_PGS_PES)

Ratio_Transport_int_G_C_C_E <- lm(Mean_SBP ~ Sex*Transport_SBP_PES_AVG + Age*Transport_SBP_PES_AVG + Age2*Transport_SBP_PES_AVG + Assesment_centre_month*Transport_SBP_PES_AVG + 
                              Assesment_centre*Transport_SBP_PES_AVG + PC1*Transport_SBP_PES_AVG + PC2*Transport_SBP_PES_AVG+ 
                              PC3*Transport_SBP_PES_AVG + PC4*Transport_SBP_PES_AVG + PC5*Transport_SBP_PES_AVG + PC10*Transport_SBP_PES_AVG+ PC11*Transport_SBP_PES_AVG +
                              PC12*Transport_SBP_PES_AVG + PC13*Transport_SBP_PES_AVG + PC14*Transport_SBP_PES_AVG + PC15*Transport_SBP_PES_AVG + 
                              PC16*Transport_SBP_PES_AVG + PC17*Transport_SBP_PES_AVG + PC18*Transport_SBP_PES_AVG + PC19*Transport_SBP_PES_AVG + PC20*Transport_SBP_PES_AVG + 
                              Scaled_Na*Transport_SBP_PES_AVG +
                              Sex*Scaled_Na + Age*Scaled_Na + Age2*Scaled_Na + Assesment_centre_month*Scaled_Na + 
                              Assesment_centre*Scaled_Na + PC1*Scaled_Na + PC2*Scaled_Na+ PC3*Scaled_Na + PC4*Scaled_Na + PC5*Scaled_Na + PC10*Scaled_Na+ PC11*Scaled_Na +
                              PC12*Scaled_Na + PC13*Scaled_Na + PC14*Scaled_Na + PC15*Scaled_Na + PC16*Scaled_Na + PC17*Scaled_Na + PC18*Scaled_Na + PC19*Scaled_Na + PC20*Scaled_Na,
                            data = Merged_PGS_PES)

summary(Ratio_Transport_int_G_C_C_E)

x <- Merged_PGS_PES$Transport_SBP_PES_AVG

Merged_PGS_PES$PES <-
  case_when(x > mean(x)+sd(x) ~ "High",
            x < mean(x)+sd(x) & x > mean(x)-sd(x) ~ "Average",
            x < mean(x)-sd(x) ~ "Low")

Merged_PGS_PES %>% 
  ggplot() +
  aes(x = Scaled_Na, y = Mean_SBP, group = PES, color = PES) +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Spot Na:K ratio (SD)") +
  ylab("SBP (mmHg)") +
  ggtitle("Main additive effects")

Merged_PGS_PES$Scaled_K <- as.numeric(scale(Merged_PGS_PES$K_urine_millimolL))

Merged_PGS_PES %>% 
  ggplot() +
  aes(x = Scaled_K, y = Mean_SBP, group = PES, color = PES) +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("Spot Urinary K (SD)") +
  ylab("SBP (mmHg)") +
  ggtitle("Main additive effects")

## Calculate the effect of PES and PRS per decile of sodium

Merged_PGS_PES$Decile_Na <- ntile(Merged_PGS_PES$Na_urine_millimolL, 10)


Na_assoc_per_decile <- function(score, decile, Na_col, df) {
  new_df <- df[df[[Na_col]] == decile,]
  new_df$scaled_score <- as.numeric(scale(new_df[[score]]))
  mod <- lm(Mean_SBP ~ Sex + Age + Age2 + Assesment_centre_month + 
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + scaled_score, data = new_df)
  return(summary(mod))
}

Vec_decile <- c(1,2,3,4,5, 6, 7, 8, 9, 10)

PES <- sapply(Vec_decile, Na_assoc_per_decile, Na_col = "Decile_Na",
              df = Merged_PGS_PES, score = "Transport_SBP_PES_AVG")
PRS <- sapply(Vec_decile, Na_assoc_per_decile, Na_col = "Decile_Na",
              df = Merged_PGS_PES, score = "SBP_PGS_Pt_0_01_AVG")


PES_extract <- apply(PES, 2, function(x) return(as.data.frame(x$coefficients)[53, 1:4]))
PRS_extract <- apply(PRS, 2, function(x) return(as.data.frame(x$coefficients)[53, 1:4]))

PES_plot <- data.frame()

for (i in 1:length(PES_extract)) {
  PES_plot <- rbind(PES_plot, PES_extract[[i]])
}

PES_plot$Decile <- c("1st-10th percentile", "10th-20th percentile", "20th-30th percentile", "30th-40th percentile", "40th-50th percentile",
                       "50th-60th percentile", "60th-70th percentile", "70th-80th percentile", "80th-90th percentile", "Top decile")

PRS_plot <- data.frame()

for (i in 1:length(PRS_extract)) {
  PRS_plot <- rbind(PRS_plot, PRS_extract[[i]])
}

PRS_plot$Decile <- c("1st-10th percentile", "10th-20th percentile", "20th-30th percentile", "30th-40th percentile", "40th-50th percentile",
                       "50th-60th percentile", "60th-70th percentile", "70th-80th percentile", "80th-90th percentile", "Top decile")

## Test stat sig effect of top decile vs the rest

Stat_compare_PES <- function(decile) {
  Z <- (PES_plot$Estimate[[10]]-PES_plot$Estimate[[decile]])/(sqrt((PES_plot$`Std. Error`[[10]]^2)+(PES_plot$`Std. Error`[[decile]]^2)))
  P<- pnorm(abs(Z), lower.tail = F)*2
  return(P)
}

Stat_compare_PRS <- function(decile) {
  Z <- (PRS_plot$Estimate[[10]]-PRS_plot$Estimate[[decile]])/(sqrt((PRS_plot$`Std. Error`[[10]]^2)+(PRS_plot$`Std. Error`[[decile]]^2)))
  P<- pnorm(abs(Z), lower.tail = F)*2
  return(P)
}


Val <- c(1,2,3,4,5,6,7,8,9,10)
sapply(Val, Stat_compare_PES)
sapply(Val, Stat_compare_PRS)


## PES Na est per quintile

PES_plot$L_CI <- PES_plot$Estimate - (1.96*PES_plot$`Std. Error`)
PES_plot$U_CI <- PES_plot$Estimate + (1.96*PES_plot$`Std. Error`)
PES_plot$Val <- c(1,2,3,4,5,6,7,8,9,10)
mean(PES_plot$Estimate)
mean(PES_plot$Estimate)-sd(PES_plot$Estimate)
mean(PES_plot$Estimate)+sd(PES_plot$Estimate)

P_P <- ggplot(PES_plot, aes(Val, Estimate)) +
  geom_errorbar(
    aes(ymin = L_CI, ymax = U_CI, color = as.factor(Val))) +
  geom_point(aes(color = as.factor(Val)), position = position_dodge(0.3)) +
  stat_smooth(method = "loess", se = F, size = 0.5, colour = "black") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in score") +
  xlab("Urinary Sodium (deciles)") +
  geom_hline(yintercept = 0.47, lty = "dashed") +
  ggtitle("Transport PES") +
  ylim(0, 1.95) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.6038735, ymax = Inf, fill = "grey", alpha = .2, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.3420018, fill = "grey", alpha = .2, color = NA)

PES_final <- P_P + scale_x_continuous(breaks = 1:10,
                         labels = paste(c("1st", "2nd", "3rd", "4th", 
                                           "5th", "6th", "7th", "8th", "9th", "10th"), " Decile")) +
  theme(axis.text.x = element_text(angle = 90))

## PRS
PRS_plot$Val <- c(1,2,3,4,5,6,7,8,9,10)
mean(PRS_plot$Estimate)
mean(PRS_plot$Estimate)-sd(PRS_plot$Estimate)
mean(PRS_plot$Estimate)+sd(PRS_plot$Estimate)

PRS_plot$L_CI <- PRS_plot$Estimate - (1.96*PRS_plot$`Std. Error`)
PRS_plot$U_CI <- PRS_plot$Estimate + (1.96*PRS_plot$`Std. Error`)

P_PRS <- ggplot(PRS_plot, aes(Val, Estimate)) +
  geom_errorbar(
    aes(ymin = L_CI, ymax = U_CI, color = as.factor(Val))) +
  geom_point(aes(color = as.factor(Val)), position = position_dodge(0.3)) +
  stat_smooth(method = "loess", se = F, size = 0.5, colour = "black") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab(" ") +
  xlab("Urinary Sodium (deciles)") +
  geom_hline(yintercept = 1.606948, lty = "dashed") +
  ggtitle("Genome-wide PGS") +
  ylim(0, 1.95) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.697952, ymax = Inf, fill = "grey", alpha = .2, color = NA) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1.515944, fill = "grey", alpha = .2, color = NA)

PRS_final <- P_PRS + scale_x_continuous(breaks = 1:10,
                                      labels = paste(c("1st", "2nd", "3rd", "4th", 
                                                       "5th", "6th", "7th", "8th", "9th", "10th"), " Decile")) +
  theme(axis.text.x = element_text(angle = 90))

ggarrange(PES_final, PRS_final, nrow = 1)

## Test equality of BP variance between deciles of PES and PGS, analogous to vQTL approach

## Need to preprocess data

Pre_process <- function(BP_col, BP_name, sex, sex_col, df) {
  test_df <-   df[!is.na(df[[BP_col]]),]
  test_df <- test_df[test_df[[sex_col]] == sex,]
  fmla <- as.formula(paste(BP_col," ~ Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 +
                           PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 +
                           PC18 + PC19 + PC19 + PC20 + Assesment_centre + Assesment_centre_month"))
  mod <- lm(fmla, data = test_df)
  ## Scaled residuals
  Residuals_RAW <- as.data.frame(mod$residuals)
  colnames(Residuals_RAW)[1] <- "Resid"
  test_df$residuals <- Residuals_RAW$Resid
  BP_mean <- mean(test_df$residuals)
  BP_SD <- sd(test_df$residuals)
  Five_BP_SD_low = BP_mean - (BP_SD * 5)
  Five_BP_SD_high = BP_mean + (BP_SD * 5)
  test_df <- test_df[residuals > Five_BP_SD_low]
  test_df <- test_df[residuals < Five_BP_SD_high]
  test_df$Norm_resid <- as.numeric(scale(test_df$residuals))
  names(test_df)[names(test_df) == "Norm_resid"] <- BP_name
  return(test_df)
}

Females_BP <- Pre_process("Mean_SBP", "Scaled_SBP", "Female", "Sex", Merged_PGS_PES)
hist(Females_BP$Scaled_SBP)
Males_BP <- Pre_process("Mean_SBP", "Scaled_SBP", "Male", "Sex", Merged_PGS_PES)
hist(Males_BP$Scaled_SBP)

Variance_df <- rbind(Females_BP, Males_BP)

Variance_df$PES_decile <- as.factor(ntile(Variance_df$Transport_SBP_PES_AVG, 10))
Variance_df$PRS_decile <- as.factor(ntile(Variance_df$SBP_PGS_Pt_0_01_AVG, 10))

## Levene's test

leveneTest(Scaled_SBP ~ as.factor(PES_decile), data = Variance_df)
leveneTest(Scaled_SBP ~ as.factor(PRS_decile), data = Variance_df)

## Some nominal evidence for non-homogeneous variance between PES deciles in terms of SBP, moreso in terms of PGS


