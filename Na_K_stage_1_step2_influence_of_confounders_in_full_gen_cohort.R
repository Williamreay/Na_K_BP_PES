#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 1: Step 2 - test the impact of plausible confounders on blood pressure and urinary Na/K, also test impact of trait transformation

## William Reay (2022)

#################################################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)

## Read in data generated from previous script

Subcohort_1_dat <- fread("Subcohort_1_Biochem_BP_genetics_EUR.txt")

## Retrieve some additional data: i) Assessment centre and year of birth, i) Month attening assesement centre

Demographics <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Basic_demographics_UKBB_July_2021_Freeze.txt")

Demographics <- Demographics %>% select(eid, `31-0.0`, `34-0.0`, `54-0.0`,
                                        `55-0.0`, `52-0.0`)
Demographics$YofB <- Demographics$`34-0.0`
Demographics$Assesment_centre_month <-  as.factor(Demographics$`55-0.0`)
Demographics$Assesment_centre <- Demographics$`54-0.0`

Demographics <- Demographics %>% select(eid, YofB, Assesment_centre_month, Assesment_centre)

Merged_subcohort_1 <- merge(Subcohort_1_dat, Demographics, by = "eid")

## Flag and remove Na+ results if < 10 millimol/L || > 400 millimol/L && Flag K= results if < 2 millimol/L || > 200 millimol/L

Filt_Na_K_subcohort_1 <- Merged_subcohort_1 %>%
  filter(Na_urine_millimolL > 10 & Na_urine_millimolL < 400) %>% filter(K_urine_millimolL > 2 & K_urine_millimolL < 200)

table1( ~  factor(`BP Medication`) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL | factor(Sex),
        data = Filt_Na_K_subcohort_1)

## 52 individuals removed

## Define function to test variance explained by each of the following on SBP, DBP, Urine Na, Urine K (raw then log trans):
## Age (Year of Birth)
## Age at recruitment^2
## Sex
## Assessment Centre
## 20 Genetic PCs as markers of ancestry
## Month of assessment.
## Date of urine test
## Blood pressure medication – also perform analyses in unmedicated subjects only.
## BMI

Multiple_r2_func <- function(Outcome, Predictor, df) {
  fmla <- as.formula(paste(Outcome, "~", Predictor))
  mod <- lm(fmla, data = df)
  r2 <- summary(mod)$adj.r.squared
  return(r2)
}

List_of_predictors <- as.list(colnames(Filt_Na_K_subcohort_1[,c(5:27,29,32,43:44)]))

SBP <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Mean_SBP", df = Filt_Na_K_subcohort_1)
SBP <- as.data.frame(SBP)
DBP <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Mean_DBP", df = Filt_Na_K_subcohort_1)
DBP <- as.data.frame(DBP)
Na_urine <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Na_urine_millimolL", df = Filt_Na_K_subcohort_1)
Na_urine <- as.data.frame(Na_urine)
K_urine <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "K_urine_millimolL", df = Filt_Na_K_subcohort_1)
K_urine <- as.data.frame(K_urine)

Merged_regression_out <- bind_cols(SBP, DBP, Na_urine, K_urine)
Merged_regression_out$Phenotype <- unlist(List_of_predictors)

write.table(Merged_regression_out, file="Stage1_cross_sectional_results/Multiple_r2_univariate_predictors_raw_values_BP_and_urinary.txt",
            sep = "\t",row.names = F, quote = F)

Melted_raw <- melt(Merged_regression_out, id.vars=c("Phenotype"), 
                   measure.vars = c("SBP", "DBP", "Na_urine", "K_urine"))

## Plot Multiple R2

RAW<- ggplot(data = Melted_raw, aes(x=Phenotype, y=value, fill = as.factor(Phenotype))) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  facet_wrap(~variable,strip.position="top",nrow=3,scales = "free_y") +
  coord_flip() + theme_bw() +
  theme(legend.position = "none") +
  ylim(-0.001,0.15) +
  ylab(expression(R^{"2"})) +
  ggtitle("Untransformed outcome")

## Repeat above with log transformed outcomes

Filt_Na_K_subcohort_1$Log_SBP <- log(Filt_Na_K_subcohort_1$Mean_SBP)
Filt_Na_K_subcohort_1$Log_DBP <- log(Filt_Na_K_subcohort_1$Mean_DBP)
Filt_Na_K_subcohort_1$Log_Na <- log(Filt_Na_K_subcohort_1$Na_urine_millimolL)
Filt_Na_K_subcohort_1$Log_K <- log(Filt_Na_K_subcohort_1$K_urine_millimolL)

SBP <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Log_SBP", df = Filt_Na_K_subcohort_1)
SBP <- as.data.frame(SBP)
DBP <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Log_DBP", df = Filt_Na_K_subcohort_1)
DBP <- as.data.frame(DBP)
Na_urine <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Log_Na", df = Filt_Na_K_subcohort_1)
Na_urine <- as.data.frame(Na_urine)
K_urine <- sapply(List_of_predictors, Multiple_r2_func, Outcome = "Log_K", df = Filt_Na_K_subcohort_1)
K_urine <- as.data.frame(K_urine)

Merged_regression_out_LOG <- bind_cols(SBP, DBP, Na_urine, K_urine)
Merged_regression_out_LOG$Phenotype <- unlist(List_of_predictors)

write.table(Merged_regression_out_LOG, file="Stage1_cross_sectional_results/Multiple_r2_univariate_predictors_log_trans_values_BP_and_urinary.txt",
            sep = "\t",row.names = F, quote = F)

Melted_raw_LOG <- melt(Merged_regression_out_LOG, id.vars=c("Phenotype"), 
                   measure.vars = c("SBP", "DBP", "Na_urine", "K_urine"))

## Plot Multiple R2

LOG <- ggplot(data = Melted_raw_LOG, aes(x=Phenotype, y=value, fill = as.factor(Phenotype))) +
  geom_bar(stat="identity", position=position_dodge(), colour = "black") +
  facet_wrap(~variable,strip.position="top",nrow=3,scales = "free_y") +
  coord_flip() + theme_bw() +
  theme(legend.position = "none") +
  ylim(-0.001,0.15) +
  ylab(expression(R^{"2"})) +
  ggtitle("Natural log transformed outcome")

## Derive Na/K ratio and output

Filt_Na_K_subcohort_1$Na_K_ratio <- Filt_Na_K_subcohort_1$Na_urine_millimolL/Filt_Na_K_subcohort_1$K_urine_millimolL

write.table(Filt_Na_K_subcohort_1, file = "Subcohort_1_FILTERED_CORRECTLY_FOR_Na_K_RATIO_INCL.txt",
            sep = "\t", row.names = F, quote = F)
