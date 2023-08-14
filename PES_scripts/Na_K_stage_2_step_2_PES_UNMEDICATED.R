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
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + 
                                     PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}

## Run list of scores

Run_PES_BP <- mapply(PES_BP_assoc, score = List_PES, BP =  List_traits,
                     MoreArgs = list(Merged_PES),
                     SIMPLIFY = T)

## Extract relevant data

Run_BP_extract <- apply(Run_PES_BP, 2, function(x) return(as.data.frame(x$coefficients)[57, 1:4]))

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
                                     Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + 
                                     PC6 + PC7 + PC8 + PC9 + PC10 + PC11 +
                                     PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + ", PGS, "+ ", score))
  PES_mod <- lm(fmla, data = df)
  return(summary(PES_mod))
}

## Run list of scores

Run_PES_BP_PGS_adjust <- mapply(PES_BP_assoc_PGS_adjust, score = List_PES, BP =  List_traits, PGS = List_PGS,
                                MoreArgs = list(Merged_PGS_PES),
                                SIMPLIFY = T)

## Extract relevant data

Run_BP_extract_PGS_adjust <- apply(Run_PES_BP_PGS_adjust, 2, function(x) return(as.data.frame(x$coefficients)[58, 1:4]))

Run_BP_results_PGS_adjust <- data.frame()

for (i in 1:length(Run_BP_extract_PGS_adjust)) {
  Run_BP_results_PGS_adjust <- rbind(Run_BP_results_PGS_adjust, Run_BP_extract_PGS_adjust[[i]])
}

Run_BP_results_PGS_adjust$PES <- unlist(List_PES)
Run_BP_results_PGS_adjust$BP_trait <- unlist(List_traits)
Run_BP_results_PGS_adjust$PGS <- "PGS adjusted"

PES_BP_results_unmed <- rbind(Run_BP_results, Run_BP_results_PGS_adjust)

write.table(PES_BP_results_unmed, file = "UKBB_PES/PES_BP_assoc_and_int/PES_assoc_FULL_cohort_with_and_without_PGS_adjustment.txt",
            sep = "\t", row.names = F, quote = F)

PES_BP_results_unmed$L_UI <- PES_BP_results_unmed$Estimate - (1.96*PES_BP_results_unmed$`Std. Error`)
PES_BP_results_unmed$U_UI <- PES_BP_results_unmed$Estimate + (1.96*PES_BP_results_unmed$`Std. Error`)

PES_BP_results_unmed$PES_name <- c("Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport",
                                         "Absorption", "Renal", "Transport")

SBP_results <- PES_BP_results_unmed %>% filter(BP_trait == "Mean_SBP")
DBP_results <- PES_BP_results_unmed %>% filter(BP_trait == "Mean_DBP")

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
  ylim(-0.05, 0.6) +
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
  ylim(-0.05, 0.6) +
  ggtitle("DBP")

ggarrange(SBP_full, DBP_full)


