##########################

## Tuning the best performing BP Na/K PES in the HCS (C+T)

## William Reay (2022)

########################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(rcompanion)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readxl)

## Get regression automation function from other script
source("~/Documents/Na_K_BP_Erin_Clare/HCS_cohort_tuning_of_PRS_PES/PGS_tuning_HCS.R")

## Read in data

HCS_pheno <- fread("HCS_cohort_tuning_of_PRS_PES/HCS_BP_cohort.txt")

HCS_pheno <- rename(HCS_pheno, "IID"="electoralId")

BP_PES <- fread("Na_K_PES_construction/PES_BP_Na_K/SBP_DBP_Na_K_PES_HCS.txt")
Colnames_BP_PES <- make.unique(names(BP_PES))
colnames(BP_PES) <- Colnames_BP_PES

Merged_cohort <- merge(BP_PES, HCS_pheno, by ="IID")

## Split into medicated and unmedicated cohorts

Med_HCS <- Merged_cohort %>% filter(BP_med == "YES")
Unmed_HCS <- Merged_cohort %>% filter(BP_med == "NO")

## Define list of PGS to test

SBP_to_test <- as.list(colnames(Merged_cohort %>% select(contains("SBP"))))
DBP_to_test <- as.list(colnames(Merged_cohort %>% select(contains("DBP"))))

## SBP
SBP_score_med <- sapply(SBP_to_test, PGS_tuning_r2_func, df = Med_HCS, Outcome = "mbs")
SBP_score_med <- as.data.frame(SBP_score_med)
SBP_score_unmed <- sapply(SBP_to_test, PGS_tuning_r2_func, df = Unmed_HCS, Outcome = "mbs")
SBP_score_unmed <- as.data.frame(SBP_score_unmed)

SBP_out <- bind_cols(SBP_score_med, SBP_score_unmed)
SBP_out$PGS <- unlist(SBP_to_test)

write.table(SBP_out, file="Na_K_PES_construction/HCS_tuning_results_PES/SBP_C_T_tuning_PES.txt",
            sep = "\t", row.names = F, quote = F)

DBP_out$Av_r2 <- (DBP_out$DBP_score_med+DBP_out$DBP_score_unmed)/2
SBP_out$Av_r2 <- (SBP_out$SBP_score_med+SBP_out$SBP_score_unmed)/2

## DBP
DBP_score_med <- sapply(DBP_to_test, PGS_tuning_r2_func, df = Med_HCS, Outcome = "mbd")
DBP_score_med <- as.data.frame(DBP_score_med)
DBP_score_unmed <- sapply(DBP_to_test, PGS_tuning_r2_func, df = Unmed_HCS, Outcome = "mbd")
DBP_score_unmed <- as.data.frame(DBP_score_unmed)

DBP_out <- bind_cols(DBP_score_med, DBP_score_unmed)
DBP_out$PGS <- unlist(DBP_to_test)

write.table(DBP_out, file="Na_K_PES_construction/HCS_tuning_results_PES/DBP_C_T_tuning_PES.txt",
            sep = "\t", row.names = F, quote = F)

## Identify the largest R2 and check association

max(DBP_out$DBP_score_med) ## RENAL_sodium_pathways_genes_sumstats_DBP_0.005_DBP_AVG
max(DBP_out$DBP_score_unmed) ## RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG
max(SBP_out$SBP_score_med) ## Aborption_pathways_genes_sumstats_SBP_0.05_SBP_AVG
max(SBP_out$SBP_score_unmed) ## TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG



## Define function to test for statistical significance

PGS_tuning_stat_sig <- function(Outcome, Score, df) {
  Input_df <- df
  Input_df$scaled_score <- as.numeric(scale(Input_df[[Score]]))
  PGS_fmla <- as.formula(paste(Outcome, "~ + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + scaled_score"))
  mod1 <- lm(PGS_fmla, data=Input_df)
  return(summary(mod1))
}

## SBP

Stat_SBP_score_med <- sapply(SBP_to_test, PGS_tuning_stat_sig, df = Med_HCS, Outcome = "mbs")
Stat_SBP_score_med_extract <- apply(Stat_SBP_score_med, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

## No stat sig for medicated

Stat_SBP_score_unmed <- sapply(SBP_to_test, PGS_tuning_stat_sig, df = Unmed_HCS, Outcome = "mbs")
Stat_SBP_score_unmed_extract <- apply(Stat_SBP_score_unmed, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

## Some stat sig for unmedicated

## DBP

Stat_DBP_score_med <- sapply(DBP_to_test, PGS_tuning_stat_sig, df = Med_HCS, Outcome = "mbd")
Stat_DBP_score_med_extract <- apply(Stat_DBP_score_med, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))

## No stat sig for medicated

Stat_DBP_score_unmed <- sapply(DBP_to_test, PGS_tuning_stat_sig, df = Unmed_HCS, Outcome = "mbd")
Stat_DBP_score_unmed_extract <- apply(Stat_DBP_score_unmed, 2, function(x) return(as.data.frame(x$coefficients)[9,1:4]))


## Make results output

Med_SBP_results <- data.frame()

for (i in 1:length(Stat_SBP_score_med_extract)) {
  Med_SBP_results <- rbind(Med_SBP_results, Stat_SBP_score_med_extract[[i]])
}

rownames(Med_SBP_results) <- SBP_to_test
Med_SBP_results$Medicated <- "Medicated"

Unmed_SBP_results <- data.frame()

for (i in 1:length(Stat_SBP_score_unmed_extract)) {
  Unmed_SBP_results <- rbind(Unmed_SBP_results, Stat_SBP_score_unmed_extract[[i]])
}

rownames(Unmed_SBP_results) <- SBP_to_test
Unmed_SBP_results$Medicated <- "Unmedicated"

SBP_stat_sig_out <- rbind(Med_SBP_results, Unmed_SBP_results)

write.table(SBP_stat_sig_out, file="Na_K_PES_construction/HCS_tuning_results_PES/SBP_PES_stat_sig_PGS_unadjusted_HCS.txt",
            sep = "\t",row.names = T, quote = F)

Med_DBP_results <- data.frame()

for (i in 1:length(Stat_DBP_score_med_extract)) {
  Med_DBP_results <- rbind(Med_DBP_results, Stat_DBP_score_med_extract[[i]])
}

rownames(Med_DBP_results) <- DBP_to_test
Med_DBP_results$Medicated <- "Medicated"

Unmed_DBP_results <- data.frame()

for (i in 1:length(Stat_DBP_score_unmed_extract)) {
  Unmed_DBP_results <- rbind(Unmed_DBP_results, Stat_DBP_score_unmed_extract[[i]])
}

rownames(Unmed_DBP_results) <- DBP_to_test
Unmed_DBP_results$Medicated <- "Unmedicated"

DBP_stat_sig_out <- rbind(Med_DBP_results, Unmed_DBP_results)

write.table(DBP_stat_sig_out, file="Na_K_PES_construction/HCS_tuning_results_PES/DBP_PES_stat_sig_PGS_unadjusted_HCS.txt",
            sep = "\t",row.names = T, quote = F)

## Define which are stat sig

SBP_nom_sig <- SBP_stat_sig_out %>% filter(`Pr(>|t|)` < 0.05)
DBP_nom_sig <- DBP_stat_sig_out %>% filter(`Pr(>|t|)` < 0.05)

## PES to carry forward as at least nominally significant in either cohort (test in unmedicated)

## RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG1
## TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG1
## TRANSPORT_sodium_pathways_genes_sumstats_SBP_0.005_SBP_AVG1 
## RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG1 

## Read in PGS for covariation

DBP_PGS <- fread("Genetic_data/BP_HCS_scores/DBP_PGS_HCS.txt")
Colnames_DBP_PGS <- make.unique(names(DBP_PGS))
colnames(DBP_PGS) <- Colnames_DBP_PGS

SBP_PGS <- fread("Genetic_data/BP_HCS_scores/SBP_PGS_HCS.txt")
Colnames_SBP_PGS <- make.unique(names(SBP_PGS))
colnames(SBP_PGS) <- Colnames_SBP_PGS

Unmed_Merged_cohort <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                        list(DBP_PGS, SBP_PGS, Unmed_HCS))

Unmed_Merged_cohort$SBP_0_005_AVG <- as.numeric(scale(Unmed_Merged_cohort$SBP_0_005_AVG))
Unmed_Merged_cohort$SBP_0_5_AVG <- as.numeric(scale(Unmed_Merged_cohort$SBP_0_5_AVG))
Unmed_Merged_cohort$DBP_1_AVG <- as.numeric(scale(Unmed_Merged_cohort$DBP_1_AVG))
Unmed_Merged_cohort$RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG <- as.numeric(scale(Unmed_Merged_cohort$RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG))
Unmed_Merged_cohort$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG <- as.numeric(scale(Unmed_Merged_cohort$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG))
Unmed_Merged_cohort$TRANSPORT_sodium_pathways_genes_sumstats_SBP_0.005_SBP_AVG <- as.numeric(scale(Unmed_Merged_cohort$TRANSPORT_sodium_pathways_genes_sumstats_SBP_0.005_SBP_AVG))
Unmed_Merged_cohort$RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG <- as.numeric(scale(Unmed_Merged_cohort$RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG))
## Test each covaried for PGS at same threshold

## RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG1
## TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG1
## TRANSPORT_sodium_pathways_genes_sumstats_SBP_0.005_SBP_AVG1 
## RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG1 

Renal_na_k_SBP_0_5 <- lm(mbs ~  + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + SBP_0_5_AVG +
                           RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG,
                         data =  Unmed_Merged_cohort)

Transport_na_k_SBP_0_005 <- lm(mbs ~  + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + SBP_0_005_AVG +
                               TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG,
                         data =  Unmed_Merged_cohort)

cor.test(Unmed_Merged_cohort$TRANSPORT_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG, Unmed_Merged_cohort$SBP_0_005_AVG)

Transport_na_SBP_0_005 <- lm(mbs ~  + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + SBP_0_005_AVG +
                                 TRANSPORT_sodium_pathways_genes_sumstats_SBP_0.005_SBP_AVG,
                               data =  Unmed_Merged_cohort)

Renal_na_DBP_1 <- lm(mbd ~  + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + DBP_1_AVG +
                           RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG,
                         data =  Unmed_Merged_cohort)
cor.test(Unmed_Merged_cohort$RENAL_potassium_pathways_genes_sumstats_DBP_1_DBP_AVG, Unmed_Merged_cohort$DBP_1_AVG)

PGS_DBP_1 <- lm(mbd ~  + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + DBP_1_AVG,
                     data =  Unmed_Merged_cohort)

Plot_input <- read_excel("Na_K_PES_construction/HCS_tuning_results_PES/PES_plot_input.xlsx")
Plot_input$L_beta <- Plot_input$Beta - (1.96*Plot_input$SE)
Plot_input$U_beta <- Plot_input$Beta + (1.96*Plot_input$SE)

ggplot(Plot_input, aes(x=Score, y=Beta, ymin=L_beta, ymax=U_beta,colour = Trait)) +
  geom_pointrange() +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept=0, lty=2) +
  ylim(-0.1, 3) +
  theme(legend.position = "bottom") +
  ylab("Effect per SD increase in PES on BP (mmHg) adjusted for genome-wide PGS")

R2_plot_input <- read_excel("Na_K_PES_construction/HCS_tuning_results_PES/Best_PES_per_category.xlsx")

ggplot(R2_plot_input, aes(x=Score,  y=Av_r2, fill = Category)) +
  geom_bar(stat="identity",colour = "black", width = 0.75) +
  facet_wrap(~Trait, strip.position="top",nrow=2, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  geom_hline(yintercept=0, lty=2, colour = "red") +
  theme(legend.position = "bottom",plot.title = element_text(size = 10)) +
  ylim(-0.001, 0.01) +
  ylab(expression(Delta~R^{"2"})) +
  ggtitle("Average variance explained of best performing PES per category")

  theme(legend.position = "bottom") +
  ylab("Effect per SD increase in PES on BP (mmHg) adjusted for genome-wide PGS")
  