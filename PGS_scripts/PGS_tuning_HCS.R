##########################

## Tuning the best performing BP PGS in the HCS

## William Reay (2022)

########################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(rcompanion)
library(ggplot2)
library(viridis)
library(ggpubr)

## PERFORM ANALYSES IN MEDICATED AND UNMEDICATED SEPARTELY 

## Define regression automation function

PGS_tuning_r2_func <- function(Outcome, Score, df) {
  Input_df <- df
  Input_df$scaled_score <- as.numeric(scale(Input_df[[Score]]))
  baseline_fmla <- as.formula(paste(Outcome, "~ + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5"))
  PGS_fmla <- as.formula(paste(Outcome, "~ + bage + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + scaled_score"))
  mod1 <- lm(baseline_fmla, data = Input_df)
  mod2 <- lm(PGS_fmla, data=Input_df)
  mod1_r2 <- summary(mod1)$adj.r.squared
  mod2_r2 <- summary(mod2)$adj.r.squared
  Delta_r2 <- mod2_r2 - mod1_r2
  return(Delta_r2)
}

## Read in data

HCS_pheno <- fread("HCS_cohort_tuning_of_PRS_PES/HCS_BP_cohort.txt")

HCS_pheno <- rename(HCS_pheno, "IID"="electoralId")

DBP_PGS <- fread("Genetic_data/BP_HCS_scores/DBP_PGS_HCS.txt")
Colnames_DBP_PGS <- make.unique(names(DBP_PGS))
colnames(DBP_PGS) <- Colnames_DBP_PGS

SBP_PGS <- fread("Genetic_data/BP_HCS_scores/SBP_PGS_HCS.txt")
Colnames_SBP_PGS <- make.unique(names(SBP_PGS))
colnames(SBP_PGS) <- Colnames_SBP_PGS

Merged_cohort <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                        list(DBP_PGS, SBP_PGS, HCS_pheno))

## Split into medicated and unmedicated cohorts

Med_HCS <- Merged_cohort %>% filter(BP_med == "YES")
Unmed_HCS <- Merged_cohort %>% filter(BP_med == "NO")

## Define list of PGS to test

SBP_to_test <- as.list(colnames(Merged_cohort %>% select(starts_with("SBP"))))
DBP_to_test <- as.list(colnames(Merged_cohort %>% select(starts_with("DBP"))))

## SBP
SBP_score_med <- sapply(SBP_to_test, PGS_tuning_r2_func, df = Med_HCS, Outcome = "mbs")
SBP_score_med <- as.data.frame(SBP_score_med)
SBP_score_unmed <- sapply(SBP_to_test, PGS_tuning_r2_func, df = Unmed_HCS, Outcome = "mbs")
SBP_score_unmed <- as.data.frame(SBP_score_unmed)

SBP_out <- bind_cols(SBP_score_med, SBP_score_unmed)
SBP_out$PGS <- unlist(SBP_to_test)

write.table(SBP_out, file="HCS_cohort_tuning_of_PRS_PES/SBP_HCS_tuning_results.txt",
            sep = "\t", row.names = F, quote = F)

## DBP
DBP_score_med <- sapply(DBP_to_test, PGS_tuning_r2_func, df = Med_HCS, Outcome = "mbd")
DBP_score_med <- as.data.frame(DBP_score_med)
DBP_score_unmed <- sapply(DBP_to_test, PGS_tuning_r2_func, df = Unmed_HCS, Outcome = "mbd")
DBP_score_unmed <- as.data.frame(DBP_score_unmed)

DBP_out <- bind_cols(DBP_score_med, DBP_score_unmed)
DBP_out$PGS <- unlist(DBP_to_test)

write.table(DBP_out, file="HCS_cohort_tuning_of_PRS_PES/DBP_HCS_tuning_results.txt",
            sep = "\t", row.names = F, quote = F)

## Plot R2 per score in the medicated and unmedicated cohort

Plot_SBP_melt <- melt(SBP_out, id.vars=c("PGS"), 
     measure.vars = c("SBP_score_med", "SBP_score_unmed"))

Plot_SBP_melt$P <- c("0.001", "0.005", "0.05", "0.01", "0.1", "0.5", "1e-04", "1e-05", "1", "5e-08",
                     "0.001", "0.005", "0.05", "0.01", "0.1", "0.5", "1e-04", "1e-05", "1", "5e-08")

Plot_SBP_melt$P <- factor(Plot_SBP_melt$P, levels = c("5e-08",
                                                      "1e-05", 
                                                      "1e-04", 
                                                      "0.001",
                                                      "0.005",
                                                      "0.01", 
                                                      "0.05", 
                                                      "0.1", 
                                                      "0.5", 
                                                      "1"))

Names <-c("SBP_score_med"="Baseline Antihypertensives (N=981)", "SBP_score_unmed"="Unmedicated (N=852)")

R2_SBP <- ggplot(Plot_SBP_melt, aes(x = P, y = value, colour = as.factor(PGS), group = 1)) +
  facet_wrap(~variable,strip.position="top",nrow=2, scales = "free_y", labeller = as_labeller(Names)) +
  theme_bw() + geom_point(cex = 3) +
  geom_line(colour = "black", linetype = "dashed") +
  theme(legend.position = "none") +
  ylim(-0.001, 0.02) +
  ylab(expression(Delta~R^{"2"})) +
  xlab("P-value threshold") +
  scale_colour_viridis(discrete = T, option="magma") +
  ggtitle("SBP")

Plot_DBP_melt <- melt(DBP_out, id.vars=c("PGS"), 
                      measure.vars = c("DBP_score_med", "DBP_score_unmed"))

Plot_DBP_melt$P <- c("0.001", "0.005", "0.05", "0.01", "0.1", "0.5", "1e-04", "1e-05", "1", "5e-08",
                     "0.001", "0.005", "0.05", "0.01", "0.1", "0.5", "1e-04", "1e-05", "1", "5e-08")

Plot_DBP_melt$P <- factor(Plot_DBP_melt$P, levels = c("5e-08",
                                                      "1e-05", 
                                                      "1e-04", 
                                                      "0.001",
                                                      "0.005",
                                                      "0.01", 
                                                      "0.05", 
                                                      "0.1", 
                                                      "0.5", 
                                                      "1"))

D_Names <-c("DBP_score_med"="Baseline Antihypertensives (N=981)", "DBP_score_unmed"="Unmedicated (N=852)")

R2_DBP <- ggplot(Plot_DBP_melt, aes(x = P, y = value, colour = as.factor(PGS), group = 1)) +
  facet_wrap(~variable,strip.position="top",nrow=2, scales = "free_y", labeller = as_labeller(D_Names)) +
  theme_bw() + geom_point(cex = 3) +
  geom_line(colour = "black", linetype = "dashed") +
  theme(legend.position = "none") +
  ylim(-0.001, 0.02) +
  ylab(expression(Delta~R^{"2"})) +
  xlab("P-value threshold") +
  scale_colour_viridis(discrete = T, option="magma") +
  ggtitle("DBP")


ggarrange(R2_SBP, R2_DBP, nrow = 1, ncol = 2)


