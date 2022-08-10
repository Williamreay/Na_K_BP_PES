#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 3 - correlation amongst PES and PRS and mClust

## Medicated and unmedicated separately 

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
library(mclust)
library(corrplot)
library(factoextra)

set.seed(3838324)

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

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")

## Corrplot in entire cohort

Cor_all <- cor(Merged_PGS_PES[,c(2:9)])

pval_all <- psych::corr.test(Merged_PGS_PES[,c(2:9)], adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_all, tl.cex = 0.5, p.mat = pval_all, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

## Make unmedicated cohort and scale

Unmed <- Merged_PGS_PES %>% filter(BP_meds == 0)

Unmed[, c(2:9)] <-  lapply(Unmed[, c(2:9)], function(x) c(scale(x)))

Unmed_GMM <- Unmed[, c(2:9)]

## Corrplot
Cor_Unmed <- cor(Unmed_GMM)

pval_unmed <- psych::corr.test(Unmed_GMM, adjust = "none")$p

Sig_corr_plot <- corrplot(Cor_Unmed, tl.cex = 0.5, p.mat = pval_unmed, insig = "blank",
                          sig.level = 0.05, type="upper", method="number",
                          tl.col="black", tl.srt=45)

Unmed_BIC <- mclustBIC(Unmed_GMM)
plot(Unmed_BIC)

## SBP only

Unmed_GMM_SBP <- Unmed_GMM %>% select(contains("SBP"))
Unmed_BIC_SBP <- mclustBIC(Unmed_GMM_SBP)
plot(Unmed_BIC_SBP)
SBP_mod <- Mclust(Unmed_GMM_SBP, x =  Unmed_BIC_SBP)
summary(SBP_mod, parameters = T)

GMM_classification_SBP <- as.data.frame(SBP_mod$classification)
SBP_p <- bind_cols(Unmed_GMM_SBP, GMM_classification_SBP)
SBP_p_2 <- bind_cols(SBP_p, Unmed)

ggplot(data = SBP_p_2, aes(x=`SBP_mod$classification`, y=Mean_SBP, group = `SBP_mod$classification`,
       fill = as.factor(`SBP_mod$classification`))) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 136.4025, lty = 2, col = "red")

TE_SBP <- lm(Mean_SBP ~ as.factor(`SBP_mod$classification`), data = SBP_p_2)

## DBP only - one cluster only

Unmed_GMM_DBP <- Unmed_GMM %>% select(contains("DBP"))
Unmed_BIC_DBP <- mclustBIC(Unmed_GMM_DBP)
plot(Unmed_BIC_DBP)
DBP_mod <- Mclust(Unmed_GMM_DBP, x =  Unmed_BIC_DBP)
summary(DBP_mod, parameters = T)
