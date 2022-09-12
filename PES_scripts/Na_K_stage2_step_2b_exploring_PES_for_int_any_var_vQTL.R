#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2c - investigating if any of the SNPs in the transport PES are vQTLs for SBP

## William Reay (2022)

#################################################

library(dplyr)
library(data.table)
library(ggplot2)

setwd("~/Documents/Na_K_BP_Erin_Clare/")

PES_SNPs_unclumped <- fread("Na_K_PES_construction/Pathway_specific_variants/SBP_var/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP.txt")
PES_SNPs_unclumped  <- PES_SNPs_unclumped[,c(1:10)]
PES_SNPs_unclumped  <- rename(PES_SNPs_unclumped, "SNP"="dbSNP")

PES_SNPs_unclumped_thres <- PES_SNPs_unclumped %>% filter(P < 0.005)

## Import vQTL data (SBP - unmedicated)

SBP_vQTL <- fread("~/Desktop/Blood_pressure_target_ranking/Blood_pressure_UKBB_vQTL/OSCA_vQTL_output/SBP_meds_excl/SBP_meds_excl_all_CHR_Levene_median.txt")

## Investigate Z score dist in each

PES_vQTL <- merge(PES_SNPs_unclumped_thres, SBP_vQTL, by = "SNP")

PES_vQTL %>% filter(P.y < 0.01)

## Not much there 




