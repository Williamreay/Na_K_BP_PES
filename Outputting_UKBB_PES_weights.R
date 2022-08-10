#################################################

## Outputting PES weights for the UKBB - full sample

## William Reay (2022)

#################################################

library(dplyr)
library(data.table)

setwd("~/Documents/Na_K_BP_Erin_Clare/")


SBP_sumstats <- fread("Genetic_data/Cleaned_GERA_sumstats_SBP.txt.gz")
DBP_sumstats <- fread("Genetic_data/Cleaned_GERA_sumstats_DBP.txt.gz")

## RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_AVG

RENAL_DBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP_PES_HCS_weights.txt", header = T)
Merged_RENAL_DBP <- merge(DBP_sumstats, RENAL_DBP_weights, by ="SNP")
Merged_RENAL_DBP$Renal_DBP_PES <- Merged_RENAL_DBP$RENAL_sodium_pathways_genes_sumstats_DBP_1_DBP

Merged_RENAL_DBP <- Merged_RENAL_DBP %>% select(dbSNP, EA.y, Renal_DBP_PES)
Merged_RENAL_DBP <- unique(Merged_RENAL_DBP)

write.table(Merged_RENAL_DBP, file="UKBB_PES/Renal_DBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_DBP_1_DBP_AVG

TRANSPORT_DBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_DBP_1_DBP_PES_HCS_weights.txt", header = T)
Merged_TRANSPORT_DBP <- merge(DBP_sumstats, TRANSPORT_DBP_weights, by ="SNP")
Merged_TRANSPORT_DBP$Transport_DBP_PES <- Merged_TRANSPORT_DBP$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_DBP_1_DBP

Merged_TRANSPORT_DBP <- Merged_TRANSPORT_DBP %>% select(dbSNP, EA.y, Transport_DBP_PES)
Merged_TRANSPORT_DBP <- unique(Merged_TRANSPORT_DBP)

write.table(Merged_TRANSPORT_DBP, file="UKBB_PES/Transport_DBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## Aborption_pathways_genes_sumstats_DBP_0.005_DBP_AVG

ABSORPTION_DBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/Aborption_pathways_genes_sumstats_DBP_0.005_DBP_PES_HCS_weights.txt", header = T)
Merged_ABSORPTION_DBP <- merge(DBP_sumstats, ABSORPTION_DBP_weights, by ="SNP")
Merged_ABSORPTION_DBP$Absorption_DBP_PES <- Merged_ABSORPTION_DBP$Aborption_pathways_genes_sumstats_DBP_0.005_DBP

Merged_ABSORPTION_DBP <- Merged_ABSORPTION_DBP %>% select(dbSNP, EA.y, Absorption_DBP_PES)
Merged_ABSORPTION_DBP <- unique(Merged_ABSORPTION_DBP)

write.table(Merged_ABSORPTION_DBP, file="UKBB_PES/Absorption_DBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_AVG

RENAL_SBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_PES_HCS_weights.txt", header = T)
Merged_RENAL_SBP <- merge(SBP_sumstats, RENAL_SBP_weights, by ="SNP")
Merged_RENAL_SBP$Renal_SBP_PES <- Merged_RENAL_SBP$RENAL_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP

Merged_RENAL_SBP <- Merged_RENAL_SBP %>% select(dbSNP, EA.y, Renal_SBP_PES)
Merged_RENAL_SBP <- unique(Merged_RENAL_SBP)

write.table(Merged_RENAL_SBP, file="UKBB_PES/Renal_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_AVG

TRANSPORT_SBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP_PES_HCS_weights.txt", header = T)
Merged_TRANSPORT_SBP <- merge(SBP_sumstats, TRANSPORT_SBP_weights, by ="SNP")
Merged_TRANSPORT_SBP$Transport_SBP_PES <- Merged_TRANSPORT_SBP$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.005_SBP

Merged_TRANSPORT_SBP <- Merged_TRANSPORT_SBP %>% select(dbSNP, EA.y, Transport_SBP_PES)
Merged_TRANSPORT_SBP <- unique(Merged_TRANSPORT_SBP)

write.table(Merged_TRANSPORT_SBP, file="UKBB_PES/Transport_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## Aborption_pathways_genes_sumstats_SBP_0.05_SBP_AVG

ABSORPTION_SBP_weights <- fread("Na_K_PES_construction/PES_weights_CHR_BP/Aborption_pathways_genes_sumstats_SBP_0.05_SBP_PES_HCS_weights.txt", header = T)
Merged_ABSORPTION_SBP <- merge(SBP_sumstats, ABSORPTION_SBP_weights, by ="SNP")
Merged_ABSORPTION_SBP$Absorption_SBP_PES <- Merged_ABSORPTION_SBP$Aborption_pathways_genes_sumstats_SBP_0.05_SBP

Merged_ABSORPTION_SBP <- Merged_ABSORPTION_SBP %>% select(dbSNP, EA.y, Absorption_SBP_PES)
Merged_ABSORPTION_SBP <- unique(Merged_ABSORPTION_SBP)

write.table(Merged_ABSORPTION_SBP, file="UKBB_PES/Absorption_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)


## For sensitivity analyses of int - transport PES at other thresholds

## P < 1
TRANSPORT_SBP_weights_1 <- fread("Na_K_PES_construction/PES_weights_CHR_BP/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_1_SBP_PES_HCS_weights.txt", header = T)
Merged_TRANSPORT_SBP_1 <- merge(SBP_sumstats, TRANSPORT_SBP_weights_1, by ="SNP")
Merged_TRANSPORT_SBP_1$Transport_SBP_PES_1 <- Merged_TRANSPORT_SBP_1$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_1_SBP

Merged_TRANSPORT_SBP_1 <- Merged_TRANSPORT_SBP_1 %>% select(dbSNP, EA.y, Transport_SBP_PES_1)
Merged_TRANSPORT_SBP_1 <- unique(Merged_TRANSPORT_SBP_1)

write.table(Merged_TRANSPORT_SBP_1, file="UKBB_PES/Transport_PES_sensitivity_scores/1_Transport_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)

## P < 0_5
TRANSPORT_SBP_weights_0_5 <- fread("Na_K_PES_construction/PES_weights_CHR_BP/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP_PES_HCS_weights.txt", header = T)
Merged_TRANSPORT_SBP_0_5 <- merge(SBP_sumstats, TRANSPORT_SBP_weights_0_5, by ="SNP")
Merged_TRANSPORT_SBP_0_5$Transport_SBP_PES_0_5 <- Merged_TRANSPORT_SBP_0_5$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.5_SBP

Merged_TRANSPORT_SBP_0_5 <- Merged_TRANSPORT_SBP_0_5 %>% select(dbSNP, EA.y, Transport_SBP_PES_0_5)
Merged_TRANSPORT_SBP_0_5 <- unique(Merged_TRANSPORT_SBP_0_5)

write.table(Merged_TRANSPORT_SBP_0_5, file="UKBB_PES/Transport_PES_sensitivity_scores/0_5_Transport_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)


## P < 0_05
TRANSPORT_SBP_weights_0_05 <- fread("Na_K_PES_construction/PES_weights_CHR_BP/TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.05_SBP_PES_HCS_weights.txt", header = T)
Merged_TRANSPORT_SBP_0_05 <- merge(SBP_sumstats, TRANSPORT_SBP_weights_0_05, by ="SNP")
Merged_TRANSPORT_SBP_0_05$Transport_SBP_PES_0_05 <- Merged_TRANSPORT_SBP_0_05$TRANSPORT_sodium_AND_potassium_pathways_genes_sumstats_SBP_0.05_SBP

Merged_TRANSPORT_SBP_0_05 <- Merged_TRANSPORT_SBP_0_05 %>% select(dbSNP, EA.y, Transport_SBP_PES_0_05)
Merged_TRANSPORT_SBP_0_05 <- unique(Merged_TRANSPORT_SBP_0_05)

write.table(Merged_TRANSPORT_SBP_0_05, file="UKBB_PES/Transport_PES_sensitivity_scores/0_05_Transport_SBP_PES_weights_dbSNP.txt",
            sep = "\t", row.names = F, quote = F, col.names = T)
