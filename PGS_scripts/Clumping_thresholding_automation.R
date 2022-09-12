##########################

## Clumping and thresholding for scoring using 1000G reference (EUR)

## William Reay (2022)

########################

library(dplyr)
library(data.table)
library(ieugwasr)
library(optparse)
library(genetics.binaRies)

setwd("~/Documents/Na_K_BP_Erin_Clare/Genetic_data/")

SBP_input <- fread("Cleaned_GERA_sumstats_SBP.txt.gz")

DBP_input <- fread("Cleaned_GERA_sumstats_DBP.txt.gz")

## Remove the MHC for scoring without MHC

MHC_removed_SBP <- SBP_input %>% filter(chromosome != "6" | (chromosome == "6" & position > 28477797 & position < 33448354))

MHC_removed_DBP <- DBP_input %>% filter(chromosome != "6" | (chromosome == "6" & position > 28477797 & position < 33448354))


HCS_Clumping_thresholding_func <- function(df, threshold, MHC, BP) {
  Sumstats <- rename(df, "rsid"="dbSNP", "pval"="P")
  C_T <- ld_clump_local(dplyr::tibble(rsid=Sumstats$rsid, pval=Sumstats$pval), 
                        clump_r2 =  0.1, clump_p = threshold, clump_kb = 250,
                        plink_bin = genetics.binaRies::get_plink_binary(),
                        bfile = "/Users/williamreay/cloudstor/Cross_disorder_PES_2019/GEUVADIS_gene_expression/g1000_eur"
  )
  Merged <- merge(C_T, Sumstats, by ="rsid")
  Merged_out <- Merged %>% select(SNP, EA, Beta)
  names(Merged_out)[names(Merged_out) == "Beta"] <- threshold
  Merged_FINAL <- unique(Merged_out)
  write.table(Merged_FINAL, file = paste0(threshold,"_",MHC,"_", BP,"_HCS_weights.txt", sep=""),
              sep = "\t", row.names = F, quote = F, col.names = T)
}

List_of_thresholds <- as.vector(c(5e-08, 1e-05, 1e-04, 1e-03, 0.05, 0.005, 0.01, 0.1, 0.5, 1))

## SBP PGS
sapply(List_of_thresholds, HCS_Clumping_thresholding_func, df = MHC_removed_SBP, MHC = "exMHC", BP = "SBP")

## DBP PGS
sapply(List_of_thresholds, HCS_Clumping_thresholding_func, df = MHC_removed_DBP, MHC = "exMHC", BP = "DBP")

