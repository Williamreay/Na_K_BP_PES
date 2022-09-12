#############################

## Processing and cleaning GERA SBP and DBP GWAS

## William Reay (2022)

#############################


setwd("~/Documents/Na_K_BP_Erin_Clare/Genetic_data/")

library(data.table)
library(dplyr)

SBP_raw <- fread("gera-sbp.tsv.gz")
names(SBP_raw)<-make.names(names(SBP_raw),unique = TRUE)
DBP_raw <- fread("gera-dbp.tsv.gz")
names(DBP_raw)<-make.names(names(DBP_raw),unique = TRUE)

setnames(SBP_raw, c("Allele.1","Allele.2", "Sample.size",
                  "Estimate.Effect", "P.value", "Effect.allele..EA.",
                  "Effect.allele.frequency..EAF."), c("A1", "A2","N", "Beta", "P", "EA",
                                                     "EAF"))

setnames(DBP_raw, c("Allele.1","Allele.2", "Sample.size",
                    "Estimate.Effect", "P.value", "Effect.allele..EA.",
                    "Effect.allele.frequency..EAF."), c("A1", "A2","N", "Beta", "P", "EA",
                                                        "EAF"))

## QC - MAF < 0.01

SBP_QC_1 <- SBP_raw %>% filter(EAF > 0.01 & EAF < 0.99)
DBP_QC_1 <- DBP_raw %>% filter(EAF > 0.01 & EAF < 0.99)

## QC duplicated variants

SBP_QC_2 <- as.data.frame(distinct(SBP_QC_1, SNP_ID, .keep_all = T))
DBP_QC_2 <- as.data.frame(distinct(DBP_QC_1, SNP_ID, .keep_all = T))

## QC ambigious SNPs

SBP_QC_3 <- SBP_QC_2 %>% filter(!(A1 =="A" & A2=="T")) %>%
                        filter(!(A1 == "T" & A2 == "A")) %>%
                        filter(!(A1 == "G" & A2 == "C")) %>%
                        filter(!(A1 == "C" & A2 == "G"))

DBP_QC_3 <- DBP_QC_2 %>% filter(!(A1 =="A" & A2=="T")) %>%
  filter(!(A1 == "T" & A2 == "A")) %>%
  filter(!(A1 == "G" & A2 == "C")) %>%
  filter(!(A1 == "C" & A2 == "G"))


## Retain only SNPs tested in > 90% of the sample (N > 89807)

SBP_QC_N <- SBP_QC_3 %>% filter(N > 89807)
DBP_QC_N <- DBP_QC_3 %>% filter(N > 89807)

## Add CHR:BP column for HCS cohort

SBP_QC_N$SNP <- paste("",SBP_QC_N$chromosome, ":", SBP_QC_N$position, sep="")
DBP_QC_N$SNP <- paste("",SBP_QC_N$chromosome, ":", DBP_QC_N$position, sep="")

SBP_QC_N <- rename(SBP_QC_N, "dbSNP"="SNP_ID")
DBP_QC_N <- rename(DBP_QC_N, "dbSNP"="SNP_ID")

## EA and A1 same 

## Output one df for HCS cohort use (CHR:BP ID)

write.table(SBP_QC_N, file="Cleaned_GERA_sumstats_SBP.txt",
            sep = "\t", row.names = F, quote = F)

write.table(DBP_QC_N, file="Cleaned_GERA_sumstats_DBP.txt",
            sep = "\t", row.names = F, quote = F)
