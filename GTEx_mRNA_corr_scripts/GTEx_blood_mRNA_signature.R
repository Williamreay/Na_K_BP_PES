############################################

## Transcriptomic signature of the different scores in blood (GTEx)

## William Reay (2022)

############################################

setwd("~/Documents/Na_K_BP_Erin_Clare/")

library(dplyr)
library(data.table)
library(ggplot2)
library(grex)
library(ggpubr)

## Load in mRNA ~ PES correlation

RAW_PES <- fread("GTEx_blood_mRNA_corr/Aggregated_mRNA_results/PES_mRNA_RAW_aggregated_results.txt")
RAW_PES <- RAW_PES %>% select(Gene, Score, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)

## Check for duplicates - none present
RAW_PES <- RAW_PES %>% filter(duplicated(Gene)==FALSE)

RAW_PES$FWER <- p.adjust(RAW_PES$`Pr(>|t|)`, method = "bonferroni")
RAW_PES$FDR <- p.adjust(RAW_PES$`Pr(>|t|)`, method = "fdr")

Corr_input_PES <- RAW_PES

write.table(RAW_PES, file="GTEx_blood_mRNA_corr/FDR_corr_results/PES_mRNA_signature_FDR_corr.txt",
            sep = "\t", row.names = F)

## Load in mRNA ~ PGS correlation

RAW_PGS <- fread("GTEx_blood_mRNA_corr/Aggregated_mRNA_results/PGS_mRNA_RAW_aggregated_results.txt")
RAW_PGS <- RAW_PGS %>% filter(Score == "PGS")
RAW_PGS <- RAW_PGS %>% select(Gene, Score, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)

## Check for duplicates - none present
RAW_PGS <- RAW_PGS %>% filter(duplicated(Gene)==FALSE)

RAW_PGS$FWER <- p.adjust(RAW_PGS$`Pr(>|t|)`, method = "bonferroni")
RAW_PGS$FDR <- p.adjust(RAW_PGS$`Pr(>|t|)`, method = "fdr")

Corr_input_PGS <- RAW_PGS

write.table(RAW_PGS, file="GTEx_blood_mRNA_corr/FDR_corr_results/PGS_mRNA_signature_FDR_corr.txt",
            sep = "\t", row.names = F)

## Test correlation between t values

cor.test(Corr_input_PES$`t value`, Corr_input_PGS$`t value`)

Plot_in_PES_PGS <- merge(Corr_input_PES, Corr_input_PGS, by = "Gene")

PES_PGS_SP <- ggplot(data = Plot_in_PES_PGS, aes(x=`t value.x`, y=`t value.y`)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  theme_bw() +
  xlab("mRNA ~ Transport PES t value") +
  ylab("mRNA ~ SBP PGS t value") +
  ggtitle("Transcriptome wide")

## Paired t test

d_full <- with(Plot_in_PES_PGS, 
          abs(`t value.x`) - abs(`t value.y`))

Nom <- Plot_in_PES_PGS %>% filter(`Pr(>|t|).x` < 0.05 | `Pr(>|t|).y` < 0.05)



T1 <- t.test(abs(Plot_in_PES_PGS$`t value.x`), abs(Plot_in_PES_PGS$`t value.y`), data = Plot_histo_in, paired = T)

##t = 27.045, df = 20246, p-value < 2.2e-16

Full_plot_1 <- Plot_in_PES_PGS %>% select(`t value.x`)
Full_plot_1 <- rename(Full_plot_1, "t value"="t value.x")
Full_plot_1$Score <- "Na+/K+ Transport PES"
Full_plot_2 <- Plot_in_PES_PGS %>% select(`t value.y`)
Full_plot_2 <- rename(Full_plot_2, "t value"="t value.y")
Full_plot_2$Score <- "SBP PGS"

Full_histo_in <- rbind(Full_plot_1, Full_plot_2)

KE_1 <- ggplot(Full_histo_in, aes(x=abs(`t value`), fill = Score)) +
  geom_density(alpha=0.4) +
  theme(legend.position = "bottom") +
  xlab("Absolute value of regression t value") +
  ylab("Density") +
  scale_fill_brewer(palette="Accent") +
  theme_bw()

## Assign genes HGNC IDs via "grex" package

data("gtexv7")

ID_df <- grex(gtexv7)
tail(ID_df)

Plot_in_PES_PGS$ensembl_id <- cleanid(gtex_id = Plot_in_PES_PGS$Gene)

Annot_out_merged <- merge(ID_df, Plot_in_PES_PGS, by = "ensembl_id")

## Look within PES gene-set itself

Transport_PES <- fread("Na_K_PES_construction/Gene_lists/TRANSPORT_sodium_AND_potassium_pathways_genes.txt")
Transport_PES$hgnc_symbol <- Transport_PES$TRANSPORT_sodium_and_potassium

Transport_PES_merged <- merge(Transport_PES, Annot_out_merged, by="hgnc_symbol")

Transport_plot_1 <- Transport_PES_merged %>% select(`t value.x`)
Transport_plot_1 <- rename(Transport_plot_1, "t value"="t value.x")
Transport_plot_1$Score <- "Na+/K+ Transport PES"
Transport_plot_2 <- Transport_PES_merged %>% select(`t value.y`)
Transport_plot_2 <- rename(Transport_plot_2, "t value"="t value.y")
Transport_plot_2$Score <- "SBP PGS"

Plot_histo_in <- rbind(Transport_plot_1, Transport_plot_2)

KE <- ggplot(Plot_histo_in, aes(x=abs(`t value`), fill = Score)) +
  geom_density(alpha=0.4) +
  theme(legend.position = "bottom") +
  xlab("Absolute value of regression t value") +
  ylab("Density") +
  scale_fill_brewer(palette="Dark2") +
  theme_bw()

Pathway_PGS_SP <- ggplot(data = Transport_PES_merged, aes(x=`t value.x`, y=`t value.y`)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", colour = "red") +
  theme_bw() +
  xlab("mRNA ~ Transport PES t value") +
  ylab("mRNA ~ SBP PGS t value") +
  ggtitle("Na+/K+ transport gene-set")

cor.test(Transport_PES_merged$`t value.x`, Transport_PES_merged$`t value.y`)

ggarrange(PES_PGS_SP, KE_1, Pathway_PGS_SP, KE, nrow = 2, ncol = 2)

## Test for statistically significant difference in pathway association

## Check if difference is normally distributed
d <- with(Plot_histo_in, 
          abs(`t value`)[Score == "Transport PES"] - abs(`t value`)[Score == "SBP PGS"])
shapiro.test(d)

T2 <- t.test(abs(`t value`) ~ Score, data = Plot_histo_in, paired = T)

## t = -4.5562, df = 286, p-value = 7.721e-06

## Urinary Na+ PGS comparison

Urinary_PGS <- fread("GTEx_blood_mRNA_corr/Aggregated_mRNA_results/Urinary_PGS_mRNA_aggregated.txt")
Urinary_PGS <- Urinary_PGS %>% select(Gene, Score, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`)

## Merge with other df

Final_full_all_three_scores <- merge(Urinary_PGS, Plot_in_PES_PGS, by ="Gene")

## PES ~ Urinary

cor.test(Final_full_all_three_scores$`t value`, Corr_input_PES$`t value`)

## PGS ~ urinary

cor.test(Final_full_all_three_scores$`t value`, Corr_input_PGS$`t value`)

plot(Final_full_all_three_scores$`t value`, Corr_input_PGS$`t value`)

## Subset by the transport pathway


Final_full_all_three_scores$ensembl_id <- cleanid(gtex_id = Final_full_all_three_scores$Gene)

All_three_merge <- merge(ID_df, Final_full_all_three_scores, by = "ensembl_id")

ALL_Transport_PES_merged <- merge(Transport_PES, All_three_merge, by="hgnc_symbol")

cor.test(ALL_Transport_PES_merged$`t value`, ALL_Transport_PES_merged$`t value.x`)
cor.test(ALL_Transport_PES_merged$`t value`, ALL_Transport_PES_merged$`t value.y`)

