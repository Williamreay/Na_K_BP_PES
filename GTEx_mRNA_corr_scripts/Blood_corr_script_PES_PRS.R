#####################################

## Correlating PES and PRS with normalised mRNA in GTEx v8

## Whole Blood

## William Reay (2021)

#####################################

suppressMessages(library(data.table))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

##Specify command line inputs
option_list = list(
  make_option("--score_name", action="store", default=NA, type='character',
              help="Name of PGS or PES [required]"),
  make_option("--score_col_id", action="store", default=NA, type='character',
              help="Column id for the PGS or PES [required]"),
  make_option("--tissue_name", action="store", default=NA, type='character',
              help="Name of the tissue for the PGS or PES [required]"),
  make_option("--gene_name", action="store", default=NA, type='character',
              help="Name of the gene to be tested")
  
)

opt = parse_args(OptionParser(option_list=option_list))

print(opt)


setwd("~/data/users/william/Na_K_GTEx_corr/")

## Read in mRNA data

mRNA_dat <- fread("Formatted_blood_mRNA_df_cov_incl.txt")

## Read in scores

Score_dat_1 <- fread("SBP_PGS_transport_PES_GTEx.txt")

## Merged

Merged_df_1 <- merge(Score_dat_1, mRNA_dat, by = "IID")

rm(mRNA_dat)

List_genes <- as.list(opt$gene_name)

## Function for correlation

mRNA_PES_PRS <- function(gene, score_col, df) {
  cat("\n")
  cat(paste("Starting analysis of effect of ", score_col, " on transcript ",  gene, sep=""))
  cat("\n")
  df$scaled_score <- as.numeric(scale(df[[score_col]]))
  mod <- lm(glue::glue('{gene} ~ sex + PC1 + PC2 + PC3 + PC4 + PC5 +
            InferredCov1 + InferredCov2 + InferredCov3 + InferredCov4 + InferredCov5 +
            InferredCov6 + InferredCov7 + InferredCov8 + InferredCov9 + InferredCov10 + pcr + platform + scaled_score'), data = df)
  return(summary(mod))
}

Run_lm <- sapply(List_genes, mRNA_PES_PRS,  score_col=opt$score_col_id, df=Merged_df_1)

Run_lm_extract <- apply(Run_lm, 2, function(x) return(as.data.frame(x$coefficients)[20, 1:4]))

Output <- data.frame()
for (i in 1:length(Run_lm_extract)) {
  Output <- rbind(Output, Run_lm_extract[[i]])
}

Output$Gene <- unlist(List_genes)

Output$FDR <- p.adjust(Output$'Pr(>|t|)', method="fdr")
Output$Score <- opt$score_name

cat("\n")
cat("#########################")
cat("\n")
cat("Writing output")
cat("\n")
cat("#########################")
cat("\n")

write.table(Output, file = paste("Gene_corr/BLOOD_", opt$score_name, "_", opt$gene_name, ".txt", sep=""), sep = "\t", row.names = F, quote = F)

#Clear environment after run
rm(list = ls())





