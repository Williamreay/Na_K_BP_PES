#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 2b - exploring transport PES x sodium interaction - UNMEDICATED ONLY

## Effects per percentile - test different curves

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
library(viridis)
library(car)

Subcohort_1_dat <- fread("Creatinine_adjusted_subcohort_1.txt")

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

PGS_BP <- fread("UKBB_PGS/SBP_DBP_PGS_UKBB_exTrained.txt")

Merged_PGS_PES <- merge(PGS_BP, Merged_PES, by = "IID")
Merged_PGS_PES$Transport_SBP_PES_AVG <- as.numeric(scale(Merged_PGS_PES$Transport_SBP_PES_AVG))
Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG <- as.numeric(scale(Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG))

## Test effect of scaled UNa:UCr, Na:K ratio and scaled Na
Merged_PGS_PES$Scaled_Na <- as.numeric(scale(Merged_PGS_PES$Na_urine_millimolL))
Merged_PGS_PES$Scaled_Na_Cr <- as.numeric(scale(Merged_PGS_PES$Na_Cr))
Merged_PGS_PES$Scaled_Na_K <- as.numeric(scale(Merged_PGS_PES$Na_K_ratio))

## Calculate the effect of sodium per score percentile

Merged_PGS_PES$Percentile_PES <- as.numeric(ntile(Merged_PGS_PES$Transport_SBP_PES_AVG, 100))
Merged_PGS_PES$Percentile_PGS <- as.numeric(ntile(Merged_PGS_PES$SBP_PGS_Pt_0_01_AVG, 100))

## Simplified model as per percentile
Na_assoc_per_percentile <- function(percentile, urinary, score, df) {
  new_df <- df[df[[score]] == percentile,]
  new_df$scaled_urinary <- as.numeric(scale(new_df[[urinary]]))
  mod <- lm(Mean_SBP ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + scaled_urinary, data = new_df)
  return(summary(mod))
}

Vec_percentile <- c(1:100)
Vec_percentile <- rep(Vec_percentile, each=3)

Trait_names <- as.list(c("Na_Cr", "Na_urine_millimolL", "Na_K_ratio"))

Trait_names_in <- rep(Trait_names, 100)

## Run models for UNa:UCr, UNa, and UNa:UK

PES <- mapply(Na_assoc_per_percentile, percentile = Vec_percentile, urinary = Trait_names_in, 
              MoreArgs = list("Percentile_PES", Merged_PGS_PES), SIMPLIFY = T)

PGS <- mapply(Na_assoc_per_percentile, percentile = Vec_percentile, urinary = Trait_names_in, 
              MoreArgs = list("Percentile_PGS", Merged_PGS_PES), SIMPLIFY = T)


PES_extract <- apply(PES, 2, function(x) return(as.data.frame(x$coefficients)[10, 1:4]))
PRS_extract <- apply(PGS, 2, function(x) return(as.data.frame(x$coefficients)[10, 1:4]))

## Df for each

PES_plot <- data.frame()

for (i in 1:length(PES_extract)) {
  PES_plot <- rbind(PES_plot, PES_extract[[i]])
}

PES_plot$Outcome <- unlist(Trait_names_in)
PES_plot$Percentile <- Vec_percentile

## Group percentile values together into < 20th, 20th-80th, and >80th
PES_plot$Group <- ifelse(PES_plot$Percentile < 21, "1st-20th percentile", 0)
PES_plot$Group <- ifelse(PES_plot$Percentile > 20 & PES_plot$Percentile < 80, "20th-80th percentile", PES_plot$Group)
PES_plot$Group <- ifelse(PES_plot$Percentile > 79, "80th-100th percentile", PES_plot$Group)

PES_plot$UCI <- PES_plot$Estimate + (1.96*PES_plot$`Std. Error`)
PES_plot$LCI <- PES_plot$Estimate - (1.96*PES_plot$`Std. Error`)

PRS_plot <- data.frame()

for (i in 1:length(PRS_extract)) {
  PRS_plot <- rbind(PRS_plot, PRS_extract[[i]])
}

PRS_plot$Percentile <- Vec_percentile
PRS_plot$Outcome <- unlist(Trait_names_in)
PRS_plot$Group <- ifelse(PRS_plot$Percentile < 21, "1st-20th percentile", 0)
PRS_plot$Group <- ifelse(PRS_plot$Percentile > 20 & PRS_plot$Percentile < 80, "20th-80th percentile", PES_plot$Group)
PRS_plot$Group <- ifelse(PRS_plot$Percentile > 79, "80th-100th percentile", PRS_plot$Group)

PRS_plot$UCI <- PRS_plot$Estimate + (1.96*PRS_plot$`Std. Error`)
PRS_plot$LCI <- PRS_plot$Estimate - (1.96*PRS_plot$`Std. Error`)

## Test stat sig effect of top decile vs the rest

## First extract each outcome trait individually

Na_Cr_plot <- PES_plot %>% filter(Outcome == "Na_Cr")
Na_plot <- PES_plot %>% filter(Outcome == "Na_urine_millimolL")
Na_K_plot <- PES_plot %>% filter(Outcome == "Na_K_ratio")

## PGS
PGS_Na_Cr_plot <- PRS_plot %>% filter(Outcome == "Na_Cr")
PGS_Na_plot <- PRS_plot %>% filter(Outcome == "Na_urine_millimolL")
PGS_Na_K_plot <- PRS_plot %>% filter(Outcome == "Na_K_ratio")


## Make plots showing gam curve for effect per decile

library(ggsci)

## Na_Cr_plot
mean(Na_Cr_plot$Estimate)
mean(Na_Cr_plot$Estimate)-sd(Na_Cr_plot$Estimate)
mean(Na_Cr_plot$Estimate)+sd(Na_Cr_plot$Estimate)

P1_Na_Cr <- ggplot(Na_Cr_plot, aes(Percentile, Estimate)) +
  #stat_smooth(method = "loess", se = F, size = 0.5, colour = "black") +
  stat_smooth(method = "loess", se = T, size = 0.5, colour = "steelblue") +
  #geom_smooth(method = "gam", formula = y ~ s(x, bs="cs"), colour="steelblue") +
  theme_bw() +
  geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 2.088, lty = "dashed") +
  ggtitle("UNa+:UCr")

P1_Na_Lm <- ggplot(Na_plot, aes(Percentile, Estimate, ymax=UCI, ymin=LCI)) +
  geom_smooth(method = "lm", formula = y ~ x, colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 1.26, lty = "dashed") +
  ggtitle("GLM")

P2_Na_Quad <- ggplot(Na_plot, aes(Percentile, Estimate, ymax=UCI, ymin=LCI)) +
  geom_smooth(method = "lm", formula=y~x+I(x^2), colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 1.26, lty = "dashed") +
  ggtitle("GLM - Quadratic")

P3_Na_GAM <- ggplot(Na_plot, aes(Percentile, Estimate, ymax=UCI, ymin=LCI)) +
  geom_smooth(method = "loess", colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 1.26, lty = "dashed") +
  ggtitle("LOESS")

Out1 <- ggarrange(P1_Na_Lm, P2_Na_Quad, P3_Na_GAM, nrow=1)

LM <- lm(Estimate ~ Percentile, data = Na_plot)
LM_2 <- lm((Estimate/`Std. Error`) ~ Percentile, data = Na_plot)

write.table(Na_plot, file="Urinary_Na_K_divided_by_creatinine/Final_percentile_results_and_plots/Na_per_PES_percentile.txt",
            sep = "\t", row.names = F, quote = F)


Na_plot$Z <- Na_plot$Estimate/Na_plot$`Std. Error`

## Make boxplots for grouped effects

BP_1 <- ggplot(Na_plot, aes(x=Group, y=Estimate, fill=Group)) +
  geom_boxplot() +
  geom_point(size = 1.2) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size=11)) +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("") +
  labs(y=expression("SBP"~beta~"per SD in UNa+")) +
  geom_hline(yintercept = 1.26, lty = "dashed") +
  ggtitle(expression(beta~"coefficient"))

BP_2 <- ggplot(Na_plot, aes(x=Group, y=Z, fill=Group)) +
  geom_boxplot() +
  geom_point(size = 1.2) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size=11)) +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("") +
  labs(y=expression("SBP Z per SD in UNa+")) +
  geom_hline(yintercept = 3.45, lty = "dashed") +
  ggtitle("Z-value")

BP <- ggarrange(BP_1, BP_2)

ggarrange(BP, Out1, nrow=2)

## Make correlation plot of Na effect vs Na:Cr effect and then show some Na:Cr effects

Comp <- merge(Na_Cr_plot, Na_plot, by = "Percentile")


C1 <- ggplot(Comp, aes(Estimate.y, Estimate.x, colour=Group.x)) +
  geom_point() +
  geom_smooth(method = "lm", colour="firebrick") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "bottom") +
  ylab("Per PES percentile effect - UNa+:UCr") +
  xlab("Per PES percentile effect - UNa+") +
  ggtitle("PES stratification: UNa+ vs UNa+:UCr") +
  scale_colour_brewer(palette = "Paired") +
  labs(colour="Group")

C2 <- ggplot(Na_Cr_plot, aes(x=Group, y=Estimate, fill=Group)) +
  geom_boxplot() +
  geom_point(size = 1.2) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size=11)) +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("") +
  labs(y=expression("SBP"~beta~"per SD in UNa+:UCr")) +
  geom_hline(yintercept = 2.089, lty = "dashed") +
  ggtitle("UNa+:UCr PES percentiles")

## Corr between is 0.3238291 [95% CI: 0.1360656, 0.4891353]
cor.test(Comp$Estimate.x, Comp$Estimate.y)
Comp$Z_Na_cr <- Comp$Estimate.x/Comp$`Std. Error.x`

C3 <- ggplot(Comp) +
  geom_smooth(aes(x=Percentile, y=Z), method="loess", colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("SBP Z per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  ggtitle("UNa+") +
  geom_hline(yintercept = 3.44, lty="dashed")

C4 <- ggplot(Comp) +
  geom_smooth(aes(x=Percentile, y=Z_Na_cr), method="loess", colour="darkorchid") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("SBP Z per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  ggtitle("UNa+:UCr") +
  geom_hline(yintercept = 5.87, lty="dashed")

ggarrange(C1, C2, C3, C4)

## Tests for diff

pairwise.t.test(Na_Cr_plot$Estimate, Na_Cr_plot$Group, p.adjust.method = "none")

pairwise.t.test(Na_plot$Estimate, Na_plot$Group, p.adjust.method = "none")

## Plot for supplementary for PES vs PGS for both

Comp_2 <- merge(Na_Cr_plot, PGS_Na_Cr_plot, by="Percentile")

PG1 <- ggplot(Comp_2, aes(Percentile, Estimate.x)) +
  geom_smooth(method = "loess", colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 1.26, lty = "dashed") +
  ggtitle("LOESS")

PG2 <- ggplot(Comp_2, aes(Percentile, Estimate.y)) +
  geom_smooth(method = "loess", colour="steelblue") +
  theme_bw() +
  #geom_point() +
  theme(legend.position = "none") +
  ylab("Effect on SBP (mmHg) per SD in UNa+") +
  xlab("Transport PES (percentiles)") +
  geom_hline(yintercept = 2.2, lty = "dashed") +
  ggtitle("LOESS")


