#################################################

## Na/K and hypertension: interplay of diet and genetics

## STAGE 2: Step 4 - Relationship with BP and interaction with PES profiles

## Medicated and unmedicated separately 

## William Reay (2022)

#################################################

## Read in files, set relevant variables as categorical

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

## Scale PES - medicated

Merged_PES[, c(2:7)] <-  lapply(Merged_PES[, c(2:7)], function(x) c(scale(x)))

## Make profile that adds up all scaled PES values

Merged_PES$PES_profile <- rowSums(Merged_PES[,c(2:7)])

## Test assoc with SBP and DBP

SBP_PES_profile_baseline <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                         Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                         PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                           PES_profile, data = Merged_PES)

DBP_PES_profile_baseline <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                 Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                 PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                 PES_profile, data = Merged_PES)

Merged_PES$Scaled_PES_profile <- as.numeric(scale(Merged_PES$PES_profile))
Merged_PES$Scaled_Na <- as.numeric(scale(Merged_PES$Na_urine_millimolL))
Merged_PES$Scaled_K <- as.numeric(scale(Merged_PES$K_urine_millimolL))

SBP_PES_profile_Na_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                                 Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                                 PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                                 Scaled_PES_profile*Scaled_Na, data = Merged_PES)

SBP_PES_profile_K_int <- lm(Mean_SBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                               Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                               PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                               Scaled_PES_profile*Scaled_K, data = Merged_PES)

DBP_PES_profile_Na_int <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                               Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                               PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                               Scaled_PES_profile*Scaled_Na, data = Merged_PES)

DBP_PES_profile_K_int <- lm(Mean_DBP ~ Sex + Age + Age2 + BP_meds + Assesment_centre_month + 
                              Assesment_centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC10 + PC11 +
                              PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 +
                              Scaled_PES_profile*Scaled_K, data = Merged_PES)
