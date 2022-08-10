############################################

## Pre-processing files for Erin and Clare

## SBP/DBP

## William Reay (2022)

############################################

library(data.table)
library(dplyr)
library(easyGgplot2)
library(table1)

## Read in basic demographics and cardiac data

setwd("~/Documents/Na_K_BP_Erin_Clare/")

Demo <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Basic_demographics_UKBB_July_2021_Freeze.txt")

Cardiac <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Quantitative_cardiovascular_data_UKBB_July_2021_Freeze.txt")

Diet_FFQ <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Diet_physical_activity_UKBB_July_2021_Freeze.txt")
  
Diet_additional <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Erin_diet.txt")

Biochem <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Biochem_immunology_data_UKBB_July_2021_Freeze.txt")

Merged_df1 <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "eid"),
                     list(Demo, Cardiac, Diet_FFQ, Diet_additional, Biochem))

## Retain SBP and DBP

Merged_BP <- Merged_df1 %>% select(eid, `4080-0.0`, `4080-0.1`,
                                   `4079-0.0`, `4079-0.1`)

## Retain only individuals with non-missing measurements

Merged_BP <- Merged_BP %>% filter(!is.na(`4080-0.0`) & !is.na(`4080-0.1`) &
                                    !is.na(`4079-0.0`) & !is.na(`4079-0.1`))

## Take mean SBP and DBP

Merged_BP$Mean_SBP <- ((Merged_BP$`4080-0.0` + Merged_BP$`4080-0.1`)/2)

Merged_BP$Mean_DBP <- ((Merged_BP$`4079-0.0` + Merged_BP$`4079-0.1`)/2)

Merged_BP <- rename(Merged_BP, "IID"="eid")

Merged_BP <- Merged_BP %>% select(IID, Mean_SBP, Mean_DBP)

## Merge with EUR IDs with genetic data available

Gen_covar <- fread("~/Desktop/Blood_pressure_target_ranking/Blood_pressure_UKBB_vQTL/Covariates_pneumonia_UKBB_white_British.tab.txt")

Gen_merged <- merge(Merged_BP, Gen_covar, by = "IID")

## Derive blood pressure medication data (6153 = females, 6177 = males)

Med_self_report <- fread("/Volumes/Wills_new_backup/UKBB_July_2021_version/Subsetted_phenotype_fields/Self_report_blood_pressure_meds_at_basline.txt")

select_Meds_1<- Med_self_report %>% select(starts_with("6153-0"))

select_Meds_2 <- Med_self_report %>% select(starts_with("6177-0"))

select_eid <- Med_self_report %>% select(eid)

Meds_df <- cbind(select_eid, select_Meds_1, select_Meds_2)

## Extract code 2

Med_select <- which(Meds_df==2, arr.ind = T)

Med_select_row <- Med_select[,1]

Med_select_row <- unique(Med_select_row)

Meds_df$BP_meds <- 0
Meds_df[Med_select_row, ]$BP_meds <- 1

table(Meds_df$BP_meds)

Meds_df <- Meds_df %>% select(eid, BP_meds)

Meds_df <- rename(Meds_df, "IID"="eid")

## Merge

Gen_meds_merge <- merge(Gen_merged, Meds_df, by = "IID")

Gen_meds_merge$BP_Medication <- ifelse(Gen_meds_merge$BP_meds == 1, "YES", "NO")


## Select relevant biochem data from df

Subset_biochem <- Biochem %>% select(eid, `30510-0.0`, `30512-0.0`, `30513-0.0`,
                                     `30520-0.0`, `30522-0.0`, `30523-0.0`,
                                     `30530-0.0`, `30532-0.0`, `30533-0.0`)

Subset_biochem <- na.omit(Subset_biochem)

Subset_biochem <- rename(Subset_biochem, "Creatinine_urine_micromolL"=`30510-0.0`,
                         "Creatinine_urine_acqusition_time"=`30512-0.0`, "Creatinine_device_ID"=`30513-0.0`,
                         "K_urine_millimolL"=`30520-0.0`, "K_acquisition_time"=`30522-0.0`, "K_device_ID"=`30523-0.0`,
                         "Na_urine_millimolL"=`30530-0.0`, "Na_acquisition_time"=`30532-0.0`, "Na_device_ID"=`30533-0.0`)


BMI <- Demo %>% select(eid, `21001-0.0`)

BMI_df <- BMI %>% filter(!is.na(`21001-0.0`))

BMI_df <- rename(BMI_df, "BMI"=`21001-0.0`)

## Biochem_BP_genetics only

Gen_meds_merge$eid <- Gen_meds_merge$IID

Biochem_BP_gen_out <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "eid"),
       list(Gen_meds_merge, BMI_df, Subset_biochem))

write.table(Biochem_BP_gen_out, file = "Biochem_BP_genetics_EUR.txt",sep = "\t",
            row.names = F, quote = F)


Histo_Na <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="Na_urine_millimolL",
                                    fill = "white", color = "black", addMeanLine = T,
                                    meanLineType = "dashed", meanLineColor = "black",
                                    addDensityCurve = TRUE, densityFill = '#33FF66',
                                    xtitle = "Urinary Sodium (millimole/L)", bins = 50,
                                    ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                                    axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                                    yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                                    mainTitle = "Urinary sodium (spot) - N = 296,527")


Histo_K <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="K_urine_millimolL",
                              fill = "white", color = "black", addMeanLine = T,
                              meanLineType = "dashed", meanLineColor = "black",
                              addDensityCurve = TRUE, densityFill = '#FF9900',
                              xtitle = "Urinary Potassium (millimole/L)", bins = 50,
                              ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                              axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                              yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                              mainTitle = "Urinary potassium (spot) - N = 296,527")



Histo_meds_SBP <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="Mean_SBP",
                                    fill = "white", color = "black", addMeanLine = T,
                                    meanLineType = "dashed", meanLineColor = "black",
                                    addDensityCurve = TRUE, densityFill = '#FF3333',
                                    xtitle = "Raw SBP (mmHg)", bins = 50,
                                    ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                                    axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                                    yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                                    mainTitle = "SBP - N = 296,527", groupName = "BP_medication", groupColors=c('#999999','#FF3333'),
                                    legendTitle ="Blood pressure medication")


Histo_meds_DBP <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="Mean_DBP",
                                    fill = "white", color = "black", addMeanLine = T,
                                    meanLineType = "dashed", meanLineColor = "black",
                                    addDensityCurve = TRUE, densityFill = '#FF3333',
                                    xtitle = "Raw DBP (mmHg)", bins = 50,
                                    ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                                    axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                                    yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                                    mainTitle = "DBP - N = 296,527", groupName = "BP_medication", groupColors=c('#999999','#0066FF'),
                                    legendTitle ="Blood pressure medication")

Histo_raw_SBP <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="Mean_SBP",
                                   fill = "white", color = "black", addMeanLine = T,
                                   meanLineType = "dashed", meanLineColor = "black",
                                   addDensityCurve = TRUE, densityFill = '#FF3333',
                                   xtitle = "Raw SBP (mmHg)", bins = 50,
                                   ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                                   axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                                   yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                                   mainTitle = "SBP - N = 296,527")

Histo_raw_DBP <- ggplot2.histogram(data = Biochem_BP_gen_out, xName="Mean_DBP",
                                   fill = "white", color = "black", addMeanLine = T,
                                   meanLineType = "dashed", meanLineColor = "black",
                                   addDensityCurve = TRUE, densityFill = '#0066FF',
                                   xtitle = "Raw DBP (mmHg)", bins = 50,
                                   ytitle="Smoothed Density", xtitleFont=c(12, "bold", "black"), ytitleFont=c(12, "bold", "black"), 
                                   axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                                   yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE,
                                   mainTitle = "DBP - N = 296,527")

table1( ~  factor(`BP Medication`) + Age + Mean_SBP + Mean_DBP + BMI + Na_urine_millimolL +
          K_urine_millimolL + Creatinine_urine_micromolL | factor(Sex),
        data = Biochem_BP_gen_out)

## Select relevant dietary data

Diet_FFQ <- Diet_FFQ %>% select(eid, starts_with("1000")) %>% select(eid, ends_with("0.0"))

Diet_FFQ_cleaned <- na.omit(Diet_FFQ)

Erin_additional <- Diet_additional %>% select(eid, `1289-0.0`, `1299-0.0`,
                                              `1309-0.0`, `1319-0.0`, `1478-0.0`,
                                              `1558-0.0`)
Erin_additional_cleaned <- na.omit(Erin_additional)

Merged_diet <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "eid"),
                      list(Biochem_BP_gen_out, Diet_FFQ_cleaned, Erin_additional_cleaned))

write.table(Merged_diet, file ="All_phenotypes_including_diet.txt",
            sep = "\t", row.names = F, quote = F)

