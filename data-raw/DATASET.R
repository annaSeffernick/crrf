#########################################################
# Code to prepare ACT Data
# From Clay et al. 2019 Precision Oncology
# https://doi.org/10.1200/PO.19.00163
# AES
# 2024-07-02
###########################################################

# Set up
library(openxlsx)

# Load data
# supplementary table
act.dat <- read.xlsx("Y:/Anna/CRR/Data/Clay-JCO-PO-Supplement-DS_po.19.00163-1.xlsx",  rows = 3:51)

# full dataset from Stan Pounds from original project
full.dat <- read.xlsx("Y:/Anna/CRR/Data/Adrenocortical Tumor Master Clinical Excel 1.30.2019 no MRNs for Anna.xlsx",
                      sheet="ACT Active Excel 6.13.2018", detectDates=TRUE)

# Compare the datasets
length(which(act.dat$Sentrix_ID %in% full.dat$Sentrix_ID)) # all 48 act cases in full.dat

# Calculate time to relapse for the patients that relapsed
rel.ids <- act.dat[which(act.dat$Relapse_YN=="YES"),]$Sentrix_ID # ids of patients who relapsed
rel.dat <- full.dat[which(full.dat$Sentrix_ID %in% rel.ids),] # filter the full dataset to relapse pts
rel.dat$Relapse_Time <- difftime(as.Date(rel.dat$Relapse_Date,format="%Y-%m-%d"), as.Date(rel.dat$Collection.Date, format="%Y-%m-%d"), unit="days")

# Create variables needed for crrf package
act.dat$RFS <- act.dat$Follow_Up_Days # Define relapse-free survival variable: time until relapse or death w/o relapse or last follow-up if censored
rel.dat.or <- rel.dat[match(rel.ids, rel.dat$Sentrix_ID),] # reorder rel.dat
act.dat[which(act.dat$Relapse_YN=="YES"),]$RFS <- rel.dat.or$Relapse_Time # replase follow-up time with time to relapse for those who relapse
act.dat$RFS_status <- act.dat$Mortality
act.dat[which(act.dat$Relapse_YN=="YES"),]$RFS_status <- "Relapse"
act.dat$RFS_status <- factor(act.dat$RFS_status)
table(act.dat$RFS_status)

# Modify some colnames
colnames(act.dat) <- c("Onco_Case", "Sentrix_ID", "Gender",
                       "Age_years", "Clinical_Presentation",
                       "Stage", "Wieneke_Classification",
                       "Initial_Treatment", "Relapse_YN",
                       "Mortality", "Deceased_From_Disease",
                       "Follow_Up_Days", "Tumor_Weight_Grams",
                       "Tumor_Size_cm", "Ki67_Percent",
                       "P53_IHC", "P53_IHC_Percent",
                       "BC_IHC", "BC_IHC_Percent",
                       "ATRX_IHC", "TP53_Mutation_Category",
                       "CTNNB1_Mutation_Category",
                       "ATRX_Mutation_Category",
                       "X11p15_status", "Methylation_Group",
                       "RFS", "RFS_status")


usethis::use_data(act.dat, overwrite = TRUE)
