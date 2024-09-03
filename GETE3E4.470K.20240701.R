#!/share/apps/R-3.6.1/bin/R

# script to add E3/E4 APOE dose to 470K data and then run linear regression with DEMMax2 for funsies
# E4DOSE=L7
# E3DOSE=2-L7-L8

APOEFile <- '/home/lgibbons/APOE/APOE470Kgenotypes.txt'
APOEGenotypes <- read.table(APOEFile, header = TRUE)
DEMFile <- '/home/lgibbons/FAMDEM/UKBB.allFAMDEM.20240528.txt'
DEM <- read.table(DEMFile, header = TRUE, sep = "\t")

#Renaming columns
colnames(APOEGenotypes)[7] <- 'rs429358' #19:44908684:T:C_C=rs429358
colnames(APOEGenotypes)[8] <- 'rs7412' #19:44908822:C:T_T=rs7412
#or colnames(APOEGenotypes)<- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE", "rs429358", "rs7412")


#Adding new columns
APOEGenotypes$E3DOSE <- (2 - APOEGenotypes$rs429358 - APOEGenotypes$rs7412)
APOEGenotypes$E4DOSE <- (APOEGenotypes$rs429358)


#Merge tables
APOEDEMData <- merge(APOEGenotypes, DEM, by = "IID")
write.table(APOEDEMData, file = "APOEDEMData.470K.csv", col.names=TRUE,row.names = FALSE)

#run linear regression
model <- lm(DEMMax2 ~ E4DOSE, data = APOEDEMData)
summary(model)


