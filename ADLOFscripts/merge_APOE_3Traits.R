#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)

# Merge Data of APOE and sex; APOE, Sex, HL
apoe_data <- read.table("UKBB.LOF.20230809.APOE.sco", header = FALSE, sep = "", fill = TRUE, stringsAsFactors = FALSE)
colnames(apoe_data) <- c("IID", "sex", "APOE")
hyperlipidemia_data <- read.table("UKBB.HL.20230915.txt", header=TRUE, sep="\t")
sex_data<- read.table("UKBB.Sex.txt", header=TRUE, sep="\t")

apoe_data_clean <- na.omit(apoe_data)
hyperlipidemia_data_clean <- na.omit(hyperlipidemia_data)
sex_data_clean <- na.omit(sex_data)

merged_APOE_Sex <- merge(sex_data_clean, apoe_data_clean, by="IID", all=FALSE)

merged_APOE_Sex1 <- merged_APOE_Sex[ ,-3]

merged_APOE_HL <- merge(merged_APOE_Sex1, hyperlipidemia_data_clean, by="IID", all=FALSE)
final_APOE_HL <- na.omit(merged_APOE_HL)


write.table(final_APOE_HL, file="Merged_APOE_Hyperlipidemia.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(merged_APOE_Sex1, file="Merged_APOE_Sex.txt", sep="\t", row.names=FALSE, quote=FALSE)

cat("Data cleaning and merging completed. File saved as Merged_APOE_Hyperlipidemia.txt\n")

#Merge data of APOE and General health; APOE and Age Edu
GeneralHealth_data <- read.table("UKBB.GeneralHealth.txt", header=TRUE, sep="\t")
AgeEdu_data <- read.table("UKBB.AgeEdu.txt", header=TRUE, sep="\t")

GeneralHealth_data_clean <- na.omit(GeneralHealth_data)
AgeEdu_data_clean <- na.omit(AgeEdu_data)

GeneralHealth_data_clean <- GeneralHealth_data_clean[rowSums(GeneralHealth_data_clean < 0) == 0, ]
AgeEdu_data_clean <- AgeEdu_data_clean[rowSums(AgeEdu_data_clean < 0) == 0, ]


merged_APOE_Health <- merge(merged_APOE_Sex1, GeneralHealth_data_clean, by="IID", all=FALSE)
merged_APOE_Edu <- merge(merged_APOE_Sex1, AgeEdu_data_clean, by="IID", all=FALSE)

final_APOE_Health <- na.omit(merged_APOE_Health)
final_APOE_Edu <- na.omit(merged_APOE_Edu)

write.table(final_APOE_Health, file="Merged_APOE_GeneralHealth.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(final_APOE_Edu, file="Merged_APOE_EducationAge.txt", sep="\t", row.names=FALSE, quote=FALSE)


#Merge the files with PCs
PC_data <- read.table("/SAN/ugi/UGIbiobank/data/downloaded/ukb23158.common.all.20230806.eigenvec.txt", header=TRUE, sep="\t")
PC_data <- PC_data[ ,-1]
Merged_APOE_Sex_PC <- merge(merged_APOE_Sex1, PC_data, by="IID", all=FALSE)
Merged_APOE_HL_PC <- merge(final_APOE_HL, PC_data, by="IID", all=FALSE)
Merged_APOE_Health_PC <- merge(final_APOE_Health, PC_data, by="IID", all=FALSE)
Merged_APOE_Edu_PC <- merge(final_APOE_Edu, PC_data, by="IID", all=FALSE)

write.table(Merged_APOE_Sex_PC, file="APOE.Sex.PC.20241023.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(Merged_APOE_HL_PC, file="APOE.HL.PC.20241023.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(Merged_APOE_Health_PC, file="APOE.Health.PC.20241023.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(Merged_APOE_Edu_PC, file="APOE.Education.PC.20241023.txt", sep="\t", row.names=FALSE, quote=FALSE)

#quick stats

#Mean years of education
mean_year_of_education <- mean(Merged_APOE_Edu_PC$AgeEdu, na.rm = TRUE)
mean_education_APOE_less_5 <- mean(Merged_APOE_Edu_PC$AgeEdu[Merged_APOE_Edu_PC$APOE < 5], na.rm = TRUE)
mean_education_APOE_greater_5 <- mean(Merged_APOE_Edu_PC$AgeEdu[Merged_APOE_Edu_PC$APOE > 5], na.rm = TRUE)

mean_year_of_education
mean_education_APOE_less_5
mean_education_APOE_greater_5

#Mean health rate
mean_health_rate <- mean(Merged_APOE_Health_PC$GeneralHealth, na.rm = TRUE)
mean_health_rate_APOE_less_5 <- mean(Merged_APOE_Health_PC$GeneralHealth[Merged_APOE_Health_PC$APOE < 5], na.rm = TRUE)
mean_health_rate_APOE_greater_5 <- mean(Merged_APOE_Health_PC$GeneralHealth[Merged_APOE_Health_PC$APOE > 5], na.rm = TRUE)

mean_health_rate
mean_health_rate_APOE_less_5
mean_health_rate_APOE_greater_5

#Proportion of HL
total_population <- nrow(Merged_APOE_HL_PC)
total_HL <- sum(Merged_APOE_HL_PC$HL == 1, na.rm = TRUE)

# Proportion of HL in the whole dataset (fraction)
prop_HL_total <- total_HL / total_population

# For individuals with APOE less than 5
APOE_less_5_population <- sum(Merged_APOE_HL_PC$APOE < 5, na.rm = TRUE)
HL_APOE_less_5 <- sum(Merged_APOE_HL_PC$HL[Merged_APOE_HL_PC$APOE < 5] == 1, na.rm = TRUE)
prop_HL_APOE_less_5 <- HL_APOE_less_5 / APOE_less_5_population

# For individuals with APOE more than 5
APOE_greater_5_population <- sum(Merged_APOE_HL_PC$APOE > 5, na.rm = TRUE)
HL_APOE_greater_5 <- sum(Merged_APOE_HL_PC$HL[Merged_APOE_HL_PC$APOE > 5] == 1, na.rm = TRUE)
prop_HL_APOE_greater_5 <- HL_APOE_greater_5 / APOE_greater_5_population

# Print the results
prop_HL_total
prop_HL_APOE_less_5
prop_HL_APOE_greater_5

