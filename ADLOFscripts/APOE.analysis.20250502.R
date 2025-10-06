# Load required datasets
apoe_data <- read.table("UKBB.LOF.20230809.APOE.sco", header = FALSE, sep = "", fill = TRUE, stringsAsFactors = FALSE)[,-2]
colnames(apoe_data) <- c("IID", "APOE")
sex_data <- read.table("UKBB.Sex.txt", header = TRUE, sep = "\t")
PC_data <- read.table("ukb23158.common.all.20230806.eigenvec.txt", header = TRUE, sep = "\t")[, -1]

# Remove NA
apoe_data <- na.omit(apoe_data)
sex_data <- na.omit(sex_data)

# Merge APOE, Sex and PC data
merged_APOE_Sex <- merge(sex_data, apoe_data, by = "IID", all = FALSE)

# List of datasets with their file names and output names for files that include PCs
data_files <- list("UKBB.HL.20230915.txt" = "APOE.HL.PC.20241028.txt",
                   "UKBB.GeneralHealth.txt" = "APOE.Health.PC.20241028.txt",
                   "UKBB.AgeEdu.txt" = "APOE.Education.PC.20241028.txt")

# Loop through each data file, clean, merge, and then merge with PCs before saving
for (file_name in names(data_files)) {
  data <- read.table(file_name, header = TRUE, sep = "\t")
  data_clean <- na.omit(data)
  data_clean <- data_clean[rowSums(data_clean < 0) == 0, ]
  
  merged_data <- merge(merged_APOE_Sex, data_clean, by = "IID", all = FALSE)
  
  final_data_PC <- merge(merged_data, PC_data, by = "IID", all = FALSE)
  
  write.table(final_data_PC, file = data_files[[file_name]], sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("Final merged files including PCs saved.\n")

#------------Linear Regression----------
files <- list(
  AgeEdu = "APOE.Education.PC.20241028.txt",
  GeneralHealth = "APOE.Health.PC.20241028.txt",
  HL = "APOE.HL.PC.20241028.txt"
)

results <- list()
results_P_value <- list ()
results_chi <- list()

for (var_name in names(files)) {
  # Read the corresponding file
  data <- read.table(files[[var_name]], header = TRUE, sep = "\t")
  
  if (var_name == "HL") {
    # Logistic regression for HL
    N0 <- glm(HL ~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
              data = data, family = binomial)
    N1 <- glm(HL ~ APOE + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, 
              data = data, family = binomial)
  } else {
    # Linear regression for continuous outcomes 
    N0 <- lm(as.formula(paste(var_name, "~ Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                                          PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), 
             data = data)
    N1 <- lm(as.formula(paste(var_name, "~ APOE + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + 
                                          PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20")), 
             data = data)
  }
  
  # Calculate log-likelihoods and chi-square test
  logLik_N0 <- logLik(N0)
  logLik_N1 <- logLik(N1)
  chi_sq_stat <- 2 * (logLik_N1 - logLik_N0)
  p_value <- pchisq(chi_sq_stat, df = 1, lower.tail = FALSE)
  
  coef_N0 <- summary(N0)$coefficients
  coef_N1 <- summary(N1)$coefficients
  
  results[[var_name]] <- list(
    Model_N0_Coefficients = coef_N0,
    Model_N1_Coefficients = coef_N1,
    Chi_Square_Statistic = chi_sq_stat,
    P_Value = p_value
  )
  
  results_chi[[var_name]] <- list(
    Chi_Square_Statistic = chi_sq_stat
  )
  results_P_value[[var_name]] <- list(
    P_Value = p_value
  )
  
  # Extract coefficients from N0 and N1 models and save them in a readable format
  coefficients_N0 <- as.data.frame(summary(N0)$coefficients)
  coefficients_N0$Model <- "N0 (Null Model)"
  
  coefficients_N1 <- as.data.frame(summary(N1)$coefficients)
  coefficients_N1$Model <- "N1 (Alternative Model)"
  
  coefficients_table <- rbind(coefficients_N0, coefficients_N1)
  
  file_name <- paste0("coefficients.", var_name, ".20241030.txt")
  write.table(coefficients_table, file = file_name, sep = "\t", row.names = TRUE, quote = FALSE)
  }
  

# Print results
print(results_P_value)
