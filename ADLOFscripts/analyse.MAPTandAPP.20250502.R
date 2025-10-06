### Step 1: Extracting LOF Variant Data for MAPT and APP
# This section loads the whole-exome sequencing LOF score files for MAPT and APP.
# Scores above 5 are considered rare LOF variants and dichotomised to 1 (carrier) or 0 (non-carrier).

# Define the genes to process
genes <- c("MAPT", "APP")
geneTables <- list()

for (gene in genes) {
    file_path <- paste0("/home/mayuan/LOFresults/UKBB.LOF.20230809.", gene, ".sco")
    raw_data <- readLines(file_path)
    split_data <- strsplit(raw_data, "\\s+")
    tab <- do.call(rbind, split_data)
    tab <- as.data.frame(tab, stringsAsFactors = FALSE)
    if (ncol(tab) == 3) {
        colnames(tab) <- c("IID", "pheno", "score")
    } else {
        next
    }
    tab$IID <- as.character(tab$IID)
    tab$pheno <- as.numeric(tab$pheno)
    tab$score <- as.numeric(tab$score)
    tab$LOF <- ifelse(tab$score > 5, 1, 0)
    geneTables[[gene]] <- tab
}

for (gene in names(geneTables)) {
  lof_count <- sum(geneTables[[gene]]$LOF == 1, na.rm = TRUE)
  cat(sprintf("Gene %s has %d LOF individuals \n", gene, lof_count))
}



### Step 2: Extracting Epilepsy Phenotype using Self-report and ICD-10 Codes
# Uses UKBB exome file, record self-report cases and ICD-10 codes-related cases to classify epilepsy cases.

# exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.txt"
# the raw UK Biobank file has a leading tab, FFS
# so add 2 to all the column numbers below
# change to using text file with only data for exome-sequenced subjects, so only need to offset by 1
##Extract ICD-10 codes reported cases
exomesFile="/SAN/ugi/UGIbiobank/data/downloaded/ukb43357.470Kexomes.txt
diagnosisFile="/home/mayuan/scripts/epilepsyDiagnosis.txt" # ICD10 codes for various epilepsy diagnoses
##specific ICD10 codes: 
##G400	G40.0 Localisation-related (focal) (partial) idiopathic epilepsy and epileptic syndromes with seizures of localised onset	33000	32990	Y
##G401	G40.1 Localisation-related (focal) (partial) symptomatic epilepsy and epileptic syndromes with simple partial seizures	33010	32990	Y
##G402	G40.2 Localisation-related (focal) (partial) symptomatic epilepsy and epileptic syndromes with complex partial seizures	33020	32990	Y
##G403	G40.3 Generalised idiopathic epilepsy and epileptic syndromes	33030	32990	Y
##G404	G40.4 Other generalised epilepsy and epileptic syndromes	33040	32990	Y
##G405	G40.5 Special epileptic syndromes	33050	32990	Y
##G406	G40.6 Grand mal seizures, unspecified (with or without petit mal)	33060	32990	Y
##G407	G40.7 Petit mal, unspecified, without grand mal seizures	33070	32990	Y
##G408	G40.8 Other epilepsy	33080	32990	Y
##G409	G40.9 Epilepsy, unspecified	33090	32990	Y
##G41	G41 Status epilepticus	33100	920	N
##G410	G41.0 Grand mal status epilepticus	33110	33100	Y
##G411	G41.1 Petit mal status epilepticus	33120	33100	Y
##G412	G41.2 Complex partial status epilepticus	33130	33100	Y
##G418	G41.8 Other status epilepticus	33140	33100	Y
##G419	G41.9 Status epilepticus, unspecified	33150	33100	Y  
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19
ICD10SourcesFile="/home/mayuan/scripts/ICD10Sources.20230811.txt" # sources of ICD10 codes - hospital admissions, causes of death
## ICD source codes
##First	Last	Source
##9542	9563	ExternalCauses
##9564	9629	DiagnosesMain
##9658	9841	DiagnosesSecondary
##10413	10625	DiagnosesHES
##9318	9319	CauseOfDeath
wd="/home/mayuan/Epilepsy.20241014"
# http://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=6
# Data Coding 6
# 1264	epilepsy
selfReportedCodes=c(1264)
selfReportFirst=5634	
selfReportLast=5769

femaleRxCols=c(5325:5340)
maleRxCols=c(5437:5448)
prescribedRxFirst=5770
prescribedRxLast=5961
setwd(wd)
ICD10Sources=data.frame(read.table(ICD10SourcesFile,header=TRUE,sep="\t"))
cmd=sprintf("tail -n +2 %s | cut -f 1",exomesFile)
for (r in 1:nrow(ICD10Sources)) {
  cmd=sprintf("%s,%d-%d",cmd,ICD10Sources$First[r]+1,ICD10Sources$Last[r]+1)
}
cmd=sprintf("%s > /home/mayuan/Epilepsy.20241014/ICD10Codes.txt",cmd)
system(cmd)

ICD10=data.frame(read.table("/home/mayuan/Epilepsy.20241014/ICD10Codes.txt",header=FALSE,sep="\t"))
epilepsy$EpilepsyICD10=0

diagnoses=data.frame(read.table(diagnosisFile,header=FALSE,sep="\t")) # ICD10 diagnoses
for (r in 1:nrow(diagnoses)) {
  epilepsy$EpilepsyICD10[rowSums(ICD10==as.character(diagnoses[r,1]), na.rm = TRUE)>0]=1
}
## Extract self-reported epilepsy diagnosis
cmd=sprintf("tail -n +2 %s | cut -f 1,%d-%d > /home/mayuan/Epilepsy.20241014/selfDiagnoses.txt",exomesFile,selfReportFirst+1,selfReportLast+1)
system(cmd)
selfReport=data.frame(read.table("/home/mayuan/Epilepsy.20241014/selfDiagnoses.txt",header=FALSE,sep="\t"))
epilepsy=data.frame(matrix(ncol=4,nrow=nrow(selfReport)))
colnames(epilepsy)=c("IID","Epilepsy","EpilepsySelf","EpilepsyICD10")
epilepsy$IID=selfReport[,1]
epilepsy$EpilepsySelf=0
for (r in 1:length(selfReportedCodes)) {
  epilepsy$EpilepsySelf[rowSums(selfReport==selfReportedCodes[r], na.rm = TRUE)>0]=1
}
##Combine both ICD10 codes and self-diagnosis data for a full epilepsy case profile 
epilepsy$Epilepsy=0
epilepsy$Epilepsy[rowSums(epilepsy[,3:4]==1, na.rm = TRUE)>0]=1
for (c in 2:ncol(epilepsy)) {
  toWrite=epilepsy[,c(1,c)]
write.table(toWrite,sprintf("/home/mayuan/Epilepsy.20241014/UKBB.%s.txt",colnames(epilepsy)[c]),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
write.table(epilepsy,"/home/mayuan/Epilepsy.20241014/UKBB.EpilepsyAll.20241014.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

phenos=as.matrix(epilepsy[,2:4])
print(cor(phenos))
### Step 3: Extracting Psychiatric Appointment Phenotype
# Extracts the field 2091 corresponding to psychiatric appointments data reported from the first visit and saves it to target phenotype directory.

wd <- "/home/mayuan/results"
setwd(wd)
dataFile <- "/SAN/ugi/UGIbiobank/data/downloaded/ukb41465.txt"
psychiatristCol <- 2091
ukbData <- read.table(dataFile, header=TRUE, sep="\t", colClasses="character")
psychiatristData <- ukbData[, c("eid", psychiatristCol)]
colnames(psychiatristData) <- c("IID", "SeenPsychiatrist")
outputFile <- "UKBB_SeenPsychiatrist.txt"
write.table(psychiatristData, file=outputFile, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
print(head(psychiatristData))
# cp UKBB_SeenPsychiatrist.txt /home/mayuan/phenos/UKBB_SP.txt



### Step 4: Extracting Educational Attainment and Overall Health Phenotype 
# Extracted using shell scripts and copied to target phenotype directory.

# bash /home/mayuan/scripts/extract.UKBB.var.41465.20230807.sh EA 435
# cp UKBB.EA.txt /home/mayuan/phenos
# bash /home/mayuan/scripts/extract.UKBB.var.41465.20230807.sh OH 642
# cp UKBB.OH.txt /home/mayuan/phenos



### Step 5: Extracting and Preprocessing Phenotype Files into a table
# Phenotype files are loaded from a specified directory and formatted for merging.

pheno_files <- list.files(path = "/home/mayuan/phenos/", full.names = TRUE)
phenoTables <- list()

for (pheno_file in pheno_files) {
  pheno_name <- gsub(".txt", "", basename(pheno_file))
  pheno_data <- read.table(pheno_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  pheno_data$IID <- as.character(pheno_data$IID)
  phenoTables[[pheno_name]] <- pheno_data
}



###Step 6: Calculate the summary statistics for each test 

cat("=== Summary Statistics by Gene and Phenotype ===\n\n")

for (gene_name in names(geneTables)) {
  gene_data <- geneTables[[gene_name]]
  gene_data$IID <- as.character(gene_data$IID)
  
  for (pheno_name in names(phenoTables)) {
    # Only process phenotypes defined in the mapping
    if (!(pheno_name %in% names(pheno_name_map))) next
    
    pheno_data <- phenoTables[[pheno_name]]
    pheno_data$IID <- as.character(pheno_data$IID)
    
    # Merge gene and phenotype data on IID (raw, unadjusted values)
    merged_data <- merge(gene_data, pheno_data, by = "IID")
    
    # Get the actual phenotype column name from your mapping
    actual_pheno_col <- pheno_name_map[[pheno_name]]
    
    # Determine if the trait is binary (logistic) or quantitative
    is_logistic <- pheno_name %in% c("UKBB.Epilepsy", "UKBB.SP")
    
    if (is_logistic) {
      # Convert to binary: any value > 0 is a case (1), otherwise control (0)
      merged_data[[actual_pheno_col]] <- ifelse(merged_data[[actual_pheno_col]] > 0, 1, 0)
      
      # Separate carriers and non-carriers
      carriers <- merged_data[merged_data$LOF == 1, ]
      noncarriers <- merged_data[merged_data$LOF == 0, ]
      
      # Compute % cases (unadjusted)
      pct_cases_carriers <- 100 * sum(carriers[[actual_pheno_col]] == 1, na.rm = TRUE) / nrow(carriers)
      pct_cases_noncarriers <- 100 * sum(noncarriers[[actual_pheno_col]] == 1, na.rm = TRUE) / nrow(noncarriers)
      
      cat(sprintf("Gene: %s, Phenotype: %s\n", gene_name, pheno_name))
      cat(sprintf("  %% Cases in Carriers: %.2f%%\n", pct_cases_carriers))
      cat(sprintf("  %% Cases in Non-Carriers: %.2f%%\n\n", pct_cases_noncarriers))
      
    } else {
      # For quantitative traits: compute mean and SD
      carriers <- merged_data[merged_data$LOF == 1, ]
      noncarriers <- merged_data[merged_data$LOF == 0, ]
      
      mean_sd_carriers <- sprintf("%.2f (%.2f)", 
                                  mean(carriers[[actual_pheno_col]], na.rm = TRUE),
                                  sd(carriers[[actual_pheno_col]], na.rm = TRUE))
      mean_sd_noncarriers <- sprintf("%.2f (%.2f)", 
                                     mean(noncarriers[[actual_pheno_col]], na.rm = TRUE),
                                     sd(noncarriers[[actual_pheno_col]], na.rm = TRUE))
      
      cat(sprintf("Gene: %s, Phenotype: %s\n", gene_name, pheno_name))
      cat(sprintf("  Mean (SD) in Carriers: %s\n", mean_sd_carriers))
      cat(sprintf("  Mean (SD) in Non-Carriers: %s\n\n", mean_sd_noncarriers))
    }
  }



### Step 7: Merging Sex and Principal Component Data into Covariates Table
# Loads and cleans PCA and sex data into a single covariates table.

PCFile <- "/home/rejudcu/UKBB/RAPfiles/covars/covars/ukb23158.common.all.20230806.eigenvec.txt"
PCTable <- read.table(PCFile, header = FALSE)
colnames(PCTable) <- c("FID", "IID", paste0("PC", 1:10))

SexFile <- "/home/mayuan/phenos/UKBB.ST.txt"
SexTable <- read.table(SexFile, header = TRUE, sep = "\t")
colnames(SexTable)[colnames(SexTable) == "IID.sex.20230807"] <- "IID"
SexTable$IID <- trimws(as.character(SexTable$IID))
PCTable$IID <- trimws(as.character(PCTable$IID))

CovarsTable <- merge(SexTable, PCTable, by = "IID")
colnames(CovarsTable) <- gsub("^sex\\.[0-9]+", "sex", colnames(CovarsTable))
na_indices <- which(grepl("^NA", colnames(CovarsTable)))
if (length(na_indices) > 0) {
    colnames(CovarsTable)[na_indices] <- paste0("PC", seq(11, 10 + length(na_indices)))
}
print(colnames(CovarsTable))
print(nrow(CovarsTable))
print(head(CovarsTable))
write.table(CovarsTable, file = "/home/mayuan/scripts/CovarsTable_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)



### Step 8: Performing Regression Analyses for LOF and Phenotype Associations
# Merges all data and runs multiple regression models per phenotype.
# Logistic regression is used for binary phenotypes, linear for quantitative ones.

# LOAD COVARIATE DATA
covars_file <- "/home/mayuan/scripts/CovarsTable_output.txt"
CovarsTable <- read.table(covars_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# LRT helper function
LRT_by_hand <- function(full_model, null_model) {
  ll_full <- as.numeric(logLik(full_model))
  ll_null <- as.numeric(logLik(null_model))
  test_stat <- 2 * (ll_full - ll_null)
  p_val <- pchisq(test_stat, df = 1, lower.tail = FALSE)
  return(list(test_stat = test_stat, p_val = p_val))
}

pcs <- paste(paste0("PC", 1:10), collapse = " + ")
pheno_name_map <- list("UKBB.EA" = "EA", "UKBB.OH" = "OH", "UKBB.Epilepsy" = "Epilepsy", "UKBB.SP" = "SP")
results <- list()

for (gene_name in names(geneTables)) {
  gene_data <- geneTables[[gene_name]]
  gene_data$IID <- as.character(gene_data$IID)
  
  for (pheno_name in names(phenoTables)) {
    if (!(pheno_name %in% names(pheno_name_map))) next
    pheno_data <- phenoTables[[pheno_name]]
    pheno_data$IID <- as.character(pheno_data$IID)
    test_data <- merge(gene_data, pheno_data, by = "IID")
    test_data <- merge(test_data, CovarsTable, by = "IID")
    actual_pheno_col <- pheno_name_map[[pheno_name]]
    full_formula_str <- paste(actual_pheno_col, "~ LOF + sex +", pcs)
    null_formula_str <- paste(actual_pheno_col, "~ sex +", pcs)
    full_formula <- as.formula(full_formula_str)
    null_formula <- as.formula(null_formula_str)
    is_logistic <- pheno_name %in% c("UKBB.Epilepsy", "UKBB.SP")
    if (is_logistic) {
      test_data[[actual_pheno_col]] <- ifelse(test_data[[actual_pheno_col]] > 0, 1, 0)
    }
    if (is_logistic) {
      model_full <- glm(full_formula, data = test_data, family = binomial)
      model_null <- glm(null_formula, data = test_data, family = binomial)
    } else {
      model_full <- lm(full_formula, data = test_data)
      model_null <- lm(null_formula, data = test_data)
    }
    lrt <- LRT_by_hand(model_full, model_null)
    b <- coef(model_full)["LOF"]
    se <- sqrt(vcov(model_full)["LOF", "LOF"])
    if (is_logistic) {
      OR <- exp(b)
      OR_low <- exp(b - 2 * se)
      OR_high <- exp(b + 2 * se)
      estimate_str <- sprintf("OR = %.2f (95%% CI: %.2f - %.2f)", OR, OR_low, OR_high)
    } else {
      beta_low <- b - 2 * se
      beta_high <- b + 2 * se
      estimate_str <- sprintf("Beta = %.2f (95%% CI: %.2f - %.2f)", b, beta_low, beta_high)
    }
    model_summary <- capture.output(summary(model_full))
    result_key <- paste(gene_name, pheno_name, sep = "_")
    results[[result_key]] <- list(
      summary = model_summary,
      LRT_stat = sprintf("LRT statistic: %.4f", lrt$test_stat),
      LRT_p_val = sprintf("LRT p-value: %.4f", lrt$p_val),
      effect = estimate_str
    )
  }
}

output_file <- "/home/mayuan/results/multiple_regression_results.txt"
sink(output_file)
for (key in names(results)) {
  cat("\n===== Regression Results for", key, "=====\n")
  cat(paste(results[[key]]$summary, collapse = "\n"), "\n")
  cat(results[[key]]$LRT_stat, "\n")
  cat(results[[key]]$LRT_p_val, "\n")
  cat(results[[key]]$effect, "\n\n")
}
sink()
cat("Results saved to:", output_file, "\n")


