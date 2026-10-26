# 01_download.R  -  Download TCGA-CESC gene expression data
# Run this line by line in RStudio: put the cursor on a line and press Ctrl+Enter.
# Read the comment above each block before running it.

# ---- 0. Set the working folder (where results get saved) ----
setwd("D:/DES/NGS/Project/analysis")
options(timeout = 600)
# ---- 1. Load the package (you already installed it) ----
library(TCGAbiolinks)

# ---- 2. Ask the TCGA portal what data exists (downloads NOTHING yet) ----
query <- GDCquery(
  project = "TCGA-CESC",                        # cervical cancer
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"               # current format; ma'am's old "HTSeq - Counts" was retired
)

# How many samples matched? Expect around 300+.
length(getResults(query, cols = "cases"))

# ---- 3. Download. 20-60 min, about 2 GB. ----
# If internet drops, just run this same line again - it resumes.
GDCdownload(query, method = "api", files.per.chunk = 50)

# ---- 4. Merge all downloaded files into one object ----
data <- GDCprepare(query)

# ---- 5. First look at your data ----
dim(data)                        # genes x samples
table(data$sample_type)          # tumour vs anything else
table(data$figo_stage)           # stage groups  -> comparison 1
table(data$primary_diagnosis)    # histology     -> comparison 2
summary(data$age_at_diagnosis)   # ages (in days!) -> comparison 3

# ---- 6. Save, so you never have to download again ----
saveRDS(data, "data/CESC_data.rds")
# From tomorrow, start sessions with:
# data <- readRDS("data/CESC_data.rds")

#file.remove("data/CESC_data.rds")
#install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks")
a


setwd("D:/DES/NGS/Project/analysis")
library(TCGAbiolinks)
data <- readRDS("data/CESC_data.rds")
data <- data[, data$sample_type == "Primary Tumor"]
dim(data)   # should now say 60660 x 304

# Comparison 1: stage. Collapse substages: anything III or IV = Late, rest = Early
s <- data$figo_stage
data$stage_group <- NA
data$stage_group[grepl("III|IV", s)] <- "Late"
data$stage_group[!is.na(s) & !grepl("III|IV", s)] <- "Early"
table(data$stage_group, useNA = "ifany")

# Comparison 2: histology. Careful trap: "Adenosquamous" contains the word
# "squamous", so we must exclude it from both groups expl
icitly.
h <- data$primary_diagnosis
data$hist_group <- NA
data$hist_group[grepl("squamous", h, ignore.case = TRUE) & !grepl("Adenosquamous", h)] <- "Squamous"
data$hist_group[grepl("adenocarcinoma", h, ignore.case = TRUE) & !grepl("Adenosquamous", h)] <- "Adeno"
table(data$hist_group, useNA = "ifany")

# Comparison 3: age. Values are in DAYS; 50 years = 18263 days.
data$age_group <- ifelse(data$age_at_diagnosis >= 18263, "Older", "Younger")
table(data$age_group, useNA = "ifany")

BiocManager::install(c("EDASeq", "edgeR"))

# 4a. Remove outlier samples (correlation below 0.6 = suspicious sample)
dataPrep <- TCGAanalyze_Preprocessing(object = data, cor.cut = 0.6)
dim(dataPrep)

# 4b. Gene IDs currently have versions (ENSG00000000003.15); the normalization
# table uses them without versions - strip the ".15" part so they match.
rownames(dataPrep) <- gsub("\\.\\d+$", "", rownames(dataPrep))

# 4c. Normalize for sequencing depth and GC content (ma'am's method)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep, geneInfo = geneInfoHT, method = "gcContent")
dim(dataNorm)

# 4d. Drop the lowest 25% of genes (near-zero noise)
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut = 0.25)
dim(dataFilt)

data <- data[, data$sample_type == "Primary Tumor"]
dim(data)        # 60660 304  <- if not, stop and tell me

table(data$stage_group, useNA = "ifany")   # paste me these three
table(data$hist_group, useNA = "ifany")
table(data$age_group, useNA = "ifany")

dataPrep <- TCGAanalyze_Preprocessing(object = data, cor.cut = 0.6)
dim(dataPrep)    # second number must be 304 now

rownames(dataPrep) <- gsub("\\.\\d+$", "", rownames(dataPrep))

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep, geneInfo = geneInfoHT, method = "gcContent")
dim(dataNorm)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut = 0.25)
dim(dataFilt)


dim(dataFilt)   # should be ~45,000 x 304 - confirm the 304
groups <- data.frame(
  barcode = colnames(data),
  stage_group = data$stage_group,
  hist_group = data$hist_group,
  age_group = data$age_group
)
saveRDS(list(counts = dataFilt, groups = groups), "data/CESC_prepared.rds")
