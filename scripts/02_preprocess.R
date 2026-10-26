# ============================================================
# 02_preprocess.R
# Turns the raw downloaded data (CESC_data.rds) into a clean,
# normalized, filtered count table ready for DE analysis.
#
# NOTE: You already ran all of this in the console on day 1.
# The output file data/CESC_prepared.rds already exists, so you
# do NOT need to run this again. This file is the written record
# of HOW that file was made (needed for your report + GitHub).
# ============================================================

# ---- Block 1: setup ----
setwd("D:/DES/NGS/Project/analysis")   # work from the analysis folder (R paths use / not \)
library(TCGAbiolinks)                  # TCGA functions: Preprocessing, Normalization, Filtering
library(EDASeq)                        # engine behind TCGAanalyze_Normalization
library(edgeR)                         # engine behind filtering and DEA

# ---- Block 2: load data and keep only primary tumours ----
data <- readRDS("data/CESC_data.rds")  # load the object saved by 01_download.R

data <- data[, data$sample_type == "Primary Tumor"]  # drop normal/metastatic samples (309 -> 304)
dim(data)                              # check: should print 60660 304

# ---- Block 3: build the three patient groups ----
# Group 1: FIGO stage, Early (I-II) vs Late (III-IV)
s <- data$figo_stage                   # stage text for each patient, e.g. "Stage IIIB"
stage_group <- ifelse(grepl("III|IV", s), "Late", "Early")  # any III or IV in the text = Late
stage_group[is.na(s)] <- NA            # patients with no recorded stage stay NA (left out of DEA)
table(stage_group)                     # check: Early 196, Late 54

# Group 2: histology, Squamous vs Adeno
h <- data$primary_diagnosis            # diagnosis text for each patient
hist_group <- rep(NA, length(h))       # start with all NA
hist_group[grepl("squamous", h, ignore.case = TRUE) &
           !grepl("adenosquamous", h, ignore.case = TRUE)] <- "Squamous"  # squamous, but not the mixed type
hist_group[grepl("adenocarcinoma", h, ignore.case = TRUE)] <- "Adeno"     # adenocarcinoma types
table(hist_group)                      # check: Adeno 50, Squamous 246

# Group 3: age at diagnosis, under 50 vs 50 and over
# age_at_diagnosis is stored in DAYS, so 50 years = 50 * 365.25 = 18263 days
age_group <- ifelse(data$age_at_diagnosis >= 18263, "Older", "Younger")
table(age_group)                       # check: the two groups we used were 178 and 121

# ---- Block 4: remove outlier samples ----
dataPrep <- TCGAanalyze_Preprocessing(object = data, cor.cut = 0.6)
# correlates every sample with every other; a sample that looks unlike
# all the rest (correlation < 0.6) is treated as a bad library and dropped

# strip version numbers off gene IDs: ENSG00000141510.16 -> ENSG00000141510
# needed so our IDs match the geneInfoHT table used in the next step
rownames(dataPrep) <- gsub("\\.\\d+$", "", rownames(dataPrep))

# ---- Block 5: normalization ----
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfoHT,   # gene lengths + GC content table
                                      method = "gcContent")    # correct GC bias + library size
# after this, expression values are comparable between samples

# ---- Block 6: filter low-expression genes ----
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut = 0.25)   # drop the bottom 25% lowest-expressed genes (mostly noise)
dim(dataFilt)                          # check: should print 45169 304

# ---- Block 7: save everything 03_dea.R needs ----
groups <- data.frame(
  barcode     = colnames(data),        # sample barcode for each patient
  stage_group = stage_group,           # Early / Late / NA
  hist_group  = hist_group,            # Squamous / Adeno / NA
  age_group   = age_group              # Younger / Older / NA
)

saveRDS(list(counts = dataFilt, groups = groups),  # one bundle: clean counts + group labels
        "data/CESC_prepared.rds")                  # this is the file 03_dea.R reads
