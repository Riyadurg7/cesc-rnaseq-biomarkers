# ---- 08_deseq2.R : re-check our DEGs with a second method ----

# ---- Block 1: install DESeq2 (one-time; BIG install, 10-20 min) ----
BiocManager::install("DESeq2")


# ---- run this first, and PASTE IT AT THE TOP of 08_deseq2.R too ----
setwd("D:/DES/NGS/Project/analysis")   # point R back at the project folder
getwd()                                # confirm - should print the analysis path

library(DESeq2)                        # reload the tools (fresh session = empty memory)
library(SummarizedExperiment)



# ---- Block 2: load RAW counts + group labels ----
library(DESeq2)                     # the second calculator
library(SummarizedExperiment)       # helper to pull the count table out

data <- readRDS("data/CESC_data.rds")                    # the ORIGINAL raw download (script 01)
data <- data[, data$sample_type == "Primary Tumor"]      # keep the 304 tumours, as always

raw <- assay(data, "unstranded")                         # raw counts: genes x 304, untouched numbers
rownames(raw) <- gsub("\\.\\d+$", "", rownames(raw))     # strip gene id versions (same as script 02)

prep   <- readRDS("data/CESC_prepared.rds")              # our saved group labels live here
groups <- prep$groups

idx2  <- match(colnames(raw), groups$barcode)            # pair labels by NAME, never position (the lesson!)
stage <- groups$stage_group[idx2]                        # Early/Late for each raw column

dim(raw)             # check: 60660 304
table(stage)         # check: Early 196, Late 54





# ---- Block 3: DESeq2 on the stage comparison ----
keep <- !is.na(stage)                     # drop the 54 patients with no stage recorded

dds <- DESeqDataSetFromMatrix(
  countData = raw[, keep],                # raw counts, labeled patients only
  colData   = data.frame(stage = factor(stage[keep], levels = c("Early", "Late"))),  # the labels; Early = baseline
  design    = ~ stage                     # the question: do counts depend on stage?
)

dds <- dds[rowSums(counts(dds)) >= 10, ]  # drop genes with almost zero reads everywhere (pure noise)

dds <- DESeq(dds)                         # THE big computation - normalize + fit + test every gene

res_stage <- results(dds, contrast = c("stage", "Late", "Early"))  # positive logFC = up in LATE (same direction rule as always)
summary(res_stage)                        # headline numbers: how many genes up / down


# ---- Block 4: compare DESeq2 vs edgeR on the stage contrast ----
deseq_sig <- rownames(res_stage)[which(res_stage$padj < 0.05 &
                                         abs(res_stage$log2FoldChange) > 2)]  # SAME cutoffs as our edgeR run
length(deseq_sig)                                 # how many genes DESeq2 calls significant

edger     <- read.csv("results/DEG_stage_significant.csv", row.names = 1)  # our edgeR list (the 453)
edger_sig <- rownames(edger)
length(edger_sig)                                 # should say 453

both <- intersect(deseq_sig, edger_sig)           # genes BOTH methods agree on
length(both)                                      # the confirmed core
round(length(both) / length(edger_sig) * 100)     # % of our edgeR DEGs that survived the cross-check

writeLines(both, "results/stage_DEGs_confirmed_both.txt")   # save the confirmed list

# full-circle moment: did DESeq2 believe ALDH8A1, the one-patient wonder?
"ENSG00000164816" %in% deseq_sig                  # TRUE = believed, FALSE = rejected as outlier-driven




# ---- Block 5: is the disagreement biology or thresholds? ----
sub <- res_stage[rownames(res_stage) %in% edger_sig, ]   # DESeq2's results for OUR 453 genes
nrow(sub)                                                # how many of the 453 DESeq2 even tested

same_dir <- sign(sub$log2FoldChange) == sign(edger[rownames(sub), "logFC"])  # same up/down direction?
sum(same_dir, na.rm = TRUE)                              # directional agreement count

sum(sub$padj < 0.05, na.rm = TRUE)                       # significant in DESeq2 (ignoring logFC size)

# and the question we actually care about: do our 12 stage HUBS survive?
hubs          <- readLines("results/hubs_stage.txt")               # SST, GRIA2, KRT20...
hub_rows      <- edger[edger$gene_name %in% hubs, ]                # their ENSG ids
hub_rows$gene_name[rownames(hub_rows) %in% deseq_sig]              # hubs confirmed at FULL strictness
hub_rows$gene_name[rownames(hub_rows) %in% rownames(sub)[which(sub$padj < 0.05)]]  # hubs significant in DESeq2 at all



# ---- Block 6: DESeq2 on the histology contrast ----
hist_g <- groups$hist_group[idx2]                        # Squamous/Adeno labels, name-matched
keep   <- !is.na(hist_g)                                 # drop unlabeled patients

dds_h <- DESeqDataSetFromMatrix(
  countData = raw[, keep],
  colData   = data.frame(hist = factor(hist_g[keep], levels = c("Squamous", "Adeno"))),  # Squamous = baseline
  design    = ~ hist
)
dds_h <- dds_h[rowSums(counts(dds_h)) >= 10, ]           # drop no-signal genes
dds_h <- DESeq(dds_h)                                    # long stop sign again (296 patients)

res_hist <- results(dds_h, contrast = c("hist", "Adeno", "Squamous"))  # positive logFC = up in ADENO (same rule as 03)

deseq_hist <- rownames(res_hist)[which(res_hist$padj < 0.05 &
                                         abs(res_hist$log2FoldChange) > 2)]   # our strict cutoffs
edger_h    <- read.csv("results/DEG_histology_significant.csv", row.names = 1)  # edgeR's 3,415

length(deseq_hist); nrow(edger_h)                        # each method's count
length(intersect(deseq_hist, rownames(edger_h)))         # agreed by both
writeLines(intersect(deseq_hist, rownames(edger_h)),
           "results/histology_DEGs_confirmed_both.txt")  # save confirmed list

res_hist["ENSG00000146648", c("log2FoldChange", "padj")] # EGFR's verdict from the second doctor




# ---- Block 7: DESeq2 on the age contrast ----
age_g <- groups$age_group[idx2]                          # Younger/Older labels, name-matched
keep  <- !is.na(age_g)                                   # drop unlabeled

dds_a <- DESeqDataSetFromMatrix(
  countData = raw[, keep],
  colData   = data.frame(age = factor(age_g[keep], levels = c("Younger", "Older"))),  # Younger = baseline
  design    = ~ age
)
dds_a <- dds_a[rowSums(counts(dds_a)) >= 10, ]           # drop no-signal genes
dds_a <- DESeq(dds_a)                                    # last long stop sign of 08

res_age <- results(dds_a, contrast = c("age", "Older", "Younger"))   # positive logFC = up in OLDER

deseq_age <- rownames(res_age)[which(res_age$padj < 0.05 &
                                       abs(res_age$log2FoldChange) > 2)]
edger_a   <- read.csv("results/DEG_age_significant.csv", row.names = 1)   # edgeR's 772

length(deseq_age); nrow(edger_a)                         # each method's count
length(intersect(deseq_age, rownames(edger_a)))          # agreed by both
writeLines(intersect(deseq_age, rownames(edger_a)),
           "results/age_DEGs_confirmed_both.txt")        # save confirmed list




# ---- Block 8: save tables + final hub check ----
write.csv(as.data.frame(res_stage), "results/DESeq2_stage_full.csv")      # full results, all 3 contrasts
write.csv(as.data.frame(res_hist),  "results/DESeq2_histology_full.csv")
write.csv(as.data.frame(res_age),   "results/DESeq2_age_full.csv")

hubs_a <- readLines("results/hubs_age.txt")                # the 12 age hubs
rows_a <- edger_a[edger_a$gene_name %in% hubs_a, ]         # their ENSG ids
rows_a$gene_name[rownames(rows_a) %in% deseq_age]    # which survive at full strictness






# ---- 09_geo_validation.R : check our genes in NEW patients (not TCGA) ----

# ---- Block 1: setup ----
setwd("D:/DES/NGS/Project/analysis")    # project folder, as always
BiocManager::install("GEOquery")        # downloader for GEO datasets (one-time; type n if asked)
library(GEOquery)                       # load it

# ---- Block 2: download GSE63514 ----
options(timeout = 600)                  # allow slow downloads (same trick as day 1)
gse  <- getGEO("GSE63514", GSEMatrix = TRUE)   # download the processed data
eset <- gse[[1]]                               # unpack it
dim(exprs(eset))                               # check: probes x samples, expect ~54675 x 128


gse  <- getGEO("GSE63514", GSEMatrix = TRUE, getGPL = FALSE)   # retry, without the huge annotation file
eset <- gse[[1]]
dim(exprs(eset))                                               # expect ~54675 x 128






# ---- Block 3: find the patient groups ----
ph <- pData(eset)        # the description table: one row per sample
colnames(ph)             # print all column names - we're hunting for one about disease/group



table(ph$`tissue type:ch1`)    # count the samples in each group






# ---- Block 4: short group labels in progression order ----
g <- ph$`tissue type:ch1`                      # the long names
group_geo <- rep(NA, length(g))                # start empty
group_geo[grepl("normal", g)]     <- "Normal"  # healthy cervix
group_geo[grepl("low grade", g)]  <- "CIN1"    # mild pre-cancer
group_geo[grepl("moderate", g)]   <- "CIN2"    # middle pre-cancer
group_geo[grepl("high grade", g)] <- "CIN3"    # severe pre-cancer
group_geo[grepl("cancer", g)]     <- "Cancer"  # the real thing

group_geo <- factor(group_geo, levels = c("Normal", "CIN1", "CIN2", "CIN3", "Cancer"))  # fix the order
table(group_geo)                               # check: 24 14 22 40 28



# ---- Block 5: dictionary - which probes measure our genes? ----
BiocManager::install("hgu133plus2.db")   # dictionary for this exact microarray (type n if asked)
library(hgu133plus2.db)

probes <- AnnotationDbi::select(hgu133plus2.db,
                                keys    = c("EGFR", "KRT20", "CDH17"),  # our three genes
                                keytype = "SYMBOL",
                                columns = "PROBEID")                    # give me their probe ids
probes                                   # one gene can have several probes - that's normal




library(hgu133plus2.db)                                    # load the dictionary

probes <- AnnotationDbi::select(hgu133plus2.db,
                                keys    = c("EGFR", "KRT20", "CDH17"),
                                keytype = "SYMBOL",
                                columns = "PROBEID")
probes                                                     # the gene -> probe list




# ---- Block 6: one expression row per gene ----
ex <- exprs(eset)                                    # probe x sample table (already log scale)

pick_probe <- function(ids) {                        # recipe: from a gene's probes, keep the strongest
  ids <- ids[ids %in% rownames(ex)]                  # only probes present in this dataset
  ids[which.max(rowMeans(ex[ids, , drop = FALSE]))]  # the one with highest average signal
}

p_egfr  <- pick_probe(probes$PROBEID[probes$SYMBOL == "EGFR"])
p_krt20 <- pick_probe(probes$PROBEID[probes$SYMBOL == "KRT20"])
p_cdh17 <- pick_probe(probes$PROBEID[probes$SYMBOL == "CDH17"])
p_egfr; p_krt20; p_cdh17                             # which probes won




# ---- Block 7: the payoff - our genes across Normal -> CIN1 -> CIN2 -> CIN3 -> Cancer ----
png("results/GEO_GSE63514_boxplots.png", width = 12, height = 5, units = "in", res = 300)  # open image file
par(mfrow = c(1, 3))                                 # three plots side by side
boxplot(ex[p_egfr, ]  ~ group_geo, main = "EGFR",  xlab = "", ylab = "expression (log)", col = "lightblue")
boxplot(ex[p_krt20, ] ~ group_geo, main = "KRT20", xlab = "", ylab = "", col = "lightpink")
boxplot(ex[p_cdh17, ] ~ group_geo, main = "CDH17", xlab = "", ylab = "", col = "lightgreen")
dev.off()                                            # save the file

par(mfrow = c(1, 3))                                 # draw again, this time in the Plots pane
boxplot(ex[p_egfr, ]  ~ group_geo, main = "EGFR",  xlab = "", ylab = "expression (log)", col = "lightblue")
boxplot(ex[p_krt20, ] ~ group_geo, main = "KRT20", xlab = "", ylab = "", col = "lightpink")
boxplot(ex[p_cdh17, ] ~ group_geo, main = "CDH17", xlab = "", ylab = "", col = "lightgreen")
par(mfrow = c(1, 1))                                 # back to one plot at a time

# does each gene differ across the 5 groups? (simple test, no assumptions)
kruskal.test(ex[p_egfr, ]  ~ group_geo)$p.value      # EGFR
kruskal.test(ex[p_krt20, ] ~ group_geo)$p.value      # KRT20
kruskal.test(ex[p_cdh17, ] ~ group_geo)$p.value      # CDH17
