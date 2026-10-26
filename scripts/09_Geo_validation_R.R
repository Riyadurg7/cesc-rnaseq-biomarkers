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






setwd("D:/DES/NGS/Project/analysis")                                  # project folder
table(read.csv("results/DEG_age_significant.csv", row.names = 1)$status)   # count UP and DOWN in the age list
