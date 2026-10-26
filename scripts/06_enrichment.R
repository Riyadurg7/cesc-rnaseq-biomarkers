# ---- Block 1: install the enrichment packages (one-time only) ----
# clusterProfiler = the tool that does GO and KEGG enrichment
# org.Hs.eg.db    = the human gene database (maps gene symbols to IDs)
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

library(clusterProfiler)   # if this loads without "Error", it's installed
library(org.Hs.eg.db)      # same check for the gene database
BiocManager::install("GO.db")   # install the missing helper (the GO terms database)
library(clusterProfiler)   # should load now


# ---- Block 2: load the packages ----
setwd("D:/DES/NGS/Project/analysis")    # make sure R is working from the analysis folder
library(clusterProfiler)                # the enrichment tool (GO + KEGG tests)
library(org.Hs.eg.db)                   # human gene database (symbol -> ID translation)

# ---- Block 3: read the gene lists we made in 04_overlap.R ----
genes_stage <- readLines("results/genes_stage.txt")      # stage DEG symbols (should be ~453)
genes_hist  <- readLines("results/genes_histology.txt")  # histology DEG symbols (~3415)
genes_age   <- readLines("results/genes_age.txt")        # age DEG symbols (~772)
genes_core  <- readLines("results/genes_core111.txt")    # the 111 genes shared by all three

# quick count check - do the numbers match what we expect?
length(genes_stage); length(genes_hist); length(genes_age); length(genes_core)



# ---- Block 4: translate gene symbols to Entrez IDs ----
# bitr = "biological ID translator" - looks up each symbol in the human database
ids_stage <- bitr(genes_stage, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ids_hist  <- bitr(genes_hist,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ids_age   <- bitr(genes_age,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ids_core  <- bitr(genes_core,  fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# how many survived translation? (each ids_ object is a 2-column table: SYMBOL, ENTREZID)
nrow(ids_stage); nrow(ids_hist); nrow(ids_age); nrow(ids_core)



# ---- Block 5: GO enrichment for the stage genes ----
ego_stage <- enrichGO(
  gene          = ids_stage$ENTREZID,   # our stage genes, as Entrez IDs
  OrgDb         = org.Hs.eg.db,         # the human database
  ont           = "BP",                 # "BP" = Biological Process (what the genes DO)
  pAdjustMethod = "BH",                 # multiple-testing correction (same FDR idea as in DEA)
  pvalueCutoff  = 0.05,                 # keep processes significant at 5%
  readable      = TRUE                  # show gene SYMBOLS in results, not numbers
)

# how many significant processes, and what are the top 10?
nrow(ego_stage)                                                        # count of enriched GO terms
head(as.data.frame(ego_stage)[, c("Description", "Count", "p.adjust")], 10)  # top 10: name, gene count, FDR


# ---- Block 5b: save stage GO results + dotplot ----
write.csv(as.data.frame(ego_stage), "results/GO_stage.csv")    # full table -> results folder

library(ggplot2)                                # plotting (for ggsave)
p <- dotplot(ego_stage, showCategory = 15)      # dotplot: dot size = no. of genes, colour = FDR
ggsave("results/GO_stage_dotplot.png", p, width = 9, height = 7, dpi = 300)  # save as image
p                                               # show it in the Plots pane too





# ---- Block 6: GO enrichment for histology, age, and core-111 ----
ego_hist <- enrichGO(gene = ids_hist$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)   # histology genes
ego_age  <- enrichGO(gene = ids_age$ENTREZID,  OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)   # age genes
ego_core <- enrichGO(gene = ids_core$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)   # the shared 111

nrow(ego_hist); nrow(ego_age); nrow(ego_core)   # how many processes each found

write.csv(as.data.frame(ego_hist), "results/GO_histology.csv")   # save all three tables
write.csv(as.data.frame(ego_age),  "results/GO_age.csv")
write.csv(as.data.frame(ego_core), "results/GO_core111.csv")

ggsave("results/KEGG_stage_dotplot.png",     dotplot(kk_stage, showCategory = 15), width = 9, height = 7, dpi = 300)  # save plots
ggsave("results/KEGG_histology_dotplot.png", dotplot(kk_hist,  showCategory = 15), width = 9, height = 7, dpi = 300)
ggsave("results/KEGG_age_dotplot.png",       dotplot(kk_age,   showCategory = 15), width = 9, height = 7, dpi = 300)
ggsave("results/KEGG_core111_dotplot.png",   dotplot(kk_core,  showCategory = 15), width = 9, height = 7, dpi = 300)




# ---- Block 8: master candidate table ----
# One row per hub gene; one column per piece of evidence; score = how many ticks

hubs_stage <- readLines("results/hubs_stage.txt")        # the 12 stage hubs
hubs_age   <- readLines("results/hubs_age.txt")          # the 12 age hubs
hubs_hist  <- readLines("results/hubs_histology.txt")    # the 8 histology hubs

all_hubs <- sort(unique(c(hubs_stage, hubs_age, hubs_hist)))   # every hub gene, listed once

master <- data.frame(
  gene       = all_hubs,
  hub_stage  = all_hubs %in% hubs_stage,    # TRUE if it's a stage hub
  hub_age    = all_hubs %in% hubs_age,      # TRUE if it's an age hub
  hub_hist   = all_hubs %in% hubs_hist,     # TRUE if it's a histology hub
  in_core111 = all_hubs %in% genes_core     # TRUE if it's also in the shared 111
)
master$score <- rowSums(master[, 2:5])      # count the TRUEs (1 to 4)
master <- master[order(-master$score, master$gene), ]   # highest score on top

master                                      # print the scoreboard
write.csv(master, "results/candidate_master_table.csv", row.names = FALSE)   # save it





# ---- fix: recreate the missing stage hubs file ----
m <- read.csv("results/stage_node_metrics.csv", check.names = FALSE)     # the Cytoscape export
deg20 <- head(m[order(-m$Degree), "display name"], 20)                   # top 20 by Degree
bet20 <- head(m[order(-m$BetweennessCentrality), "display name"], 20)    # top 20 by Betweenness
hubs_stage <- intersect(deg20, bet20)                                    # in BOTH lists = hub
hubs_stage                                                               # should print your 12 (SST, GRIA2, KRT20...)
writeLines(hubs_stage, "results/hubs_stage.txt")                         # save the file this time