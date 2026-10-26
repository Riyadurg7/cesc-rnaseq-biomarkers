# Load the metrics table; check.names=FALSE keeps column names with spaces intact
m <- read.csv("results/stage_node_metrics.csv", check.names = FALSE)

# See the exact column names R found (so we use the right ones)
colnames(m)

# Top 20 genes by Degree
head(m[order(-m$Degree), "display name"], 20)

# Top 20 genes by Betweenness
head(m[order(-m$BetweennessCentrality), "display name"], 20)



# Recompute both top-20 lists and keep their overlap = stage hub genes
deg20 <- head(m[order(-m$Degree), "display name"], 20)
bet20 <- head(m[order(-m$BetweennessCentrality), "display name"], 20)
stage_hubs <- intersect(deg20, bet20)
stage_hubs                                   # should show the 12 names
writeLines(stage_hubs, "results/hubs_stage.txt")


#Age_Node_Metrics

# Age hubs = overlap of top-20 Degree and top-20 Betweenness
m <- read.csv("results/age_node_metrics.csv", check.names = FALSE)
deg20 <- head(m[order(-m$Degree), "display name"], 20)
bet20 <- head(m[order(-m$BetweennessCentrality), "display name"], 20)
age_hubs <- intersect(deg20, bet20)
age_hubs
writeLines(age_hubs, "results/hubs_age.txt")


#For Histology

# Load histology significant genes, keep the 1000 with smallest FDR, export symbols
deg_hist <- read.csv("results/DEG_histology_significant.csv", row.names = 1)
top_hist <- head(deg_hist[order(deg_hist$FDR), ], 1000)
writeLines(unique(na.omit(top_hist$gene_name)), "results/genes_histology_top1000.txt")



# Histology hubs = overlap of top-20 Degree and top-20 Betweenness
m <- read.csv("results/hist_node_metrics.csv", check.names = FALSE)
deg20 <- head(m[order(-m$Degree), "display name"], 20)
bet20 <- head(m[order(-m$BetweennessCentrality), "display name"], 20)
hist_hubs <- intersect(deg20, bet20)
hist_hubs
writeLines(hist_hubs, "results/hubs_histology.txt")
