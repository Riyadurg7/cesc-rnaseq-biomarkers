#Step 4

# ---- 04_overlap.R : which genes appear in more than one comparison? ----

# Work from the project folder
setwd("D:/DES/NGS/Project/analysis")

# Load the significant-gene tables that script 03 saved
# (row.names = 1 means: first column of the CSV = gene IDs)
deg_stage <- read.csv("results/DEG_stage_significant.csv", row.names = 1)
deg_hist  <- read.csv("results/DEG_histology_significant.csv", row.names = 1)
deg_age   <- read.csv("results/DEG_age_significant.csv", row.names = 1)

# Count the genes in each list (should match your three status tables)
nrow(deg_stage); nrow(deg_hist); nrow(deg_age)

# Take just the gene IDs from each list
stage_genes <- rownames(deg_stage)
hist_genes  <- rownames(deg_hist)
age_genes   <- rownames(deg_age)

# How many genes are shared between each PAIR of lists?
# intersect() = "give me what is in both"
length(intersect(stage_genes, hist_genes))
length(intersect(stage_genes, age_genes))
length(intersect(hist_genes, age_genes))

# Genes present in ALL THREE lists - the core candidates
core <- Reduce(intersect, list(stage_genes, hist_genes, age_genes))
length(core)

# Look at the core genes by name, with their stage-contrast stats
deg_stage[core, c("gene_name", "logFC", "FDR")]



# One-time install of the Venn diagram package
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

# Draw the three-circle Venn: each circle = one comparison's DEG list
p <- ggVennDiagram(list(Stage = stage_genes, Histology = hist_genes, Age = age_genes)) +
  scale_fill_gradient(low = "white", high = "steelblue")

# Save for the report
ggsave("results/venn_DEGs.png", p, width = 7, height = 6, dpi = 300)
p




# One gene symbol per line, no header - the format STRING accepts
# na.omit drops genes that have no symbol; unique removes duplicates
writeLines(unique(na.omit(deg_stage$gene_name)), "results/genes_stage.txt")
writeLines(unique(na.omit(deg_hist$gene_name)),  "results/genes_histology.txt")
writeLines(unique(na.omit(deg_age$gene_name)),   "results/genes_age.txt")




# Save the full core-111 table (stage-contrast stats) for the report
write.csv(deg_stage[core, ], "results/core_111_genes.csv")

# Save the core gene SYMBOLS for STRING (protein-coding ones will map, ncRNAs won't - that's expected)
writeLines(unique(na.omit(deg_stage[core, "gene_name"])), "results/genes_core111.txt")
