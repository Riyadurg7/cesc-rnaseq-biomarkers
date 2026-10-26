setwd("D:/DES/NGS/Project/analysis")
library(TCGAbiolinks)
library(ggplot2)

prep   <- readRDS("data/CESC_prepared.rds")
counts <- prep$counts
groups <- prep$groups

early_bc <- groups$barcode[groups$stage_group == "Early" & !is.na(groups$stage_group)]
late_bc  <- groups$barcode[groups$stage_group == "Late"  & !is.na(groups$stage_group)]
length(early_bc); length(late_bc)   # check: 196 and 54

deg_stage <- TCGAanalyze_DEA(
  mat1 = counts[, early_bc],
  mat2 = counts[, late_bc],
  Cond1type = "Early",
  Cond2type = "Late",
  fdr.cut = 1, logFC.cut = 0,   # keep all genes; we filter in the next block
  method = "glmLRT"
)
dim(deg_stage)
head(deg_stage)

head(deg_stage[order(deg_stage$FDR), ], 10)

top_gene <- rownames(deg_stage)[which.max(abs(deg_stage$logFC))]
deg_stage[top_gene, ]
boxplot(counts[top_gene, early_bc], counts[top_gene, late_bc],
        names = c("Early", "Late"), main = top_gene)

#Redo and check the most significant gene instead
top_gene <- rownames(deg_stage)[which.min(deg_stage$FDR)]
deg_stage[top_gene, ]
boxplot(log2(counts[top_gene, early_bc] + 1),
        log2(counts[top_gene, late_bc] + 1),
        names = c("Early", "Late"), main = top_gene, ylab = "log2(counts + 1)")
deg_stage[top_gene, ]


# ---- Block 5: label every gene as UP / DOWN / NotSig ----

# Start by marking ALL genes as "NotSig" (not significant)
deg_stage$status <- "NotSig"

# A gene is "UP" if it changes more than 4-fold (logFC > 2)
# AND the corrected p-value (FDR) is below 0.05
deg_stage$status[deg_stage$logFC >  2 & deg_stage$FDR < 0.05] <- "UP"

# A gene is "DOWN" if it changes more than 4-fold in the other
# direction (logFC < -2) AND FDR < 0.05
# (this is the line where ma'am's script had the bug: logFC < 2 instead of < -2)
deg_stage$status[deg_stage$logFC < -2 & deg_stage$FDR < 0.05] <- "DOWN"

# Count how many genes fall into each category - paste me this
table(deg_stage$status)

# Save the full table (all 45,058 genes) for the record
write.csv(deg_stage, "results/DEG_stage_full.csv")

# Save only the significant genes - this is the list that goes to Cytoscape later
write.csv(deg_stage[deg_stage$status != "NotSig", ], "results/DEG_stage_significant.csv")



# ---- Contrast 2: Squamous vs Adeno ----

# Pick the barcodes (patient IDs) of each histology group
squam_bc <- groups$barcode[groups$hist_group == "Squamous" & !is.na(groups$hist_group)]
adeno_bc <- groups$barcode[groups$hist_group == "Adeno"    & !is.na(groups$hist_group)]

# Check the group sizes: should be 246 and 50
length(squam_bc); length(adeno_bc)

# Run the differential expression test (keep all genes, filter later)
deg_hist <- TCGAanalyze_DEA(
  mat1 = counts[, squam_bc],      # group 1 = Squamous
  mat2 = counts[, adeno_bc],      # group 2 = Adeno
  Cond1type = "Squamous",
  Cond2type = "Adeno",
  fdr.cut = 1, logFC.cut = 0,     # no filtering inside the function
  method = "glmLRT"
)

# Direction check: most significant gene, numbers + log-scale boxplot
top_gene <- rownames(deg_hist)[which.min(deg_hist$FDR)]
deg_hist[top_gene, ]
boxplot(log2(counts[top_gene, squam_bc] + 1),
        log2(counts[top_gene, adeno_bc] + 1),
        names = c("Squamous", "Adeno"), main = top_gene, ylab = "log2(counts + 1)")
# Expected: positive logFC = higher in ADENO (group 2 again) - confirm with your eyes,
# then write the conclusion here as a comment

# Label every gene UP / DOWN / NotSig (same cutoffs as before)
deg_hist$status <- "NotSig"
deg_hist$status[deg_hist$logFC >  2 & deg_hist$FDR < 0.05] <- "UP"
deg_hist$status[deg_hist$logFC < -2 & deg_hist$FDR < 0.05] <- "DOWN"
table(deg_hist$status)   # paste me this

# Save full table + significant-only table
write.csv(deg_hist, "results/DEG_histology_full.csv")
write.csv(deg_hist[deg_hist$status != "NotSig", ], "results/DEG_histology_significant.csv")

# Volcano plot
p <- ggplot(deg_hist, aes(x = logFC, y = -log10(FDR), color = status)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = c(UP = "red", DOWN = "blue", NotSig = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Squamous vs Adenocarcinoma (TCGA-CESC)",
       x = "log2 fold change", y = "-log10 FDR")
ggsave("results/volcano_histology.png", p, width = 7, height = 5, dpi = 300)
p
                        
                            
                            
                            
                            
                            # ---- Contrast 3: Younger vs Older ----
                          
                          # Pick the barcodes of each age group
                          young_bc <- groups$barcode[groups$age_group == "Younger" & !is.na(groups$age_group)]
                          old_bc   <- groups$barcode[groups$age_group == "Older"   & !is.na(groups$age_group)]
                          
                          # Check: should be 178 and 121
                          length(young_bc); length(old_bc)
                          
                          # Differential expression test
                          deg_age <- TCGAanalyze_DEA(
                            mat1 = counts[, young_bc],
                            mat2 = counts[, old_bc],
                            Cond1type = "Younger",
                            Cond2type = "Older",
                            fdr.cut = 1, logFC.cut = 0,
                            method = "glmLRT"
                          )
                          
                          # Label UP / DOWN / NotSig
                          deg_age$status <- "NotSig"
                          deg_age$status[deg_age$logFC >  2 & deg_age$FDR < 0.05] <- "UP"
                          deg_age$status[deg_age$logFC < -2 & deg_age$FDR < 0.05] <- "DOWN"
                          table(deg_age$status)   # paste me this too
                          
                          # Save tables
                          write.csv(deg_age, "results/DEG_age_full.csv")
                          write.csv(deg_age[deg_age$status != "NotSig", ], "results/DEG_age_significant.csv")
                          
                          # Volcano plot
                          p <- ggplot(deg_age, aes(x = logFC, y = -log10(FDR), color = status)) +
                            geom_point(size = 0.8) +
                            scale_color_manual(values = c(UP = "red", DOWN = "blue", NotSig = "grey70")) +
                            geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
                            geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
                            theme_minimal() +
                            labs(title = "Younger vs Older at diagnosis (TCGA-CESC)",
                                 x = "log2 fold change", y = "-log10 FDR")
                          ggsave("results/volcano_age.png", p, width = 7, height = 5, dpi = 300)
                          p





# ---- Block 6: volcano plot ----
# Every dot = one gene. Right side = up in one group, left = up in the other,
# height = statistical significance. Colored dots = your significant genes.

p <- ggplot(deg_stage, aes(x = logFC, y = -log10(FDR), color = status)) +
  geom_point(size = 0.8) +                                          # one dot per gene
  scale_color_manual(values = c(UP = "red", DOWN = "blue", NotSig = "grey70")) +  # colors per category
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +          # the logFC cutoff lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +      # the FDR cutoff line
  theme_minimal() +                                                 # clean white background
  labs(title = "Early vs Late stage (TCGA-CESC)",
       x = "log2 fold change", y = "-log10 FDR")

# Save as a high-resolution image into results/
ggsave("results/volcano_stage.png", p, width = 7, height = 5, dpi = 300)

# Also show it in the Plots panel
p