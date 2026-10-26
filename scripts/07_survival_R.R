# ---- 07_survival.R : do our candidate genes predict patient outcome? ----

# ---- Block 1: setup ----
setwd("D:/DES/NGS/Project/analysis")    # work from the analysis folder

install.packages("survminer")           # pretty survival plots (one-time install, from CRAN not Bioconductor)

library(TCGAbiolinks)                   # to download the clinical table
library(survival)                       # the statistics engine (comes with R already)
library(survminer)                      # the nice Kaplan-Meier plot functions



# ---- Block 2: download the clinical table ----
clin <- GDCquery_clinic("TCGA-CESC", "clinical")   # one row per patient: survival, stage, age, treatment...

dim(clin)                      # expect roughly 307 patients x ~70-100 columns
table(clin$vital_status)       # how many patients are Alive vs Dead



# ---- Block 3: build time and status ----
clin$status <- ifelse(clin$vital_status == "Dead", 1, 0)   # event: 1 = died, 0 = alive (censored)

clin$time <- ifelse(clin$vital_status == "Dead",
                    clin$days_to_death,            # dead  -> we know exactly when
                    clin$days_to_last_follow_up)   # alive -> observed at least this long

summary(clin$time)       # check: should be positive numbers; note how many NA's
table(clin$status)       # should repeat 235 / 72




# ---- Block 4: first KM curve - by stage (sanity check) ----
clin$stage_group <- ifelse(grepl("III|IV", clin$figo_stage), "Late", "Early")   # same rule as our DEA
clin$stage_group[is.na(clin$figo_stage)] <- NA                                  # missing stage stays missing

fit <- survfit(Surv(time, status) ~ stage_group, data = clin)   # the Kaplan-Meier computation

ggsurvplot(fit,
           pval       = TRUE,      # print the log-rank p-value on the plot
           risk.table = TRUE,      # table below: how many patients still "at risk" over time
           legend.title = "Stage")






# ---- Block 5: load expression + gene name lookup ----
prep   <- readRDS("data/CESC_prepared.rds")    # the clean counts from script 02
counts <- prep$counts                          # genes (ENSG ids) x 304 tumours

lookup <- read.csv("results/DEG_stage_full.csv", row.names = 1)   # table with ENSG rownames + gene_name column
head(lookup$gene_name)                                             # peek - symbols should print

# ---- Block 6: KM curve for KRT20 ----
gene <- "KRT20"
ens  <- rownames(lookup)[match(gene, lookup$gene_name)]   # translate symbol -> ENSG id
ens                                                       # check: should print one ENSG000... id

expr  <- log2(counts[ens, ] + 1)                          # KRT20 expression per tumour, log scale
group <- ifelse(expr >= median(expr), "High", "Low")      # split the 304 tumours at the median

km <- data.frame(
  patient = substr(colnames(counts), 1, 12),   # cut barcode to patient id (TCGA-XX-XXXX)
  group   = group                              # High or Low KRT20
)
km <- merge(km, clin[, c("submitter_id", "time", "status")],   # attach survival columns
            by.x = "patient", by.y = "submitter_id")
nrow(km)                                       # how many patients matched - expect ~304

fit_g <- survfit(Surv(time, status) ~ group, data = km)   # KM: does High KRT20 differ from Low?
ggsurvplot(fit_g, pval = TRUE, legend.title = gene)



# ---- Block 7: KM for all six candidates ----
genes_to_test <- c("KRT20", "CDH17", "GABRA2", "GRIN2B", "EGFR", "ERBB2")

pvals <- data.frame(gene = genes_to_test, p_median = NA)    # empty results table to fill

for (g in genes_to_test) {                                  # repeat for each gene:
  ens   <- rownames(lookup)[match(g, lookup$gene_name)]     #   symbol -> ENSG id
  expr  <- log2(counts[ens, ] + 1)                          #   expression, log scale
  group <- ifelse(expr >= median(expr), "High", "Low")      #   median split
  km <- data.frame(patient = substr(colnames(counts), 1, 12), group = group)
  km <- merge(km, clin[, c("submitter_id", "time", "status")],
              by.x = "patient", by.y = "submitter_id")      #   attach survival
  fit_g <- survfit(Surv(time, status) ~ group, data = km)   #   KM fit
  pvals$p_median[pvals$gene == g] <- surv_pvalue(fit_g, data = km)$pval   # store p-value
  pl <- ggsurvplot(fit_g, pval = TRUE, legend.title = g)    #   make the plot
  ggsave(paste0("results/KM_", g, "_median.png"), pl$plot,  #   save as KM_GENE_median.png
         width = 7, height = 6, dpi = 300)
}

pvals                                            # the six p-values side by side
write.csv(pvals, "results/KM_pvalues_median.csv", row.names = FALSE)   # save the summary



# ---- Block 8: top-third vs bottom-third test (ma'am's 0.67/0.33 thresholds) ----
test_tertile <- function(g) {                             # define a recipe named test_tertile
  ens  <- rownames(lookup)[match(g, lookup$gene_name)]    # symbol -> ENSG id
  expr <- log2(counts[ens, ] + 1)                         # expression per tumour
  hi <- expr >= quantile(expr, 0.67)                      # TRUE for top 33% of tumours
  lo <- expr <= quantile(expr, 0.33)                      # TRUE for bottom 33%
  km <- data.frame(patient = substr(colnames(counts), 1, 12),
                   group   = ifelse(hi, "High", ifelse(lo, "Low", NA)))  # middle third = NA
  km <- merge(km, clin[, c("submitter_id", "time", "status")],
              by.x = "patient", by.y = "submitter_id")    # attach survival
  km <- km[!is.na(km$group), ]                            # drop the middle third
  fit <- survfit(Surv(time, status) ~ group, data = km)   # KM on the extremes
  pl  <- ggsurvplot(fit, pval = TRUE, legend.title = g)
  ggsave(paste0("results/KM_", g, "_tertile.png"), pl$plot, width = 7, height = 6, dpi = 300)
  surv_pvalue(fit, data = km)$pval                        # the function hands back the p-value
}

test_tertile("GABRA2")   # the borderline one - moment of truth
test_tertile("CDH17")
test_tertile("EGFR")




# ---- Block 9: Cox model - does EGFR matter AFTER stage and histology? ----
ens_egfr <- rownames(lookup)[match("EGFR", lookup$gene_name)]    # EGFR's ENSG id

cox_df <- data.frame(
  patient = substr(colnames(counts), 1, 12),           # patient id from barcode
  EGFR    = as.numeric(log2(counts[ens_egfr, ] + 1)),  # EGFR expression, continuous (no splitting!)
  stage   = prep$groups$stage_group,                   # Early/Late (made in script 02)
  hist    = prep$groups$hist_group                     # Squamous/Adeno (made in script 02)
)
cox_df <- merge(cox_df, clin[, c("submitter_id", "time", "status")],
                by.x = "patient", by.y = "submitter_id")   # attach survival

cox_fit <- coxph(Surv(time, status) ~ EGFR + stage + hist, data = cox_df)   # one model, three factors
summary(cox_fit)                                           # the verdict table


#stageLate    exp(coef) = 0.69    p = 0.28

# ---- Block 10: detective work - KM and Cox disagree on stage ----

# Check 1: do the stage labels used in Cox (from expression data)
# agree with the stage labels used in the KM (from clinical table)?
chk <- merge(
  data.frame(patient    = substr(colnames(counts), 1, 12),
             stage_expr = prep$groups$stage_group),      # label that went into Cox
  data.frame(patient    = clin$submitter_id,
             stage_clin = clin$stage_group),             # label that went into the KM
  by = "patient")
table(chk$stage_expr, chk$stage_clin, useNA = "ifany")   # want: big numbers ONLY on the diagonal

# Check 2: stage alone, in the same 248 patients the Cox used
summary(coxph(Surv(time, status) ~ stage, data = cox_df))



# ---- Block 11: fix the Cox - pair by NAME, never by position ----

# proof of the bug: same samples, different order?
all(sort(prep$groups$barcode) == sort(colnames(counts)))   # TRUE  = same sample set
all(prep$groups$barcode == colnames(counts))               # FALSE = but different order!

idx <- match(colnames(counts), prep$groups$barcode)   # for each counts column, find its row in groups

cox_df <- data.frame(
  patient = substr(colnames(counts), 1, 12),           # patient id
  EGFR    = as.numeric(log2(counts[ens_egfr, ] + 1)),  # EGFR expression
  stage   = prep$groups$stage_group[idx],              # label fetched by matching barcode
  hist    = prep$groups$hist_group[idx]                # same for histology
)
cox_df <- merge(cox_df, clin[, c("submitter_id", "time", "status")],
                by.x = "patient", by.y = "submitter_id")

# re-check the labels agree with the clinical table now
table(cox_df$stage,
      ifelse(grepl("III|IV", clin$figo_stage[match(cox_df$patient, clin$submitter_id)]), "Late", "Early"),
      useNA = "ifany")                                 # want a clean diagonal this time

cox_fit <- coxph(Surv(time, status) ~ EGFR + stage + hist, data = cox_df)   # the honest Cox
summary(cox_fit)


# ---- Block 12: proportional hazards check ----
cox.zph(cox_fit)   # tests if each factor's effect is stable over time (Cox's key assumption)




# ---- Block 13: Cox for all six candidates ----
cox_results <- data.frame(gene = genes_to_test,            # table to fill
                          HR = NA, lower95 = NA, upper95 = NA, p = NA)

for (g in genes_to_test) {                                 # same recipe as EGFR, per gene:
  ens <- rownames(lookup)[match(g, lookup$gene_name)]      #   symbol -> ENSG
  df <- data.frame(
    patient = substr(colnames(counts), 1, 12),             #   patient id
    expr    = as.numeric(log2(counts[ens, ] + 1)),         #   expression, continuous
    stage   = prep$groups$stage_group[idx],                #   labels via idx = matched by NAME (our fix!)
    hist    = prep$groups$hist_group[idx]
  )
  df <- merge(df, clin[, c("submitter_id", "time", "status")],
              by.x = "patient", by.y = "submitter_id")     #   attach survival
  fit <- coxph(Surv(time, status) ~ expr + stage + hist, data = df)   # adjusted Cox
  s <- summary(fit)                                        #   pull out the gene's numbers:
  cox_results[cox_results$gene == g, c("HR", "lower95", "upper95", "p")] <-
    c(s$conf.int["expr", c("exp(coef)", "lower .95", "upper .95")],   # hazard ratio + 95% CI
      s$coefficients["expr", "Pr(>|z|)"])                             # adjusted p-value
}

cox_results                                                # the report table
write.csv(cox_results, "results/cox_results.csv", row.names = FALSE)   # save it



# ---- Block 14: ROC - can each gene classify the tumour groups? ----
install.packages("pROC")      # one-time install (CRAN)
library(pROC)                 # the ROC/AUC tool (also used in ma'am's 2026 paper)

roc_results <- data.frame(gene = genes_to_test, AUC_stage = NA, AUC_hist = NA)   # table to fill

for (g in genes_to_test) {
  ens   <- rownames(lookup)[match(g, lookup$gene_name)]   # symbol -> ENSG
  expr  <- as.numeric(log2(counts[ens, ] + 1))            # expression, continuous
  stage <- prep$groups$stage_group[idx]                   # by-NAME matched labels (the fix)
  hist  <- prep$groups$hist_group[idx]
  roc_results$AUC_stage[roc_results$gene == g] <-
    auc(roc(stage, expr, levels = c("Early", "Late"), quiet = TRUE))      # AUC: Early vs Late
  roc_results$AUC_hist[roc_results$gene == g] <-
    auc(roc(hist, expr, levels = c("Adeno", "Squamous"), quiet = TRUE))   # AUC: Adeno vs Squamous
}

roc_results                                               # the scorecard
write.csv(roc_results, "results/roc_results.csv", row.names = FALSE)     # save it




# ---- Block 15: figure - EGFR's histology ROC curve ----
ens_egfr <- rownames(lookup)[match("EGFR", lookup$gene_name)]    # EGFR's ENSG id
expr_e   <- as.numeric(log2(counts[ens_egfr, ] + 1))             # EGFR expression
hist_g   <- prep$groups$hist_group[idx]                          # histology labels (by-name matched)

roc_egfr <- roc(hist_g, expr_e, levels = c("Adeno", "Squamous"), quiet = TRUE)   # the ROC object

png("results/ROC_EGFR_histology.png", width = 7, height = 7, units = "in", res = 300)  # open image file
plot(roc_egfr, print.auc = TRUE, main = "EGFR: Squamous vs Adenocarcinoma")            # draw curve + AUC
dev.off()                                                                               # close/save the file

plot(roc_egfr, print.auc = TRUE, main = "EGFR: Squamous vs Adenocarcinoma")  # show it in the Plots pane too
