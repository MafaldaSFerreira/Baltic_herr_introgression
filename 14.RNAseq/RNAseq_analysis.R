# ============================================================================== #
# RNAseq analysis: Effect of temperature on gene expression in Atlantic herring
# Main effect: TEMPERATURE (SEVEN vs TEN degrees)
# Covariates: Chr12 genotype (GENO), SEX, TANK
# Additional tests: Chr12 genotype effect, Temperature x Genotype interaction
# Output: DE results + introgressed gene expression summary table
# ============================================================================== #

library(edgeR)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(biomaRt)

# ============================================================================== #
# 0. ANNOTATION ####
# ============================================================================== #

mart109 <- useMart("ensembl",
                   dataset = "charengus_gene_ensembl",
                   host = "https://feb2023.archive.ensembl.org")
attr <- c("ensembl_gene_id", "external_gene_name", "chromosome_name",
          "start_position", "end_position", "description")
regions <- getBM(attributes=attr, mart=mart109)
colnames(regions) <- c("ENSEMBLID", "NAME", "Chrom", "Start", "End", "Description")
regions$Chrom <- as.numeric(regions$Chrom)

# ============================================================================== #
# 1. LOAD DATA ####
# ============================================================================== #

# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
# Resolved to an absolute path before we setwd() below, so it stays valid
# regardless of the working directory.
FIGSHARE_ROOT <- normalizePath("../Figshare")

save.image(file.path(FIGSHARE_ROOT, "14.RNAseq/RNAanalysis_Environment_20260723.RData"))

#dir.create("results", showWarnings=FALSE, recursive=TRUE)
# Raw STAR alignment output -- not part of the Figshare deposit (too large to
# share). Set this to wherever your own STAR results live.
setwd("<path_to_your_STAR_alignment_results>")
results_dir <- file.path(FIGSHARE_ROOT, "14.RNAseq")

design <- read.table("design_genes.txt", header=TRUE, sep="\t")

# Ensure factors are correct
design$TEMP <- factor(design$TEMP, levels=c("SEVEN", "TEN"))  # SEVEN = reference
design$GENO <- factor(design$GENO)                             # Chr12 genotype
design$SEX  <- factor(design$SEX)
design$TANK <- factor(design$TANK)

# Load counts
df <- readDGE(files   = design$COUNTS,
              path    = ".",
              labels  = design$SHORTNAME,
              columns = c(1, 5),
              group   = design$TEMP)

# ============================================================================== #
# 2. FILTERING AND NORMALISATION ####
# ============================================================================== #

keep <- filterByExpr(df)
cat("Genes kept after filtering:", sum(keep), "\n")
#Genes kept after filtering: 20129 
df   <- df[keep, , keep.lib.sizes=FALSE]
df   <- calcNormFactors(df)

# ============================================================================== #
# 3. MDS PLOT ####
# ============================================================================== #

plotMDS_df <- plotMDS(df, plot=FALSE)
mds_df <- data.frame(
  X1        = plotMDS_df$x,
  X2        = plotMDS_df$y,
  SHORTNAME = design$SHORTNAME
) %>% left_join(design, by="SHORTNAME")

variance <- plotMDS_df$var.explained[1:2]

#### MDS PLOT ####
mds_plot <- ggplot(mds_df, aes(x=X1, y=X2, color=TEMP, shape=GENO)) +
  geom_point(size=2) +
  xlab(paste0("Dim 1 (", round(variance[1]*100, 1), "%)")) +
  ylab(paste0("Dim 2 (", round(variance[2]*100, 1), "%)")) +
  labs(color="Temperature", shape="Chr12 genotype") +
  theme_classic()+
  theme(axis.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9))

mds_plot_tank_sex <- ggplot(mds_df, aes(x=X1, y=X2, color=TANK, shape=SEX)) +
  geom_point(size=2) +
  xlab(paste0("Dim 1 (", round(variance[1]*100, 1), "%)")) +
  ylab(paste0("Dim 2 (", round(variance[2]*100, 1), "%)")) +
  labs(color="Tank", shape="Sex") +
  theme_classic()+
  theme(axis.title = element_text(size=9),
        legend.text = element_text(size=9),
        legend.title = element_text(size=9))


ggsave(mds_plot, filename=paste0(results_dir,"/figures/MDS_temperature_analysis.pdf"), height=5, width=7)
ggsave(mds_plot_tank_sex, filename=paste0(results_dir,"/figures/MDS_tank_sex_analysis.pdf"), height=5, width=7)

# ============================================================================== #
# 4. MAIN MODEL: TEMPERATURE + COVARIATES ####
# Model: ~ TEMP + GENO
# SEX and TANK ommitted because there is no obvious batch effect related to them
# TEMP is the main effect of interest
# GENO (Chr12), SEX, TANK are covariates
# ============================================================================== #

design_matrix <- model.matrix(~ TEMP + GENO, data=design)
cat("Design matrix columns:\n")
print(colnames(design_matrix))

# Check rank
cat("Matrix rank:", qr(design_matrix)$rank, 
    "vs columns:", ncol(design_matrix), "\n")
#Matrix rank: 3 vs columns: 3
# Estimate dispersions
df <- estimateDisp(df, design_matrix, robust=TRUE)
cat("Common BCV:", sqrt(df$common.dispersion), "\n")
#Common BCV: 0.2316585 
plotBCV(df)

# Fit GLM
fit <- glmQLFit(df, design_matrix, robust=TRUE)


# ============================================================================== #
# 5. TEST 1: MAIN EFFECT OF TEMPERATURE ####
# ============================================================================== #

qlf_TEMP <- glmQLFTest(fit, coef="TEMPTEN")

# Summary
cat("\n--- Temperature effect (SEVEN vs TEN) ---\n")
print(summary(decideTests(qlf_TEMP)))
# TEMPTEN
# Down      5181
# NotSig   10373
# Up        4575

# Full results table
res_TEMP <- topTags(qlf_TEMP, n=Inf, sort.by="PValue")$table
res_TEMP$ENSEMBLID <- rownames(res_TEMP)
res_TEMP <- left_join(res_TEMP, regions, by="ENSEMBLID")

# Significant DE genes
DE_TEMP_FDR05 <- res_TEMP %>% filter(FDR < 0.05)
cat("Significant DE genes by temperature (FDR < 0.05):", nrow(DE_TEMP_FDR05), "\n")
#Significant DE genes by temperature (FDR < 0.05): 9756 
DE_TEMP_FDR05_LFC1 <- res_TEMP %>% filter(FDR < 0.05 & abs(logFC) > 1)
cat("Significant DE genes by temperature (FDR < 0.05):", nrow(DE_TEMP_FDR05_LFC1), "\n")
#Significant DE genes by temperature (FDR < 0.05): 659 

write.csv(res_TEMP, paste0(results_dir,"/results/DE_temperature_full_results.csv"), row.names=FALSE)
write.csv(DE_TEMP_FDR05,  paste0(results_dir,"/results/DE_temperature_significant_FDR05.csv"), row.names=FALSE)
write.csv(DE_TEMP_FDR05_LFC1,  paste0(results_dir,"/results/DE_temperature_significant_FDR05_LogFC1.csv"), row.names=FALSE)

# ============================================================================== #
# 6. TEST 2: EFFECT OF CHR12 GENOTYPE (adjusting for temperature) ####
# ============================================================================== #

qlf_GENO <- glmQLFTest(fit, coef=grep("GENO", colnames(design_matrix), value=TRUE))

cat("\n--- Chr12 genotype effect ---\n")
print(summary(decideTests(qlf_GENO)))
# GENONorth
# Down         330
# NotSig     19454
# Up           345

res_GENO <- topTags(qlf_GENO, n=Inf, sort.by="PValue")$table
res_GENO$ENSEMBLID <- rownames(res_GENO)
res_GENO <- left_join(res_GENO, regions, by="ENSEMBLID")

DE_GENO_FDR05 <- res_GENO %>% filter(FDR < 0.05)
cat("Significant DE genes by Chr12 genotype (FDR < 0.05):", nrow(DE_GENO_FDR05), "\n")
#Significant DE genes by Chr12 genotype (FDR < 0.05): 675 
DE_GENO_FDR05_logFC1 <- res_GENO %>% filter(FDR < 0.05 & abs(logFC) > 1)
cat("Significant DE genes by Chr12 genotype (FDR < 0.05):", nrow(DE_GENO_FDR05_logFC1), "\n")
#Significant DE genes by Chr12 genotype (FDR < 0.05): 146

write.csv(res_GENO, paste0(results_dir,"/results/DE_chr12genotype_full_results.csv"), row.names=FALSE)
write.csv(DE_GENO_FDR05, paste0(results_dir,"/results/DE_chr12genotype_full_results_FDR05.csv"), row.names=FALSE)
write.csv(DE_GENO_FDR05_logFC1, paste0(results_dir,"/results/DE_chr12genotype_full_results_FDR05_logFC1.csv"), row.names=FALSE)

# ============================================================================== #
# 7. TEST 3: TEMPERATURE x CHR12 GENOTYPE INTERACTION ####
# Fit a separate model with the interaction term
# ============================================================================== #

design_matrix_int <- model.matrix(~ TEMP * GENO, data=design)

df_int <- estimateDisp(df, design_matrix_int, robust=TRUE)
fit_int <- glmQLFit(df_int, design_matrix_int, robust=TRUE)

# Test interaction term
int_coefs <- grep("TEMPTEN:GENO", colnames(design_matrix_int), value=TRUE)
qlf_INT <- glmQLFTest(fit_int, coef=int_coefs)

cat("\n--- Temperature x Chr12 genotype interaction ---\n")
print(summary(decideTests(qlf_INT)))
# TEMPTEN:GENONorth
# Down                   0
# NotSig             20129
# Up                     0
# There is no interaction

res_INT <- topTags(qlf_INT, n=Inf, sort.by="PValue")$table
res_INT$ENSEMBLID <- rownames(res_INT)
res_INT <- left_join(res_INT, regions, by="ENSEMBLID")

DE_INT <- res_INT %>% filter(FDR < 0.05)
cat("Significant interaction genes (FDR < 0.05):", nrow(DE_INT), "\n")
#Significant interaction genes (FDR < 0.05): 0 

write.csv(res_INT, paste0(results_dir,"/results/DE_temperature_genotype_interaction_results.csv"), row.names=FALSE)

# ============================================================================== #
# 8. INTROGRESSED GENE SUMMARY TABLE ####
# For each introgressed gene: is it expressed? is it DE by temperature?
# ============================================================================== #

# Read annotation file with ENSEMBL IDs and gene coordinates
scan1_match_df <- read.table(
  file.path(FIGSHARE_ROOT, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.modified.txt"),
  sep="\t", header=TRUE, row.names=NULL)

# Read introgression region coordinates
intro_reg <- read.table(
  file.path(FIGSHARE_ROOT, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt"),
  sep="\t", header=TRUE, row.names=NULL)

# Build region labels from coordinates (e.g. "Chr10:25.04-25.22")
intro_reg <- intro_reg %>%
  mutate(
    chr_num    = str_remove(seqnames, "chr"),
    region_key = paste0(chr_num, "_", start, "_", end),
    region     = paste0("Chr", chr_num, ":",
                        round(start/1e6, 2), "-",
                        round(end/1e6, 2))
  )

# Use GenomicRanges to map genes to introgression regions via overlap
# maxgap=20000 includes genes within 20kb of introgression region boundaries
library(GenomicRanges)
scan1_match_df_gr <- makeGRangesFromDataFrame(scan1_match_df)
seqlevelsStyle(scan1_match_df_gr) <- "NCBI"
intro_reg_gr <- makeGRangesFromDataFrame(intro_reg)
seqlevelsStyle(intro_reg_gr) <- "NCBI"

overlaps <- findOverlaps(intro_reg_gr, scan1_match_df_gr, maxgap=20000)

# Assign region label to each gene via overlap indices
scan1_match_df$region <- NA
scan1_match_df$region[overlaps@to] <- intro_reg$region[overlaps@from]

cat("Genes with region assigned:", sum(!is.na(scan1_match_df$region)), "\n")
#Genes with region assigned: 166 
cat("Genes without region (outside 20kb):", sum(is.na(scan1_match_df$region)), "\n")
#Genes without region (outside 20kb): 0 

# Build introgressed gene table:
introgressed_genes <- scan1_match_df %>%
  filter(!is.na(region)) %>%
  # Exclude: genes outside introgression regions (this is unnecessary step, but just for precaution)
  filter(!is.na(ensembl_gene_id) & ensembl_gene_id != "") %>%
  # Exclude: unannotated transcripts, pseudogenes, and housekeeping small RNAs
  # (SNORD = small nucleolar RNAs, U5 = small nuclear RNA — constitutive/housekeeping)
  # Include: lncRNA and snoRNA as these may have regulatory roles in adaptation
  filter(!grepl("novel|PSEUDOGENE",
                external_gene_name, ignore.case=TRUE)) %>%
  dplyr::rename(ENSEMBLID = ensembl_gene_id,
                gene_name = external_gene_name) %>%
  mutate(
    gene_name = toupper(gene_name),
    # Classify gene biotype for stratified reporting
    gene_class = case_when(
      grepl("lncRNA", gene_name, ignore.case=TRUE) ~ "lncRNA",
      grepl("snoRNA|snRNA|^SNORD|^U[0-9]", gene_name, ignore.case=TRUE) ~ "ncRNA",
      TRUE ~ "protein_coding"
    )
  ) %>%
  distinct(ENSEMBLID, .keep_all=TRUE) %>%
  dplyr::select(ENSEMBLID, gene_name, gene_class, region)

cat("Total genes (all classes):", nrow(introgressed_genes), "\n")
# Total genes (all classes): 154 
cat("  Protein-coding:", sum(introgressed_genes$gene_class == "protein_coding"), "\n")
#Protein-coding: 115 
cat("  lncRNA:",         sum(introgressed_genes$gene_class == "lncRNA"), "\n")
#lncRNA: 12 
cat("  ncRNA (snoRNA/snRNA):", sum(introgressed_genes$gene_class == "ncRNA"), "\n")
#ncRNA (snoRNA/snRNA): 27

# ============================================================================== #
# Build expression and DE summary
# ============================================================================== #

introgressed_summary <- introgressed_genes %>%
  mutate(
    # Initial expressed flag based on whether gene passed filterByExpr
    expressed  = ENSEMBLID %in% rownames(df),
    logFC_temp = res_TEMP$logFC[match(ENSEMBLID, res_TEMP$ENSEMBLID)],
    pval_temp  = res_TEMP$PValue[match(ENSEMBLID, res_TEMP$ENSEMBLID)],
    FDR_temp   = res_TEMP$FDR[match(ENSEMBLID, res_TEMP$ENSEMBLID)],
    DE_by_temp = ifelse(!is.na(FDR_temp) & FDR_temp < 0.05, "Yes", "No"),
  )

# Add mean CPM from normalised counts for expressed genes
cpm_matrix <- cpm(df, log=FALSE, normalized.lib.sizes=TRUE)

introgressed_summary$mean_CPM <- sapply(introgressed_summary$ENSEMBLID, function(id) {
  if (id %in% rownames(cpm_matrix)) round(mean(cpm_matrix[id, ]), 2) else NA
})

# ============================================================================== #
# Verify expression status of genes that did not pass filterByExpr
# Check whether they are truly absent or present but below quantification threshold
# ============================================================================== #

df_raw <- readDGE(
  files   = design$COUNTS,
  path    = ".",
  labels  = design$SHORTNAME,
  columns = c(1, 5),
  group   = design$TEMP
)


unexpressed_ids <- introgressed_summary %>%
  filter(expressed == FALSE) %>%
  pull(ENSEMBLID) %>%
  unique()

length(unexpressed_ids)
#38

cat("\n--- Verification of unexpressed genes ---\n")
cat("Checking whether genes are absent from raw counts or filtered due to low counts:\n\n")

# Build low_expressed lookup dynamically from raw counts
# Any gene present in raw counts but not passing filterByExpr is "Low (filtered)"
low_expressed_ids <- c()
low_raw_counts    <- c()

for (g in unexpressed_ids) {
  gname <- introgressed_summary$gene_name[introgressed_summary$ENSEMBLID == g][1]
  if (g %in% rownames(df_raw$counts)) {
    total <- sum(df_raw$counts[g, ])
    cat(gname, "(", g, ") - present in raw counts, total:", round(total, 1), "\n")
    low_expressed_ids          <- c(low_expressed_ids, g)
    low_raw_counts[g]          <- total
  } else {
    cat(gname, "(", g, ") - not found in raw counts at all\n")
  }
}

# SNORNA ( ENSCHAG00000012661 ) - present in raw counts, total: 4 
# LNCRNA ( ENSCHAG00000016191 ) - present in raw counts, total: 81.6 
# SNORNA ( ENSCHAG00000000386 ) - present in raw counts, total: 5 
# ZBTB11 ( ENSCHAG00000001046 ) - present in raw counts, total: 106 
# SNORD14 ( ENSCHAG00000001073 ) - present in raw counts, total: 125.1 
# SNORD14 ( ENSCHAG00000001080 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001084 ) - present in raw counts, total: 211.8 
# SNORD14 ( ENSCHAG00000001094 ) - present in raw counts, total: 252.9 
# SNORD14 ( ENSCHAG00000001109 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001120 ) - present in raw counts, total: 71.9 
# SNORD14 ( ENSCHAG00000001134 ) - present in raw counts, total: 355.1 
# SNORD14 ( ENSCHAG00000001146 ) - present in raw counts, total: 71.9 
# SNORD14 ( ENSCHAG00000001163 ) - present in raw counts, total: 59.1 
# SNORD14 ( ENSCHAG00000001174 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001179 ) - present in raw counts, total: 71.9 
# SNORD14 ( ENSCHAG00000001185 ) - present in raw counts, total: 211.8 
# SNORD14 ( ENSCHAG00000001192 ) - present in raw counts, total: 71.9 
# SNORD14 ( ENSCHAG00000001197 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001202 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001212 ) - present in raw counts, total: 71.9 
# SNORD14 ( ENSCHAG00000001220 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001228 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001260 ) - present in raw counts, total: 163.1 
# SNORD14 ( ENSCHAG00000001267 ) - present in raw counts, total: 75.8 
# LNCRNA ( ENSCHAG00000021695 ) - present in raw counts, total: 64 
# LNCRNA ( ENSCHAG00000022117 ) - present in raw counts, total: 26 
# GNG12 ( ENSCHAG00000022120 ) - present in raw counts, total: 216 
# LNCRNA ( ENSCHAG00000022272 ) - present in raw counts, total: 0 
# LNCRNA ( ENSCHAG00000023281 ) - present in raw counts, total: 6 
# LNCRNA ( ENSCHAG00000008040 ) - present in raw counts, total: 21 
# LNCRNA ( ENSCHAG00000011128 ) - present in raw counts, total: 45 
# PRKAR2B ( ENSCHAG00000007917 ) - present in raw counts, total: 7 
# LY86 ( ENSCHAG00000003470 ) - present in raw counts, total: 335 
# F13A1 ( ENSCHAG00000003517 ) - present in raw counts, total: 26.1 
# F13A1A.1  ( ENSCHAG00000003862 ) - present in raw counts, total: 417.9 
# U5 ( ENSCHAG00000004421 ) - present in raw counts, total: 30 
# CDH20 ( ENSCHAG00000009619 ) - present in raw counts, total: 92.6 
# LNCRNA ( ENSCHAG00000010866 ) - present in raw counts, total: 69 

cat("\nGenes classified as Low (filtered):", length(low_expressed_ids), "\n")
cat("Genes truly absent from counts:",
    length(unexpressed_ids) - length(low_expressed_ids), "\n")
#Genes classified as Low (filtered): 38
#Genes truly absent from counts: 0 

n_samples <- ncol(df_raw$counts)

# ============================================================================== #
# Refine expressed status using verification results
# ============================================================================== #

introgressed_summary <- introgressed_summary %>%
  mutate(
    expressed = case_when(
      ENSEMBLID %in% low_expressed_ids ~ "Low (filtered)",
      expressed == TRUE                 ~ "Yes",
      TRUE                              ~ "No"
    ),
    # mean_CPM is deliberately left as NA for genes that did not pass
    # filterByExpr: CPM from `df` is defined only for retained genes
    # (library sizes recomputed after filtering + TMM norm factors), so a
    # value computed any other way would not be on a comparable scale.
    # Raw counts are reported separately in total_raw_counts instead.
    total_raw_counts = ifelse(
      ENSEMBLID %in% names(low_raw_counts),
      round(unname(low_raw_counts[ENSEMBLID]), 2),
      NA_real_
    )
  )

print(introgressed_summary %>%
        dplyr::select(gene_name, region, expressed, mean_CPM, total_raw_counts,
                      logFC_temp, FDR_temp, DE_by_temp) %>%
        arrange(region))

introgressed_summary %>% filter(expressed == "Yes") %>% summarise(n=n())
# 116

cat("\n======= EXPRESSION SUMMARY BY GENE CLASS =======\n")
for (cls in c("protein_coding", "lncRNA", "ncRNA")) {
  sub <- introgressed_summary %>% filter(gene_class == cls)
  cat(sprintf("\n%s (n=%d):\n", cls, nrow(sub)))
  cat("  Expressed (Yes):",           sum(sub$expressed == "Yes"), "\n")
  cat("  Low expression (filtered):", sum(sub$expressed == "Low (filtered)"), "\n")
  cat("  Not expressed:",             sum(sub$expressed == "No"), "\n")
}
# DE tallies moved to Section 10 — they depend on the restricted FDR

# protein_coding (n=115):
#   Expressed (Yes): 108 
# Low expression (filtered): 7 
# Not expressed: 0 
# 
# lncRNA (n=12):
#   Expressed (Yes): 4 
# Low expression (filtered): 8 
# Not expressed: 0 
# 
# ncRNA (n=27):
#   Expressed (Yes): 4 
# Low expression (filtered): 23 
# Not expressed: 0 


# ============================================================================== #
# 9. VOLCANO PLOT: TEMPERATURE EFFECT ####
# Highlight introgressed genes in black
# ============================================================================== #
# Deleted, moved to section 10 


# ============================================================================== #
# 10. TARGETED DE ANALYSIS: INTROGRESSED GENES ONLY ####
#
# Motivation:
# Sections 5-7 apply Benjamini-Hochberg (BH) FDR correction across ALL expressed
# genes (thousands of tests). Because our a priori biological question concerns
# only the introgressed genes, that genome-wide correction is unnecessarily
# conservative for this gene set.
#
# Here we restrict the multiple-testing correction to the introgressed genes that
# were actually tested. This is valid because the introgressed set is defined
# INDEPENDENTLY of the expression data (from the introgression scan / population
# genetics), so subsetting introduces no selection bias w.r.t. the DE p-values.
#
# IMPORTANT - MODEL FITTING IS UNCHANGED:
# We do NOT refit the model on the subset. Dispersion estimation (estimateDisp)
# and the GLM fit (glmQLFit) still use the FULL gene set (Sections 2-4), so
# empirical-Bayes dispersion shrinkage still borrows strength across all genes.
# We only change the SET of tests over which BH is computed. Raw p-values are
# therefore IDENTICAL to the genome-wide analysis; only the FDR changes.
# ============================================================================== #

# A priori set of introgressed genes (defined by the introgression scan)
# This uses 133 genes, so after filtering out PSEUDOGENES, novel, etc
# Perhaps I should test after with all 166 introgressed genes
intro_ids <- unique(introgressed_genes$ENSEMBLID)

# Helper: subset a genome-wide topTags table to the introgressed genes that were
# actually tested (non-NA p-value), then recompute BH FDR within that subset only.
# Keeps the original genome-wide FDR alongside for direct comparison.
subset_fdr <- function(res, subset_ids) {
  res %>%
    filter(ENSEMBLID %in% subset_ids, !is.na(PValue)) %>%
    arrange(PValue) %>%
    mutate(
      FDR_genomewide   = FDR,                         # original, for comparison
      FDR_introgressed = p.adjust(PValue, method = "BH")  # subset correction
    )
}

# ------------------------------------------------------------------ #
# 10a. Temperature effect within introgressed genes (main test)
# ------------------------------------------------------------------ #
intro_TEMP <- subset_fdr(res_TEMP, intro_ids)

res_TEMP %>%
  filter(ENSEMBLID %in% intro_ids, !is.na(PValue)) %>%
  arrange(PValue) %>%
  mutate(
    FDR_genomewide   = FDR,                         # original, for comparison
    FDR_introgressed = p.adjust(PValue, method = "BH")  # subset correction
  )

cat("\n=== Targeted DE: TEMPERATURE (introgressed genes only) ===\n")
cat("Introgressed genes tested (family size)  :", nrow(intro_TEMP), "\n")
# Introgressed genes tested (family size)  : 116 
# of the 133 genes left above, only 116 were actually tested because some were
# excluded for having too low expression initially.

cat("Significant, genome-wide FDR < 0.05      :",
    sum(intro_TEMP$FDR_genomewide   < 0.05), "\n")

#Significant, genome-wide FDR < 0.05      : 65 

cat("Significant, introgressed-only FDR < 0.05:",
    sum(intro_TEMP$FDR_introgressed < 0.05), "\n")

# Significant, introgressed-only FDR < 0.05: 66 

# Genes that gain significance under the less conservative correction
gained_TEMP <- intro_TEMP %>%
  filter(FDR_genomewide >= 0.05 & FDR_introgressed < 0.05)
cat("Newly significant after subset correction:", nrow(gained_TEMP), "\n")
#Newly significant after subset correction: 1
if (nrow(gained_TEMP) > 0) {
  print(gained_TEMP %>%
          dplyr::select(dplyr::any_of(c("NAME", "logFC", "PValue",
                                        "FDR_genomewide", "FDR_introgressed"))))
}
# NAME     logFC     PValue FDR_genomewide FDR_introgressed
# 1 cd36 0.3137296 0.02688321      0.0546598       0.04724927

write.csv(intro_TEMP,
          paste0(results_dir, "/results/DE_temperature_introgressed_only.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------ #
# 10b. Chr12 genotype effect within introgressed genes
# NOTE: if GENO has >2 levels the test is an omnibus ANOVA-style test, so the
# table has multiple logFC columns (logFC.GENO*) rather than a single logFC.
# PValue / FDR still behave as expected, so the correction below is unaffected.
# ------------------------------------------------------------------ #
intro_GENO <- subset_fdr(res_GENO, intro_ids)
cat("\n=== Targeted DE: CHR12 GENOTYPE (introgressed genes only) ===\n")
cat("Introgressed genes tested (family size)  :", nrow(intro_GENO), "\n")
# 116
cat("Significant, introgressed-only FDR < 0.05:",
    sum(intro_GENO$FDR_introgressed < 0.05), "\n")
# 1
write.csv(intro_GENO,
          paste0(results_dir, "/results/DE_chr12genotype_introgressed_only.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------ #
# 10c. Temperature x genotype interaction within introgressed genes
# ------------------------------------------------------------------ #
intro_INT <- subset_fdr(res_INT, intro_ids)
cat("\n=== Targeted DE: TEMP x GENO INTERACTION (introgressed genes only) ===\n")
cat("Introgressed genes tested (family size)  :", nrow(intro_INT), "\n")
# 116
cat("Significant, introgressed-only FDR < 0.05:",
    sum(intro_INT$FDR_introgressed < 0.05), "\n")
# 0

write.csv(intro_INT,
          paste0(results_dir, "/results/DE_interaction_introgressed_only.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------ #
# 10d. Fold the introgressed-only FDR back into the summary table (Section 8)
# so downstream tables/figures can flag DE using the less conservative threshold.
# ------------------------------------------------------------------ #
introgressed_summary <- introgressed_summary %>%
  left_join(
    intro_TEMP %>%
      dplyr::select(ENSEMBLID, FDR_temp_introgressed = FDR_introgressed),
    by = "ENSEMBLID"
  ) %>%
  mutate(
    DE_by_temp_introgressed = ifelse(
      !is.na(FDR_temp_introgressed) & FDR_temp_introgressed < 0.05, "Yes", "No"
    )
  )

cat("\n--- Temperature DE within introgressed genes: genome-wide vs targeted ---\n")
cat("DE by temp (genome-wide FDR):   ",
    sum(introgressed_summary$DE_by_temp == "Yes", na.rm = TRUE), "\n")
# 65
cat("DE by temp (introgressed FDR):  ",
    sum(introgressed_summary$DE_by_temp_introgressed == "Yes", na.rm = TRUE), "\n")
# 66

# This will be the new Supplementary Table 8:
write.csv(introgressed_summary,
          paste0(results_dir,
                 "/results/introgressed_genes_expression_temperature_summary_targeted_FDR_introgression.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------ #
# 10e. DE tally by gene class, using the restricted (introgressed-only) FDR
# ------------------------------------------------------------------ #
de_by_class <- introgressed_summary %>%
  group_by(gene_class) %>%
  summarise(
    n_annotated   = n(),
    n_expressed   = sum(expressed == "Yes"),
    n_low         = sum(expressed == "Low (filtered)"),
    DE_genomewide = sum(DE_by_temp == "Yes", na.rm = TRUE),
    DE_restricted = sum(DE_by_temp_introgressed == "Yes", na.rm = TRUE),
    .groups = "drop"
  )

de_by_class <- bind_rows(
  de_by_class,
  de_by_class %>%
    summarise(across(where(is.numeric), sum)) %>%
    mutate(gene_class = "TOTAL")
)

print(as.data.frame(de_by_class))
# gene_class n_annotated n_expressed n_low DE_genomewide DE_restricted
# 1         lncRNA          12           4     8             2             2
# 2          ncRNA          27           4    23             2             2
# 3 protein_coding         115         108     7            61            62
# 4          TOTAL         154         116    38            65            66

####
idx <- rownames(df) %in% intro_ids
camera_res <- camera(df, index = list(introgressed = idx),
       design = design_matrix, contrast = "TEMPTEN")
# NGenes Direction    PValue
# introgressed    116        Up 0.0366559

# Volcano plot: 
volcano_df <- res_TEMP %>%
  left_join(introgressed_summary %>%
              dplyr::select(ENSEMBLID, DE_by_temp_introgressed),
            by = "ENSEMBLID") %>%
  mutate(
    introgressed = ENSEMBLID %in% intro_ids,
    # Background genes: genome-wide FDR (no restricted FDR exists for them)
    sig_bg   = !introgressed & FDR < 0.05,
    # Introgressed genes: restricted FDR
    intro_DE = introgressed &
      !is.na(DE_by_temp_introgressed) &
      DE_by_temp_introgressed == "Yes",
    bg_group = case_when(
      introgressed       ~ NA_character_,
      sig_bg & logFC > 0 ~ "Up",
      sig_bg & logFC < 0 ~ "Down",
      TRUE               ~ "NS"
    ),
    label = ifelse(introgressed & !is.na(NAME) & NAME != "", toupper(NAME), NA)
  )

force_label_genes <- c("THRB", "SEC16B", "IGF2BP1", "PDGFRA", "SLCO1C1")

volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = -log10(PValue))) +
  geom_point(data = filter(volcano_df, bg_group == "NS"),
             colour = "grey85", size = 0.5, alpha = 0.4) +
  geom_point(data = filter(volcano_df, bg_group == "Up"),
             colour = "#f4a582", size = 0.6, alpha = 0.5) +
  geom_point(data = filter(volcano_df, bg_group == "Down"),
             colour = "#92c5de", size = 0.6, alpha = 0.5) +
  # Introgressed genes, fill mapped so a legend is generated
  geom_point(data = filter(volcano_df, introgressed),
             aes(fill = intro_DE),
             colour = "black", shape = 21, size = 2.4, stroke = 0.45) +
  scale_fill_manual(
    values = c(`TRUE` = "#FFC107", `FALSE` = "white"),
    labels = c(`TRUE` = "DE (FDR < 0.05)", `FALSE` = "Not DE"),
    name   = "Introgressed genes"
  ) +
  geom_text_repel(data = filter(volcano_df, introgressed & !is.na(label) &
                                  !label %in% force_label_genes),
                  aes(label = label), size = 2.5, fontface = "bold.italic",
                  max.overlaps = 10, box.padding = 0.5,
                  segment.color = "grey50", segment.size = 0.3,
                  colour = "grey20") +
  geom_text_repel(data = filter(volcano_df, label %in% force_label_genes),
                  aes(label = label), size = 2.5, fontface = "bold.italic",
                  force = 8, force_pull = 0.2, nudge_x = -1, nudge_y = 1,
                  max.overlaps = Inf, box.padding = 0.6,
                  segment.color = "black", segment.size = 0.4,
                  min.segment.length = 0, colour = "black") +
  # Dashed line = genome-wide FDR threshold; applies to background genes only
  geom_hline(yintercept = -log10(max(res_TEMP$PValue[res_TEMP$FDR < 0.05],
                                     na.rm = TRUE)),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  labs(x = "log2 FC (10°C vs 7°C)",
       y = expression(-log[10](P~value))) +
  theme_classic() +
  theme(axis.text  = element_text(size = 9, colour = "black"),
        axis.title = element_text(size = 9, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

# Add a Panel to figure 20 with the results of the enrichment
# Panel c: logFC by gene set — visualises the camera result
lfc_df <- res_TEMP %>%
  filter(!is.na(logFC)) %>%
  mutate(set = factor(ifelse(ENSEMBLID %in% intro_ids,
                             "Introgressed", "Background"),
                      levels = c("Background", "Introgressed")))

# Sensible y-limits: clip the display, not the data (coord_cartesian, not ylim)
y_lim <- quantile(lfc_df$logFC, c(0.001, 0.999), na.rm = TRUE)

panel_lfc <- ggplot(lfc_df, aes(x = set, y = logFC, fill = set)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.3, outlier.colour = "grey70",
               outlier.alpha = 0.4, width = 0.55, linewidth = 0.35) +
  # Show the individual introgressed genes — n is small enough to plot
  geom_jitter(data = filter(lfc_df, set == "Introgressed"),
              width = 0.15, size = 0.9, alpha = 0.55,
              colour = "grey20") +
  # Mean as a diamond: camera tests a shift in mean, not median
  stat_summary(fun = mean, geom = "point", shape = 23,
               size = 2, fill = "white", colour = "black", stroke = 0.4) +
  scale_fill_manual(values = c("Background"   = "grey85",
                               "Introgressed" = "#FFC107")) +
  coord_cartesian(ylim = y_lim) +
  labs(x = NULL, y = "log2 FC (10°C vs 7°C)") +
  theme_classic() +
  theme(axis.text    = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none")

# Top row: MDS + logFC boxplot side by side
top_row <- plot_grid(
  mds_plot,
  mds_plot_tank_sex,
  ncol       = 2,
  labels     = c("a", "b"),
  rel_widths = c(1.1, 0.9)   # MDS wider — it carries a legend
)

# Stack: top row over the volcano
bottom_row <- plot_grid(
  volcano_plot,
  panel_lfc,
  ncol        = 2,
  labels      = c("c", "d"),   
  rel_widths= c(1.5, 0.5)   # volcano gets more vertical space
)

figure <- plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(0.9, 1.1)
)

figure

ggsave(figure,
       filename = paste0(results_dir, "/figures/rnaseq_temperature_summary_figure_20260723.pdf"),
       height = 180, width = 180, units = "mm")

ggsave(figure,
       filename = paste0(results_dir, "/figures/rnaseq_temperature_summary_figure_20260723.png"),
       height = 180, width = 180, units = "mm", dpi = 300)


# We need to update Supplementary Fig. 21 to include the extra DE gene:
library(patchwork)

plot_data <- introgressed_summary %>%
  filter(expressed == "Yes") %>%
  arrange(desc(ifelse(is.na(FDR_temp_introgressed), 1, FDR_temp_introgressed))) %>%
  mutate(
    gene_label = paste0(gene_name, " (", ENSEMBLID, ")"),
    gene_label = factor(gene_label, levels = rev(unique(gene_label))),
    direction  = case_when(
      DE_by_temp_introgressed == "Yes" & logFC_temp > 0 ~ "Up at 10°C",
      DE_by_temp_introgressed == "Yes" & logFC_temp < 0 ~ "Down at 10°C",
      TRUE                                  ~ "NS"
    ),
    # Use logFC as fill, but set NS genes to 0 so they appear grey
    fill_value = case_when(
      direction == "NS" ~ 0,
      TRUE              ~ logFC_temp
    )
  ) %>%
  mutate(display_label = sub(" \\(ENSCHAG.*\\)", "", gene_label))

# Shared fill scale limits — symmetric around 0
fill_lim <- max(abs(plot_data$fill_value), na.rm = TRUE)

make_panel <- function(df, show_legend = FALSE) {
  ggplot(df, aes(x = 0, y = gene_label, fill = fill_value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low      = "#4575b4",   # Down at 10°C
      mid      = "gray90",    # NS
      high     = "#d73027",   # Up at 10°C
      midpoint = 0,
      limits   = c(-fill_lim, fill_lim),
      name     = "log2 FC\n(10°C vs 7°C)",
      breaks   = c(-fill_lim, 0, fill_lim),
      labels   = c(
        paste0("Down (", round(-fill_lim, 1), ")"),
        "NS / 0",
        paste0("Up (", round(fill_lim, 1), ")")
      )
    ) +
    scale_y_discrete(labels = function(x) sub(" \\(ENSCHAG.*\\)", "", x)) +
    scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y       = element_text(size = 7, face = "italic", color = "black"),
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.line         = element_blank(),
      legend.position   = if (show_legend) "right" else "none",
      legend.title      = element_text(size = 8),
      legend.text       = element_text(size = 7),
      legend.key.width  = unit(0.4, "cm"),
      legend.key.height = unit(1.5, "cm")
    )
}

n        <- nrow(plot_data)
half     <- ceiling(n / 2)

df_left  <- plot_data %>% 
  dplyr::slice(1:half) %>%
  mutate(gene_label = droplevels(gene_label))

df_right <- plot_data %>% 
  dplyr::slice((half + 1):n) %>%
  mutate(gene_label = droplevels(gene_label))

p_left  <- make_panel(df_left,  show_legend = FALSE)
p_right <- make_panel(df_right, show_legend = TRUE)


# Combine side by side
combined <- p_left + p_right + plot_layout(ncol = 2)
combined

ggsave(combined,
       filename = file.path(results_dir, "figures/rnaseq_introgressed_heatmap_20260723.pdf"),
       height = 180, width = 180, units = "mm")

ggsave(combined,
       filename = file.path(results_dir, "figures/rnaseq_introgressed_heatmap_20260723.png"),
       height = 140, width = 120, units = "mm")

# ============================================================================== #
# 11. FDR CORRECTION FOR SUPPLEMENTARY TABLE 7 ####

# Restricted FDR for the ecotype dataset, matching the approach used for the
# temperature dataset: BH applied within the introgressed genes actually tested.
tab <- read.csv(file.path(results_dir, "Cheng/introgressed_genes_expression_DE_summary.csv"), stringsAsFactors = FALSE)

# Exclude unannotated / novel entries
tab <- tab %>% filter(gene_class != "other_or_unannotated")

cat("Rows after removing unannotated entries:", nrow(tab), "\n")
# 120
cat("  Expressed          :", sum(tab$expression_status == "Expressed"), "\n")
# 92
cat("  Low_or_filtered    :", sum(tab$expression_status == "Low_or_filtered"), "\n")
# 28

# # Benjamini-Hochberg applied within the introgressed genes actually tested.
# NOTE: this must run AFTER filtering - the correction family is the set of
# genes reported in the table. Genes with NA p-values were not tested (either
# they failed the expression filter, or DESeq2 flagged a Cook's distance
# outlier) and are correctly excluded from the family
tested <- !is.na(tab$pvalue)

tab$FDR_introgressed <- NA_real_
tab$FDR_introgressed[tested] <- p.adjust(tab$pvalue[tested], method = "BH")

tab$DE_FDR05_introgressed <- ifelse(
  !is.na(tab$FDR_introgressed) & tab$FDR_introgressed < 0.05, "Yes", "No"
)

cat("\nGenes tested (correction family size):", sum(tested), "\n")
cat("Significant, genome-wide FDR < 0.05  :", sum(tab$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant, restricted FDR < 0.05   :", sum(tab$FDR_introgressed < 0.05, na.rm = TRUE), "\n")
# 89
# 0 
# 0


# Values quoted in the main text
cat("\nTHRB: log2FC = %.2f, P = %.4f, FDR(genome-wide) = %.3f, FDR(introgressed) = %.3f\n" %>%
      sprintf(
        tab$log2FoldChange_Baltic_vs_Atlantic[tab$gene_name == "THRB"],
        tab$pvalue[tab$gene_name == "THRB"],
        tab$padj[tab$gene_name == "THRB"],
        tab$FDR_introgressed[tab$gene_name == "THRB"]
      ))
# THRB: log2FC = 1.31, P = 0.0064, FDR(genome-wide) = 0.267, FDR(introgressed) = 0.283
cat("\nSEC16B: log2FC = %.2f, P = %.4f, FDR(genome-wide) = %.3f, FDR(introgressed) = %.3f\n" %>%
      sprintf(
        tab$log2FoldChange_Baltic_vs_Atlantic[tab$gene_name == "SEC16B"],
        tab$pvalue[tab$gene_name == "SEC16B"],
        tab$padj[tab$gene_name == "SEC16B"],
        tab$FDR_introgressed[tab$gene_name == "SEC16B"]
      ))

#SEC16B: log2FC = -1.16, P = 0.0686, FDR(genome-wide) = 0.605, FDR(introgressed) = 0.695

# Prepare table for publication: 
# ------------------------------------------------------------------ #
# 3. Select and rename columns to published headers
# ------------------------------------------------------------------ #
final <- tab %>%
  dplyr::select(
    `Old ENSEMBL ID`                        = old_ENSEMBLID,
    `New ENSEMBL ID`                        = new_ENSEMBLID,
    `Old Gene Name`                         = old_gene_name,
    `New Gene Name`                         = new_gene_name,
    `Introgression Region`                  = region,
    `Gene Class`                            = gene_class,
    `Expression Status`                     = expression_status,
    `Raw Total Count`                       = raw_total_count,
    `Raw Samples With Counts`               = raw_samples_with_count,
    `Mean Normalized Count`                 = mean_normalized_count,
    `Log2Fold Change Baltic vs Atlantic`    = log2FoldChange_Baltic_vs_Atlantic,
    `P Value`                               = pvalue,
    `FDR (genome-wide)`                     = padj,
    `DE at FDR (genome-wide) < 0.05`        = DE_FDR05,
    `FDR (introgressed genes)`              = FDR_introgressed,
    `DE at FDR (introgressed genes) < 0.05` = DE_FDR05_introgressed
  )

cat("\nFinal table:", nrow(final), "rows x", ncol(final), "columns\n")

final

write.csv(final, paste0(results_dir,"/Cheng/Supplementary_Table_7_with_restricted_FDR.csv"), row.names = FALSE)
#install.packages("writexl"); library(writexl)
writexl::write_xlsx(final, paste0(results_dir,"/Cheng/Supplementary_Table_7_with_restricted_FDR.xlsx"))

# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Tahoe 26.3.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Stockholm
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.3.2      cowplot_1.2.0        GenomicRanges_1.60.0 GenomeInfoDb_1.44.3  IRanges_2.42.0      
# [6] S4Vectors_0.46.0     BiocGenerics_0.54.1  generics_0.1.4       biomaRt_2.64.0       ggrepel_0.9.7       
# [11] lubridate_1.9.5      forcats_1.0.1        stringr_1.6.0        dplyr_1.2.0          purrr_1.2.1         
# [16] readr_2.2.0          tidyr_1.3.2          tibble_3.3.1         tidyverse_2.0.0      ggplot2_4.0.2       
# [21] edgeR_4.6.3          limma_3.64.3        
# 
# loaded via a namespace (and not attached):
#   [1] KEGGREST_1.48.1         gtable_0.3.6            httr2_1.2.2             Biobase_2.68.0         
# [5] lattice_0.22-9          tzdb_0.5.0              vctrs_0.7.1             tools_4.5.1            
# [9] curl_7.0.0              AnnotationDbi_1.70.0    RSQLite_2.4.6           blob_1.3.0             
# [13] pkgconfig_2.0.3         dbplyr_2.5.2            RColorBrewer_1.1-3      S7_0.2.1               
# [17] lifecycle_1.0.5         GenomeInfoDbData_1.2.14 compiler_4.5.1          farver_2.1.2           
# [21] textshaping_1.0.5       Biostrings_2.76.0       progress_1.2.3          statmod_1.5.1          
# [25] pillar_1.11.1           crayon_1.5.3            cachem_1.1.0            digest_0.6.39          
# [29] tidyselect_1.2.1        locfit_1.5-9.12         stringi_1.8.7           labeling_0.4.3         
# [33] splines_4.5.1           fastmap_1.2.0           grid_4.5.1              cli_3.6.5              
# [37] magrittr_2.0.4          utf8_1.2.6              withr_3.0.2             filelock_1.0.3         
# [41] rappdirs_0.3.4          prettyunits_1.2.0       scales_1.4.0            UCSC.utils_1.4.0       
# [45] bit64_4.6.0-1           timechange_0.4.0        XVector_0.48.0          httr_1.4.8             
# [49] bit_4.6.0               otel_0.2.0              ragg_1.5.1              png_0.1-8              
# [53] hms_1.1.4               memoise_2.0.1           BiocFileCache_2.16.2    rlang_1.1.7            
# [57] Rcpp_1.1.1              glue_1.8.0              DBI_1.3.0               xml2_1.5.2             
# [61] rstudioapi_0.18.0       jsonlite_2.0.0          R6_2.6.1                systemfonts_1.3.2     