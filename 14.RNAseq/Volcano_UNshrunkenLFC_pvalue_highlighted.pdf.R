###############################################
# Volcano plot only
# Baltic vs Atlantic
# - uses UNshrunken DESeq2 results
# - highlights candidate genes
# - handles gene IDs like: ENSCHAG00020004898|pcdh2ab3
# - caps extreme values at plot boundaries
###############################################

setwd("/Users/chema212/Desktop/DATA/Uppsala2023/Herring/v34_HE/b1_RNAseq/08_DEseq2/RefGtfAll/")

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

# --------------------------- 1) Input ---------------------------
data_list <- read.table("sample_table.csv", header = TRUE, sep = "\t",
                        row.names = 1, check.names = FALSE)
mydata    <- read.table("Ref_GTF_gene_count_matrix.csv", header = TRUE, sep = ",",
                        row.names = 1, check.names = FALSE)

# Align colData rows to countData columns exactly
stopifnot(all(colnames(mydata) %in% rownames(data_list)))
data_list <- data_list[colnames(mydata), , drop = FALSE]
stopifnot(identical(colnames(mydata), rownames(data_list)))

# --------------------------- 2) DESeq2 dataset ---------------------------
dds <- DESeqDataSetFromMatrix(countData = mydata,
                              colData   = data_list,
                              design    = ~ 1)

# --------------------------- 3) Subset to Stage == 12dpf and Baltic/Atlantic ---------------------------
keep <- !is.na(dds$Stage) & dds$Stage == "12dpf" &
        !is.na(dds$Population) & dds$Population %in% c("Atlantic", "Baltic")
dds <- dds[, keep]

for (v in c("Stage", "Condition", "Population", "Tissue")) {
  if (v %in% names(colData(dds))) {
    colData(dds)[[v]] <- droplevels(as.factor(colData(dds)[[v]]))
  }
}

# replicate checks
if (ncol(dds) < 4) stop("Too few samples after filtering (need >= 2 per group).")
tbl_pop <- table(dds$Population)
if (any(tbl_pop < 2)) {
  stop(sprintf("Need >= 2 replicates per population. Counts: %s",
               paste(names(tbl_pop), tbl_pop, collapse = ", ")))
}

# --------------------------- 4) Low-count filter ---------------------------
dds <- dds[rowSums(counts(dds) >= 100) >= 3, ]

# --------------------------- 5) Model ---------------------------
dds$Population <- relevel(factor(dds$Population), ref = "Atlantic")
design(dds) <- ~ Population

# --------------------------- 6) Fit DESeq2 ---------------------------
dds <- DESeq(dds, parallel = FALSE)

print(resultsNames(dds))
coef_name <- grep("Population_Baltic_vs_Atlantic", resultsNames(dds), value = TRUE)
if (length(coef_name) != 1) stop("Could not find coefficient 'Population_Baltic_vs_Atlantic'.")

res_raw <- results(dds, name = coef_name)
res_raw_df <- as.data.frame(res_raw)

# --------------------------- 7) Prepare volcano data ---------------------------
res_raw_df$Gene_full <- rownames(res_raw_df)

# Extract Ensembl gene ID before "|"
res_raw_df$Gene_ID <- sub("\\|.*", "", res_raw_df$Gene_full)

# Extract gene symbol after "|" if present
res_raw_df$Gene_symbol <- ifelse(grepl("\\|", res_raw_df$Gene_full),
                                 sub(".*\\|", "", res_raw_df$Gene_full),
                                 res_raw_df$Gene_full)

# Remove rows without needed values
volc_raw <- res_raw_df[complete.cases(res_raw_df[, c("log2FoldChange", "pvalue")]), ]

# Compute -log10(P)
eps <- .Machine$double.xmin
volc_raw$negLog10P <- -log10(pmax(volc_raw$pvalue, eps))

# --------------------------- 8) Candidate gene list ---------------------------
highlight_genes <- c(
"ENSCHAG00020031780",
"ENSCHAG00020031788",
"ENSCHAG00020031803",
"ENSCHAG00020031930",
"ENSCHAG00020032169",
"ENSCHAG00020031932",
"ENSCHAG00020031933",
"ENSCHAG00020031935",
"ENSCHAG00020031782",
"ENSCHAG00020000304",
"ENSCHAG00020000308",
"ENSCHAG00020000309",
"ENSCHAG00020000523",
"ENSCHAG00020000310",
"ENSCHAG00020000311",
"ENSCHAG00020000312",
"ENSCHAG00020000313",
"ENSCHAG00020000314",
"ENSCHAG00020000315",
"ENSCHAG00020000316",
"ENSCHAG00020000317",
"ENSCHAG00020000318",
"ENSCHAG00020000519",
"ENSCHAG00020000320",
"ENSCHAG00020011192",
"ENSCHAG00020025844",
"ENSCHAG00020012784",
"ENSCHAG00020012800",
"ENSCHAG00020012806",
"ENSCHAG00020012809",
"ENSCHAG00020012824",
"ENSCHAG00020016177",
"ENSCHAG00020025124",
"ENSCHAG00020025125",
"ENSCHAG00020017485",
"ENSCHAG00020017489",
"ENSCHAG00020019620",
"ENSCHAG00020021010",
"ENSCHAG00020021012",
"ENSCHAG00020021015",
"ENSCHAG00020021023",
"ENSCHAG00020021026",
"ENSCHAG00020024833",
"ENSCHAG00020024951",
"ENSCHAG00020024232",
"ENSCHAG00020024234",
"ENSCHAG00020017641",
"ENSCHAG00020017644",
"ENSCHAG00020017647",
"ENSCHAG00020008802",
"ENSCHAG00020010640",
"ENSCHAG00020010873",
"ENSCHAG00020010875",
"ENSCHAG00020010849",
"ENSCHAG00020010854",
"ENSCHAG00020010857",
"ENSCHAG00020025147",
"ENSCHAG00020010860",
"ENSCHAG00020010863",
"ENSCHAG00020025148",
"ENSCHAG00020010869",
"ENSCHAG00020010926",
"ENSCHAG00020011666",
"ENSCHAG00020024711",
"ENSCHAG00020024286",
"ENSCHAG00020025151",
"ENSCHAG00020024301",
"ENSCHAG00020024505",
"ENSCHAG00020024508",
"ENSCHAG00020024529",
"ENSCHAG00020024531",
"ENSCHAG00020024534",
"ENSCHAG00020024750",
"ENSCHAG00020024753",
"ENSCHAG00020025190",
"ENSCHAG00020025205",
"ENSCHAG00020024710",
"ENSCHAG00020024733",
"ENSCHAG00020021336",
"ENSCHAG00020021525",
"ENSCHAG00020021541",
"ENSCHAG00020024822",
"ENSCHAG00020025167",
"ENSCHAG00020024786",
"ENSCHAG00020015717",
"ENSCHAG00020015722",
"ENSCHAG00020015723",
"ENSCHAG00020015725",
"ENSCHAG00020016173",
"ENSCHAG00020031855",
"ENSCHAG00020016781",
"ENSCHAG00020016784",
"ENSCHAG00020017912",
"ENSCHAG00020027915",
"ENSCHAG00020027917",
"ENSCHAG00020012312",
"ENSCHAG00020012314",
"ENSCHAG00020012321",
"ENSCHAG00020013820",
"ENSCHAG00020013821",
"ENSCHAG00020017724",
"ENSCHAG00020017795",
"ENSCHAG00020018642",
"ENSCHAG00020018673",
"ENSCHAG00020022753",
"ENSCHAG00020024838",
"ENSCHAG00020024841",
"ENSCHAG00020028068",
"ENSCHAG00020024861",
"ENSCHAG00020024939",
"ENSCHAG00020027606",
"ENSCHAG00020027699",
"ENSCHAG00020027823",
"ENSCHAG00020027825",
"ENSCHAG00020027826",
"ENSCHAG00020027885",
"ENSCHAG00020027891",
"ENSCHAG00020027899",
"ENSCHAG00020026618",
"ENSCHAG00020027437",
"ENSCHAG00020024969",
"ENSCHAG00020025362",
"ENSCHAG00020011991",
"ENSCHAG00020012020",
"ENSCHAG00020014734",
"ENSCHAG00020014627",
"ENSCHAG00020008190",
"ENSCHAG00020010190",
"ENSCHAG00020013029",
"ENSCHAG00020013027",
"ENSCHAG00020013038",
"ENSCHAG00020013055",
"ENSCHAG00020007273",
"ENSCHAG00020008296",
"ENSCHAG00020011400",
"ENSCHAG00020011411",
"ENSCHAG00020011471",
"ENSCHAG00020003854",
"ENSCHAG00020003873"
)

# Mark highlighted genes
volc_raw$highlight <- ifelse(volc_raw$Gene_ID %in% highlight_genes, "Target", "Other")

# Label using gene symbol; if no symbol available, use Gene_ID
volc_raw$Label <- ifelse(
  volc_raw$highlight == "Target" &
  volc_raw$pvalue < 0.05 &
  abs(volc_raw$log2FoldChange) > 1,
  
  ifelse(is.na(volc_raw$Gene_symbol) | volc_raw$Gene_symbol == volc_raw$Gene_full,
         volc_raw$Gene_ID,
         volc_raw$Gene_symbol),
  
  ""
)


# Check overlap
matched_genes <- intersect(volc_raw$Gene_ID, highlight_genes)
cat("Matched candidate genes:", length(matched_genes), "\n")
print(matched_genes)

# --------------------------- 9) Define significance groups ---------------------------
volc_raw$regulation <- "NS"
volc_raw$regulation[volc_raw$pvalue < 0.05 & volc_raw$log2FoldChange >  1] <- "Up (Baltic)"
volc_raw$regulation[volc_raw$pvalue < 0.05 & volc_raw$log2FoldChange < -1] <- "Down (Baltic)"
volc_raw$regulation <- factor(volc_raw$regulation,
                              levels = c("Up (Baltic)", "Down (Baltic)", "NS"))

cols <- c(
  "Up (Baltic)"   = "#D62828",
  "Down (Baltic)" = "#1D4ED8",
  "NS"            = "#9CA3AF"
)

# --------------------------- 10) Cap values to plot boundaries ---------------------------
x_min <- -5
x_max <- 7
y_min <-  0
y_max <- 15

volc_raw$log2FC_capped    <- pmin(pmax(volc_raw$log2FoldChange, x_min), x_max)
volc_raw$negLog10P_capped <- pmin(pmax(volc_raw$negLog10P, y_min), y_max)

# Optional: mark capped points
volc_raw$is_capped <- with(volc_raw,
                           log2FoldChange < x_min |
                           log2FoldChange > x_max |
                           negLog10P > y_max)

# --------------------------- 11) Output directory ---------------------------
outdir <- "DEG_12dpf_Baltic_vs_Atlantic_volcano_only"
dir.create(outdir, showWarnings = FALSE)

# Save volcano input table
write.csv(volc_raw,
          file = file.path(outdir, "volcano_plot_input_with_highlight.csv"),
          row.names = FALSE)

# Save matched target genes only
write.csv(subset(volc_raw, highlight == "Target"),
          file = file.path(outdir, "highlighted_candidate_genes_in_volcano.csv"),
          row.names = FALSE)

# --------------------------- 12) Plot volcano ---------------------------
pdf(file.path(outdir, "Volcano_UNshrunkenLFC_pvalue_highlighted.pdf"),
    width = 6, height = 3.5)

ggplot(volc_raw, aes(x = log2FC_capped, y = negLog10P_capped)) +

  # background points
  geom_point(data = subset(volc_raw, highlight == "Other"),
             aes(color = regulation),
             size = 1.2, alpha = 0.75) +

  # highlighted target genes
  geom_point(data = subset(volc_raw, highlight == "Target"),
             shape = 21, fill = "yellow", color = "black",
             size = 2.3, stroke = 0.7) +

  # labels for target genes
  geom_text_repel(data = subset(volc_raw, highlight == "Target"),
                  aes(label = Label),
                  size = 3.2,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.alpha = 0.6,
                  min.segment.length = 0) +

  scale_color_manual(values = cols,
                     name = "Regulation\n(p < 0.05; |LFC| >= 1)") +

  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +

  scale_x_continuous(limits = c(x_min, x_max),
                     breaks = seq(x_min, x_max, by = 1)) +
  scale_y_continuous(limits = c(y_min, y_max),
                     breaks = seq(y_min, y_max, by = 5)) +

  labs(x = "log2FC Baltic vs Atlantic",
       y = "-log10(P)",
       title = "Volcano plot: Baltic vs Atlantic at 12dpf") +

  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5))

dev.off()

cat("DONE\n")
cat("Output folder:", outdir, "\n")


