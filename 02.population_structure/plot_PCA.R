# PCA Plots 
# 2025-08-11
# LOGS ####

# Libraries ####
library(ggrastr)
library(tidyverse)
library(cowplot)
library(viridis)
library(biomaRt)

# STWD ####
# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
# REPO is resolved to an absolute path before we setwd() below, so it stays
# valid for any later reads/writes regardless of the working directory.
REPO <- normalizePath(".")
setwd("../Figshare/2.population_structure/PCA")


# I asked chat to simplify the code:

library(tidyverse)
library(cowplot)
library(gridExtra)

#-------------------------
# Sample info (used for all plots)
#-------------------------
# indTable <- read.table(
#   "../plotting_files/sampleinfo125.txt.ordered.mafalda",
#   header = TRUE, sep = "\t", comment.char = ""
# )

# Grab new order file:
indTable <- read.table(
  file.path(REPO, "00.revisions_relabeling/sampleinfo125_consolidated.txt"),
  header = TRUE, sep = "\t", comment.char = ""
)


#-------------------------
# File paths (your actual files)
#-------------------------
EVEC_1 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.1SNPevery10000.eigenvec"
EVAL_1 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.1SNPevery10000.eigenval"

EVEC_2 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela.eigenvec"
EVAL_2 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela.eigenval"

EVEC_3 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv.eigenvec"
EVAL_3 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv.eigenval"


EVEC_4 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.1SNPevery10000.eigenvec"
EVAL_4 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.1SNPevery10000.eigenval"

EVEC_5 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela.eigenvec"
EVAL_5 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela.eigenval"

EVEC_6 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela_woInv.eigenvec"
EVAL_6 <- "data/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela_woInv.eigenval"



# Helper: fail early with a helpful message if a file is missing
assert_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, "\nCurrent working directory: ", getwd(), call. = FALSE)
  }
}
walk(c(EVEC_1,EVAL_1,EVEC_2,EVAL_2,EVEC_3,EVAL_3), assert_file)

# ---- Pretty-print population labels (no underscores, read like sentences) --
# General rule: replace "_" with " ". A couple of combined-group labels read
# better as full phrases than as space-separated words, so they're
# special-cased here instead.

label_overrides <- c(
  "Canada_Spring_Summer"     = "Canada Spring and Summer",
  "Britain_Ireland_NorthSea" = "Britain, Ireland and North Sea"
)

clean_label <- function(x) {
  if (x %in% names(label_overrides)) return(label_overrides[[x]])
  gsub("_", " ", x)
}

#-------------------------
# Function to read PCA and return a ggplot (or just its legend)
#-------------------------
make_pca_plot <- function(evec_file, eval_file, title = "", return_legend = FALSE) {
  pca <- readr::read_table(evec_file, col_names = FALSE)
  eigenval <- scan(eval_file)
  
  # Drop nuisance first column, set names to match your code
  pca <- pca[ , -1]
  names(pca)[1] <- "ind"
  names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
  
  # % variance explained
  pve <- eigenval / sum(eigenval) * 100
  
  # Merge with metadata + ensure stable factor order for colors
  pcaWithInfo <- merge(pca, indTable, by.x = "ind", by.y = "id")
  ordered <- pcaWithInfo[order(pcaWithInfo$order_in_admixture), ]
  ordered$bars_plot <- factor(ordered$bars_plot, levels = unique(ordered$bars_plot))
  # Relabel the factor's levels in place -- order/count untouched, so this
  # doesn't affect the color matching above, just the displayed text.
  levels(ordered$bars_plot) <- vapply(levels(ordered$bars_plot), clean_label, character(1))
  
  p <- ggplot(ordered, aes(PC1, PC2, col = bars_plot, shape = species)) +
    geom_point(size = 2, alpha = 0.9) +
    scale_color_manual(values = unique(ordered$color)) +
    theme_classic(base_size = 9) +
    labs(
      x = paste0("PC1 (", signif(pve[1], 3), "%)"),
      y = paste0("PC2 (", signif(pve[2], 3), "%)"),
      title = title
    ) +
    theme(
      legend.position = if (return_legend) "bottom" else "none",
      legend.box = "horizontal",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 10)
    )
  
  if (return_legend) {
    p <- p + guides(
      shape = guide_legend(title = NULL, order = 1, ncol = 1, override.aes = list(size = 3)),
      color = guide_legend(title = NULL, order = 2, ncol = 2, override.aes = list(size = 3))
    )
    return(cowplot::get_legend(p))
  }
  
  p
}

#-------------------------
# Build the three PCA plots
#-------------------------
windows_pca <- make_pca_plot(EVEC_1, EVAL_1)
top50snpsinversion <- make_pca_plot(EVEC_2, EVAL_2)
top50snpswoinversion <- make_pca_plot(EVEC_3, EVAL_3)

windows_pca_noPac <- make_pca_plot(EVEC_4, EVAL_4)
top50snpsinversion_noPac <- make_pca_plot(EVEC_5, EVAL_5)
top50snpswoinversion_noPac <- make_pca_plot(EVEC_6, EVAL_6)


#-------------------------
# Arrange: 3 plots + legend as the 4th column
#-------------------------
# Legend (extract once from any plot with legend shown)
# legend_grob <- cowplot::get_legend(
#   make_pca_plot(EVEC_1, EVAL_1, show_legend = TRUE) +
#     guides(shape = guide_legend(override.aes = list(size = 3)))
# )
# 

# all_pcas <- gridExtra::grid.arrange(
#   windows_pca, top50snpsinversion, top50snpswoinversion, legend_grob,
#   ncol = 4,
#   widths = c(1, 1, 1, 0.6)
# )

# ggsave(all_pcas, filename="figures/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.allPCAs.pdf",
#        units = "mm",
#        height = 90,
#        width = 210)




# # Legend (extract once from any plot with legend shown)
# legend_grob <- cowplot::get_legend(
#   make_pca_plot(EVEC_4, EVAL_4, show_legend = TRUE) +
#     guides(shape = guide_legend(override.aes = list(size = 3)))
# )
# 
# all_pcas <- gridExtra::grid.arrange(
#   windows_pca_noPac, top50snpsinversion_noPac, top50snpswoinversion_noPac, legend_grob,
#   ncol = 4,
#   widths = c(1, 1, 1, 0.6)
# )
# 
# 
# ggsave(all_pcas, filename="figures/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.allPCAs_JustAtlanticherring.pdf",
#        units = "mm",
#        height = 90,
#        width = 210)


# ============================================================================
# Revisions -- 2026-08-13
# Revising figure labels for consistency across figures: combining the two
# PCA rows (top = all individuals incl. Pacific herring, bottom = Atlantic
# only) into one figure with a single shared legend (from the top row) and
# two-level column headers.
# ============================================================================

# ---- Column headers ---------------------------------------------------
# Top level: "1 SNP every 10 kb" over column 1 only; "Top 50 SNPs in regions
# under selection" spanning columns 2-3.
# Second level: blank under column 1; "With inversions" / "Without
# inversions" under columns 2-3.
header_top_1       <- grid::textGrob("1 SNP every 10 kb", gp = grid::gpar(fontsize = 8))
header_top_2       <- grid::textGrob("Top 50 SNPs in regions under selection", gp = grid::gpar(fontsize = 8))
header_sub_blank   <- grid::nullGrob()
header_sub_with    <- grid::textGrob("With inversions", gp = grid::gpar(fontsize = 8))
header_sub_without <- grid::textGrob("Without inversions", gp = grid::gpar(fontsize = 8))

# ---- Single shared legend, taken from the top row, laid out horizontally --
legend_grob <- make_pca_plot(EVEC_1, EVAL_1, return_legend = TRUE)

# ---- Panel labels "a" / "b" ------------------------------------------
panel_label_a <- grid::textGrob("a", gp = grid::gpar(fontsize = 9, fontface = "bold"))
panel_label_b <- grid::textGrob("b", gp = grid::gpar(fontsize = 9, fontface = "bold"))

grobs_list <- list(
  header_top_1, header_top_2,                                            # 1, 2
  header_sub_blank, header_sub_with, header_sub_without,                 # 3, 4, 5
  legend_grob,                                                           # 6
  windows_pca, top50snpsinversion, top50snpswoinversion,                 # 7, 8, 9
  windows_pca_noPac, top50snpsinversion_noPac, top50snpswoinversion_noPac, # 10, 11, 12
  header_top_1, header_top_2,                                            # 13, 14 (repeat of 1, 2)
  header_sub_blank, header_sub_with, header_sub_without,                 # 15, 16, 17 (repeat of 3, 4, 5)
  panel_label_a, panel_label_b                                           # 18, 19
)

layout_matrix <- rbind(
  c(18, 1,  2,  2),
  c(NA, 3,  4,  5),
  c(NA, 7,  8,  9),
  c(19, 13, 14, 14),
  c(NA, 15, 16, 17),
  c(NA, 10, 11, 12),
  c(NA, 6,  6,  6)
)

combined_pca <- gridExtra::grid.arrange(
  grobs = grobs_list,
  layout_matrix = layout_matrix,
  widths  = c(0.08, 1, 1, 1),
  heights = c(0.15, 0.15, 1, 0.15, 0.15, 1, 0.8)
)

ggsave(combined_pca,
       filename = file.path(REPO, "02.population_structure/figures/herring_sentieon_125ind_231031_PCA_combined_revised_20260813.png"),
       units = "mm",
       width = 180,
       height = 210)
