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
setwd("Manuscript/Figshare/2.population_structure/PCA")


# I asked chat to simplify the code:

library(tidyverse)
library(cowplot)
library(gridExtra)

#-------------------------
# Sample info (used for all plots)
#-------------------------
indTable <- read.table(
  "../plotting_files/sampleinfo125.txt.ordered.mafalda",
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

#-------------------------
# Function to read PCA and return a ggplot
#-------------------------
make_pca_plot <- function(evec_file, eval_file, title = "", show_legend = FALSE) {
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
  
  ggplot(ordered, aes(PC1, PC2, col = bars_plot, shape = species)) +
    geom_point(size = 2, alpha = 0.9) +
    scale_color_manual(values = unique(ordered$color)) +
    theme_classic(base_size = 9) +
    labs(
      x = paste0("PC1 (", signif(pve[1], 3), "%)"),
      y = paste0("PC2 (", signif(pve[2], 3), "%)"),
      title = title
    ) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 10)
    )
}

#-------------------------
# Build the three PCA plots
#-------------------------
windows_pca <- make_pca_plot(EVEC_1, EVAL_1, "1 SNP every 10k")
top50snpsinversion <- make_pca_plot(EVEC_2, EVAL_2, "Top 50 SNPs (with inversions)")
top50snpswoinversion <- make_pca_plot(EVEC_3, EVAL_3, "Top 50 SNPs (without inversions)")


windows_pca_noPac <- make_pca_plot(EVEC_4, EVAL_4, "1 SNP every 10k")
top50snpsinversion_noPac <- make_pca_plot(EVEC_5, EVAL_5, "Top 50 SNPs (with inversions)")
top50snpswoinversion_noPac <- make_pca_plot(EVEC_6, EVAL_6, "Top 50 SNPs (without inversions)")



#-------------------------
# Arrange: 3 plots + legend as the 4th column
#-------------------------
all_pcas <- gridExtra::grid.arrange(
  windows_pca, top50snpsinversion, top50snpswoinversion, legend_grob,
  ncol = 4,
  widths = c(1, 1, 1, 0.6)
)

ggsave(all_pcas, filename="figures/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.allPCAs.pdf",
       units = "mm",
       height = 90,
       width = 210)




# Legend (extract once from any plot with legend shown)
legend_grob <- cowplot::get_legend(
  make_pca_plot(EVEC_4, EVAL_4, show_legend = TRUE) +
    guides(shape = guide_legend(override.aes = list(size = 3)))
)

all_pcas <- gridExtra::grid.arrange(
  windows_pca_noPac, top50snpsinversion_noPac, top50snpswoinversion_noPac, legend_grob,
  ncol = 4,
  widths = c(1, 1, 1, 0.6)
)


ggsave(all_pcas, filename="figures/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.allPCAs_JustAtlanticherring.pdf",
       units = "mm",
       height = 90,
       width = 210)

