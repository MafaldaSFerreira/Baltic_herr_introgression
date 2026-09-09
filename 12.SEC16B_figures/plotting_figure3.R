# This R code allows to re-plot Figure 3 in the Baltic introgression paper
# All input files can be found either on GitHub or Figshare (if too big)

setwd("Manuscript/Figshare/")
setwd("~/Documents/Postdoc/Project_Herring/Introgression/Manuscript/Figshare/")

library(ggrastr)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggrepel)

# INPUT FILES ####

## TWISST ####

weights <- read.table("6.twisst/outputs/output_outgroup_sprat/output.wg.phyml.w100.Mi10.weights.csv.gz", skip = 3, header=T)
windows <- read.table("6.twisst/outputs/output_outgroup_sprat/output.wg.phyml.w100.Mi10.data.tsv", header=T)

weights = weights / apply(weights, 1, sum)
good_rows = which(is.na(apply(weights, 1, sum)) == F)
weights <- weights[good_rows,]
windows <- windows[good_rows,]

df<-cbind(windows,weights)
df<-mutate(df,
           topo1_frac=topo1/(topo1+topo2+topo3),
           topo2_frac=topo2/(topo1+topo2+topo3),
           topo3_frac=topo3/(topo1+topo2+topo3))

## ALLELE FREQUENCIES ####

load("5.processing_pool_seq_data/inputs/60.Neff.AF.2024-08-14.Rdata")
names<-names(pops_freq_df)
pool_order <- read.table("5.processing_pool_seq_data/plotting_files/pool_order")

## CHI-SQUARE ####
load("5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_baltic_autumn_han_pops_chiseq.output.RData")
load("5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData")

## MUTATIONS ####
moderate_mutations_50kb <- read.table(sep="\t", header = T, "5.processing_pool_seq_data/results/moderate_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")

## INTROGRESSION REGIONS ####
intro_reg<-read.table(header=T,"4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

## SPRAT GENOTYPES ####
sprat.genotype<-read.table(sep="\t", header=F, "1.mapping_variant_calling/genotypes/all_mutations_combined.sprat.genotypes.txt")

## SELECTION RESULTS ####
xpehh_results<-read.table("8.xpehh/results_xpehh/xpehh_group1+2.out", skip=1)
xpehh_results$CHROM<-str_split_fixed(xpehh_results$V2, "_", 2)[,1]
xpehh_results$POS<-as.numeric(str_split_fixed(xpehh_results$V2, "_", 2)[,2])
xpehh_results$V2<-NULL

colnames(xpehh_results)<-c("Index", "Freq", "iHH_A1", "iHH_B1", "iHH_P1", "XPEHH", "std_XPEHH", "CHROM", "POS")

xpehh_results_sig<-xpehh_results  %>% 
  arrange(desc(abs(std_XPEHH))) %>%   # Sort SNPs by XPEHH score (descending)
  mutate(
    rank=row_number(),           # number of SNPs with higher scores
    rank.perc=rank/n(),          # fraction of SNPs with higher score
    rank.log=-log10(rank.perc)   # P-value
  )

## ANNOTATIONS ####
scan1_match_df<-read.table("4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.modified.txt", sep="\t", header=T, row.names=NULL)

## REVISIONS: FST SCAN ####
#load("5.processing_pool_seq_data/fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst")
load("5.processing_pool_seq_data/fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData")


# PAPER PLOTS #####

# Define regoin to plot
padding <- 0.15
padding_mb <- 1.5e5

c<-"chr10"
s<-25040001

#s<-25000001 # trying this to extend the Twisst plot slightly

e<-25220000
chr<-str_remove(c, "chr")

# Subset the input files 
# Positive selection
tmp_xpehh<-xpehh_results_sig %>% filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

# ChiSQ
#tmp_csq_BS_AS <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

# Revision:
# Join FST values to the chi-square object by position
tmp_csq_BS_AS <- baltic_spring_vs_atlantic_spring_chiseq %>% 
  filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb) %>%
  left_join(df_spring_maf0.05 %>% filter(FST >= 0) %>% select(Chr, Pos, FST), 
            by=c("CHROM"="Chr", "POS"="Pos"))


# Allele Frequencies 
tmp_pops_freq_df <- pops_freq_df %>% 
  filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb) %>%
  pivot_longer(cols = 3:62, names_to = "populations", values_to = "frequencies")

# Define the introgression region:
candidate_matrix <- intro_reg %>% filter(seqnames == c & start >= s-padding_mb & end <= e+ padding_mb)

# Positive selection 
plot_xpehh_2<-ggplot(tmp_xpehh)+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=6,ymax=-6),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_point(aes(x=POS/1e6, y=std_XPEHH, color=factor(abs(std_XPEHH) >= 2)), size=0.5)+
  scale_color_manual(values=c("black", "red"))+
  ylim(-6,6)+
  xlab(paste0(chr," position (Mb)"))+
  ylab("XPEHH score")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_blank())

# Chi-square
# plot_csq_BS_AS <- ggplot(tmp_csq_BS_AS)+
#   geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
#             fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
#   geom_point(aes(x=POS/1e6, y=-log10(chisq_p), color=factor(-log10(chisq_p) > 9.1)), size=0.5)+
#   scale_color_manual(values=c("gray", "black"))+
#   #geom_line(aes(x=POS/1e6, y=rollmean(dAF_BS_vs_AS, 15, na.pad=TRUE))) +
#   theme_classic()+
#   xlim(s/1e6-padding, e/1e6 + padding)+
#   xlab(paste0(chr," position (Mb)"))+
#   ylab(expression(paste(-log[10], " (P value)") ) )+
#   theme(legend.position = "none",
#         axis.text.y = element_text(colour="black",size=6),
#         axis.text.x = element_text(colour="black",size=6),
#         axis.title.y = element_text(colour="black",size=6),
#         axis.title.x = element_blank())


# Plot dAF coloured by FST
plot_csq_BS_AS <- ggplot(tmp_csq_BS_AS) +
  geom_rect(data=candidate_matrix, mapping=aes(xmin=start/1e6, xmax=end/1e6, ymin=0, ymax=80),
            fill="azure3", color=NA, alpha=0.5, inherit.aes=FALSE) +
  # Grey points first so they go to the back
  geom_point(data=tmp_csq_BS_AS %>% filter(is.na(FST)),
             aes(x=POS/1e6, y=-log10(chisq_p)), 
             color="darkgreen", size=0.5) +
  # Coloured points on top
  geom_point(data=tmp_csq_BS_AS %>% filter(!is.na(FST)),
             aes(x=POS/1e6, y=-log10(chisq_p), color=FST), size=0.5) +
  scale_color_viridis(option="magma", direction=-1,
                      breaks=c(0, 0.1, 0.2, 0.25),
                      labels=c("0", "0.1", "0.2", ">0.25"),
                      limits=c(0, 0.25), oob=scales::squish,
                      name=expression(F[ST])) +
  theme_classic() +
  xlim(s/1e6-padding, e/1e6+padding) +
  xlab(paste0(chr, " position (Mb)")) +
  ylab(expression(paste(-log[10], " (P value)"))) +
  theme(legend.position="right",
        legend.key.size=unit(0.3, "cm"),
        legend.title=element_text(size=6),
        legend.text=element_text(size=5),
        axis.text.y=element_text(colour="black", size=6),
        axis.text.x=element_text(colour="black", size=6),
        axis.title.y=element_text(colour="black", size=6),
        axis.title.x=element_blank())

# Allele Frequency heat map
plot_pops_freq_df <- tmp_pops_freq_df %>% 
  ggplot()+
  rasterise(geom_tile(aes(y=factor(populations, level=pool_order$V1), x=POS/1e6, color=frequencies)), dpi=300) + 
  scale_color_viridis(direction = -1) +
  ylab("Populations")+
  xlab("Allele frequencies")+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        axis.title.y = element_text(size=6),
        legend.position = "none")

# Gene models
numrow<-scan1_match_df %>% 
  filter(seqnames==chr) %>% summarise(numrow=n())

# yy for gene names
yy<-rep(seq(1,4),(numrow$numrow/4))

if(length(yy) < numrow$numrow) {
  
  extra<-numrow$numrow-length(yy)
  
  if (extra == 1) { yy <- c(yy, 1)
  } else if (extra == 2) { yy <- c(yy, 1, 2)
  } else if (extra == 3) { yy <- c(yy, 1, 2, 3)} 
  
  
}

gene_names<-scan1_match_df %>% 
  filter(seqnames==chr) %>%
  ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_rect(aes(xmin=as.numeric(start)/1e6,xmax=as.numeric(end)/1e6
                ,ymin=1,ymax=1.25),colour="black",fill="lightblue",linewidth=0.25)+
  geom_text(aes(label=external_gene_name, x=start/1e6, y=2, angle=-45), size=2.5)+
  ylim(0,3)+
  xlim(s/1e6-padding, e/1e6+padding)+
  xlab(paste0(chr," position (Mb)"))+
  ylab("gene models")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=6),
        legend.position = "none")

# Twisst plot
twisst_plot <- df %>%
  filter(scaffold==c) %>%
  ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=topo1),col="#882E72",size=0.5)+
  geom_line(aes(x=mid/1000000,y=topo2),col="#E8601C",size=0.5)+
  geom_line(aes(x=mid/1000000,y=topo3),col="#5289C7",size=0.5)+
  xlim(s/1e6-padding, e/1e6+padding)+
  ylab("topology support (%)") +
  xlab(paste0(chr," position (Mb)"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "none",
        axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_text(colour="black",size=6))

# Add all plots together
all_plots <- plot_grid(gene_names, plot_csq_BS_AS, plot_pops_freq_df, plot_xpehh_2, twisst_plot, nrow=5, align = "v")


# Save plots
#ggsave(all_plots, filename=paste0("sec16b_panel_2025-08-09.pdf"), height = 110, width = 100, units = "mm")
#ggsave(all_plots, filename=paste0("sec16b_panel_2026-02-25.pdf"), height = 110, width = 100, units = "mm")
ggsave(all_plots, filename=paste0("~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure3/sec16b_panel_2026-05-20.pdf"), height = 110, width = 100, units = "mm")
ggsave(all_plots, filename=paste0("~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure3/sec16b_panel_2026-05-20_legend.pdf"), height = 110, width = 100, units = "mm")

# Plot heatmap independently and bigger
plot_pops_freq_df <- tmp_pops_freq_df %>% 
  ggplot()+
  rasterise(geom_tile(aes(y=factor(populations, level=rev(pool_order$V1)), x=POS/1e6, color=frequencies)), dpi=300) + 
  scale_color_viridis(direction = -1) +
  ylab("Populations")+
  xlab("Allele frequencies")+
  theme_classic()+
  theme(#axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6),
    axis.title.y = element_text(size=6),
    legend.position = "none")

ggsave(plot_pops_freq_df, filename=paste0("sec16b_AF_heatmap.pdf"), height = 50, width = 100, units = "mm")


# GENOTYPE HEAT MAP ####

# Select the missense mutations:
sec16B_snps <- moderate_mutations_50kb %>% filter(seqnames=="chr10") %>% dplyr::select(start)

# Read in the genotype information:
genotypes <- read.table("1.mapping_variant_calling/genotypes/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.intro_reg.genotypes.2023-11-30.txt")
indv_v02<-read.table("1.mapping_variant_calling/genotypes/individuals_v01_20231116.txt")
header<-read.table("1.mapping_variant_calling/genotypes/header_vcf_20231116.txt")

# Pacific and Atlantic herring genotypes
colnames(genotypes)<-header[,1]

genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/0"]<-1
genotypes[genotypes=="1/1"]<-2
genotypes[genotypes=="./."]<-NA

# Merge with Sprat genotypes
sprat.genotype$V4 <- NULL
sprat.genotype[sprat.genotype=="0/0"]<-0
sprat.genotype[sprat.genotype=="0/1"]<-1
sprat.genotype[sprat.genotype=="1/0"]<-1
sprat.genotype[sprat.genotype=="1/1"]<-2
sprat.genotype[sprat.genotype=="./."]<-NA
colnames(sprat.genotype) <- c("CHROM", "POS", "SPRAT")

HerrSprat_genotypes <- left_join(genotypes, sprat.genotype, by=c("CHROM", "POS"))

# Subset genotypes
sec16b_missense_Genos_freq <- HerrSprat_genotypes %>% 
  filter(CHROM=="chr10") %>% filter(POS %in% sec16B_snps$start) %>%
  filter(POS >=25040001 & POS <=25220000) %>%
  pivot_longer(cols=3:128, names_to="Individuals", values_to="Genotype")

# Remove Hybrid population Balsfjord
sec16b_missense_Genos_freq_filter <- sec16b_missense_Genos_freq %>% 
  filter(!Individuals %in% c("HWS61_Balsfjord_Atlantic", "HWS62_Balsfjord_Atlantic",
                                                                                               "HWS63_Balsfjord_Atlantic", "HWS64_Balsfjord_Atlantic"))

indv_v03 <- indv_v02 %>% filter(!V1 %in% c("HWS61_Balsfjord_Atlantic", "HWS62_Balsfjord_Atlantic",
                                           "HWS63_Balsfjord_Atlantic", "HWS64_Balsfjord_Atlantic"))

plot_genotypes <- sec16b_missense_Genos_freq_filter %>%
  ggplot()+
  rasterise(geom_tile(aes(y=factor(Individuals,level=indv_v03[,1]), x=as.character(POS), fill=Genotype)), dpi=300) + 
  #scale_fill_manual(values=c("#1F78B4", "#CCCCCC", "#D95F02", "black"))+
  scale_fill_viridis(discrete = TRUE) +
  ylab("Populations")+
  xlab("Genotypes")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=6),
    axis.text.x = element_text(size=6, angle=90, hjust=0.5, vjust=0.5),
    axis.title.y = element_blank(),
    legend.position = "none")

ggsave(plot_genotypes, filename="sec16B_genotypes_big_2025-09-13.pdf", height = 100, width = 190, units = "mm")


# ============================================================================= #
# FIGURE 3 -- REVISIONS ADDENDUM (2026-08-26) ####
# Append to the bottom of the existing Figure 3 script. Nothing above is edited;
# this re-derives the panels that needed fixing and re-assembles the figure.
#
# What this fixes:
#   1. Population group labels on both heatmaps, generated from the data
#      instead of added by hand in Affinity.
#   2. Reads the CORRECTED pool_order, not the stale live one.
#   3. AF heatmap sized via rel_heights + its own frequency legend on the right.
#   4. Gene name labels placed by ggrepel instead of a fixed -45 deg rotation.
#   5. Genotype heatmap folded into the figure, with a gene header above it,
#      a species strip beside it, and a four-key genotype legend below.
# ============================================================================= #

library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggrastr)
library(ggrepel)

REPO <- "~/Documents/Postdoc/Repositories/Baltic_herr_introgression"

# ----------------------------------------------------------------------------- #
# 0. Corrected pool_order + display labels
# ----------------------------------------------------------------------------- #
# The live 5.processing_pool_seq_data/plotting_files/pool_order is STILL the old
# version: Orkney sits in Ireland_Britain, N_NorthSea in NEAtlantic, and there is
# no BrackishSemilandlocked group. Using it would re-introduce the mislabels this
# revision is meant to fix.
pool_order_corr <- read.table(
  file.path(REPO, "00.revisions_relabeling/pool_order_20260813.txt"),
  header = FALSE, stringsAsFactors = FALSE, comment.char = "", sep = "\t")
colnames(pool_order_corr) <- c("population", "raw_group")

pool_group_display <- c(
  "PacificOceanPacificherring"   = "Pacific herring",
  "ArcticPacificherring"         = "Arctic Pacific herring",
  "BalsfjordHybridPopulation"    = "Balsfjord",
  "BalticSummer"                 = "Baltic Summer",
  "BalticSpring"                 = "Baltic Spring",
  "BalticAutum"                  = "Baltic Autumn",      # raw code really is "Autum"
  "BrackishSemilandlocked"       = "Semi-landlocked\nbrackish",
  "BalticAtlanticTransitionZone" = "Baltic-Atlantic\ntransition",
  "NEAtlanticFjords"             = "NE Atlantic fjords",
  "NEAtlantic"                   = "NE Atlantic",
  "NorthSea"                     = "North Sea",
  "Ireland_Britain"              = "Britain and Ireland",
  "NWAtlantic"                   = "NW Atlantic"
)
stopifnot(all(pool_order_corr$raw_group %in% names(pool_group_display)))

# ----------------------------------------------------------------------------- #
# 1. Helper: plotting order + group assignment -> axis breaks/labels
# ----------------------------------------------------------------------------- #
group_runs <- function(row_levels, row_groups) {
  r      <- rle(as.character(row_groups))
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  mids   <- floor((starts + ends) / 2)
  
  if (anyDuplicated(r$values)) {
    dup <- unique(r$values[duplicated(r$values)])
    warning("These groups are NOT contiguous in the plotting order, so each ",
            "appears as more than one label: ", paste(dup, collapse = ", "),
            call. = FALSE)
  }
  
  data.frame(group    = r$values,
             start_i  = starts,
             end_i    = ends,
             mid_row  = row_levels[mids],
             boundary = ends + 0.5,
             stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------------------- #
# 2. Allele-frequency heatmap, with group labels AND a frequency legend
# ----------------------------------------------------------------------------- #
# ORIENTATION: levels[1] is drawn at the BOTTOM, so rev() puts pool_order row 1
# (Pacific herring) at the TOP. The original script was inconsistent about this
# -- the panel in all_plots used level = pool_order$V1, the standalone export
# used rev(). Those are vertically mirrored.
pool_order_corr <- pool_order_corr %>%
  mutate(raw_group = factor(raw_group, levels = unique(raw_group))) %>%
  arrange(raw_group) %>%
  mutate(raw_group = as.character(raw_group))

af_levels <- rev(pool_order_corr$population)
af_groups <- rev(pool_order_corr$raw_group)
runs_af   <- group_runs(af_levels, af_groups)

plot_pops_freq_df_lab <- tmp_pops_freq_df %>%
  ggplot() +
  rasterise(geom_tile(aes(y = factor(populations, levels = af_levels),
                          x = POS / 1e6, color = frequencies)), dpi = 300) +
  # NOTE: white separators are hard to see against the pale end of the viridis
  # ramp. Change to "black" if the group boundaries don't read clearly.
  geom_hline(yintercept = head(runs_af$boundary, -1), colour = "white", linewidth = 0.3) +
  # legend on the right. Because every panel in the grid is aligned on its PANEL
  # area (axis = "lr"), the legend sits outside that area and the other panels
  # are simply padded to match -- the x-axes stay aligned.
  scale_color_viridis(direction = -1, name = "Allele\nfrequency",
                      breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_y_discrete(breaks = runs_af$mid_row,
                   labels = unname(pool_group_display[runs_af$group]),
                   expand = c(0, 0)) +
  ylab(NULL) +
  theme_classic() +
  theme(axis.title.x        = element_blank(),
        axis.text.x         = element_text(size = 6),
        axis.text.y         = element_text(size = 5, colour = "black", lineheight = 0.8),
        axis.ticks.y        = element_line(colour = "grey40", linewidth = 0.2),
        axis.ticks.length.y = unit(2.5, "mm"),
        axis.line.y         = element_blank(),
        legend.position     = "right",
        legend.title        = element_text(size = 5),
        legend.text         = element_text(size = 4),
        legend.key.width    = unit(2, "mm"),
        legend.key.height   = unit(3, "mm"))

# ----------------------------------------------------------------------------- #
# 3. Genotype heatmap: population labels, species strip, genotype legend
# ----------------------------------------------------------------------------- #
# Individual names in the genotype matrix are `id` values (the VCF was
# re-headered -- note "newID" in its filename), NOT `name_in_vcf`. 49 of the 125
# differ between those two columns.
indTable_f3 <- read.table(
  file.path(REPO, "00.revisions_relabeling/sampleinfo125_consolidated.txt"),
  header = TRUE, sep = "\t", comment.char = "")   # comment.char="" -> hex colours survive

# Deliberately coarser than bars_plot on the Pacific side, matching the grouping
# the pool-seq allele-frequency heatmaps use. Hardcoded on purpose: this is a
# presentation choice for these figures, not a property of the sample table.
#
# Both spellings are covered because bars_plot carries the Figshare-side labels
# unless 00.population_labels.R Step 7 has overwritten them with correct_label.
GENOTYPE_GROUP_DISPLAY <- c(
  "Vancouver"                = "Pacific Pacific herring",
  "Japan"                    = "Pacific Pacific herring",
  "Barents_Sea"              = "Arctic Pacific herring",
  "White_Sea"                = "Arctic Pacific herring",
  "Pechora_Sea"              = "Arctic Pacific herring",
  "Balsfjord"                = "Balsfjord hybrid\npopulation",
  "Canada_Autumn"            = "Canada Autumn",
  "Canada_Spring_Summer"     = "Canada Spring\nand Summer",
  "Canada_Spring"            = "Canada Spring\nand Summer",
  "Canada_Summer"            = "Canada Spring\nand Summer",
  "Britain_Ireland_NorthSea" = "Britain, Ireland\nand North Sea",
  "Ireland_Britain"          = "Britain, Ireland\nand North Sea",
  "North_Sea"                = "Britain, Ireland\nand North Sea",
  "Norway_Coastal"           = "Norway stationary",
  "Norway_Costal"            = "Norway stationary",
  "Norway_stationary"        = "Norway stationary",
  "Norway_Pelagic"           = "Norway migratory",
  "Norway_migratory"         = "Norway migratory",
  "Baltic_Sea_Autumn"        = "Baltic Autumn",
  "Baltic_Autumn"            = "Baltic Autumn",
  "Baltic_Sea_Spring"        = "Baltic Spring",
  "Baltic_Spring"            = "Baltic Spring",
  "SPRAT"                    = "European sprat"
)

species_display <- c("C.harengus" = "Atlantic\nherring",
                     "C.pallasii" = "Pacific\nherring",
                     "SPRAT"      = "European\nsprat")

geno_levels <- as.character(indv_v03[, 1])

# population AND species come from the same join, so the y labels and the
# species strip can never disagree about which individual sits where
geno_meta <- data.frame(id = geno_levels, stringsAsFactors = FALSE) %>%
  left_join(indTable_f3 %>% select(id, bars_plot, species), by = "id") %>%
  mutate(bars_plot     = case_when(id == "SPRAT" ~ "SPRAT", TRUE ~ bars_plot),
         species       = case_when(id == "SPRAT" ~ "SPRAT", TRUE ~ species),
         group_display = unname(GENOTYPE_GROUP_DISPLAY[bars_plot]))

unmapped <- unique(geno_meta$bars_plot[is.na(geno_meta$group_display)])
if (length(unmapped)) {
  stop("bars_plot values with no entry in GENOTYPE_GROUP_DISPLAY: ",
       paste(unmapped, collapse = ", "))
}

# run-length encode the COLLAPSED label, not bars_plot -- otherwise Vancouver
# and Japan come out as two adjacent runs both labelled "Pacific Pacific herring",
# with a separator line drawn between them
runs_geno <- group_runs(geno_levels, geno_meta$group_display)
runs_sp   <- group_runs(geno_levels, geno_meta$species) %>%
  mutate(label = unname(species_display[group]),
         mid_i = (start_i + end_i) / 2)
n_ind <- length(geno_levels)

# SNP order must be explicit. as.character(POS) would otherwise be sorted
# alphabetically -- which happens to equal genomic order only while every
# position has the same number of digits. The gene blocks below depend on this.
snp_levels <- as.character(sort(unique(sec16b_missense_Genos_freq_filter$POS)))
n_snp      <- length(snp_levels)

# Genotype colours: 0 = yellow, 1 = teal, 2 = dark blue, NA = white.
# These are viridis' own 3-level values with 0 and 2 swapped, so the palette
# matches the previous scale_fill_viridis(discrete = TRUE) figure. NA is made an
# explicit factor level -- that is what gives it a legend key at all.
# Kept identical to Supplementary Figure 25 so the same genotype is the same
# colour in both figures.
genotype_levels <- c("0", "1", "2", "NA")
genotype_colors <- c("0"  = "#FDE725",
                     "1"  = "#21918C",
                     "2"  = "#440154",
                     "NA" = "white")
genotype_labels <- c("Homozygote AA", "Heterozygote", "Homozygote BB", "Missing")

geno_dat <- sec16b_missense_Genos_freq_filter %>%
  mutate(Genotype_f = factor(ifelse(is.na(Genotype), "NA", as.character(Genotype)),
                             levels = genotype_levels),
         POS_f      = factor(as.character(POS), levels = snp_levels))

plot_genotypes_lab <- geno_dat %>%
  ggplot() +
  rasterise(geom_tile(aes(y = factor(Individuals, levels = geno_levels),
                          x = POS_f, fill = Genotype_f)), dpi = 300) +
  # black, not white: "Missing" tiles are white, so white separators between
  # two missing-data blocks are invisible
  geom_hline(yintercept = head(runs_geno$boundary, -1), colour = "black",
             linewidth = 0.3) +
  scale_fill_manual(values = genotype_colors, labels = genotype_labels,
                    drop = FALSE, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(breaks = runs_geno$mid_row,
                   labels = runs_geno$group,      # already display text
                   expand = c(0, 0)) +
  ylab(NULL) +
  theme_classic() +
  theme(axis.title.x        = element_blank(),
        axis.text.y         = element_text(size = 5, colour = "black", lineheight = 0.8),
        axis.text.x         = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y        = element_line(colour = "grey40", linewidth = 0.2),
        axis.ticks.length.y = unit(0.5, "mm"),
        axis.line.y         = element_blank(),
        legend.position     = "none",
        plot.margin         = margin(t = 0, r = 0, b = 1, l = 1, unit = "mm"))

# The four-key legend, extracted so it can be placed under the figure
genotype_legend <- cowplot::get_legend(
  plot_genotypes_lab +
    theme(legend.position = "bottom",
          legend.title    = element_blank(),
          legend.text     = element_text(size = 5),
          legend.key.size = unit(2.5, "mm"),
          legend.key      = element_rect(colour = "grey40", linewidth = 0.2)) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(colour = "grey40")))
)

# Species strip, drawn to the RIGHT of the heatmap (the left is taken by the
# population labels). Without it nothing on the panel identifies which
# individuals are Atlantic herring -- the population names alone don't say.
#
# y-domain must match the heatmap's: scale_y_discrete(expand = c(0,0)) puts
# n rows in a panel running 0.5 to n+0.5.
species_strip <- ggplot(runs_sp) +
  geom_segment(aes(x = 0.03, xend = 0.03, y = start_i - 0.4, yend = end_i + 0.4),
               linewidth = 0.8) +
  geom_text(aes(x = 0.18, y = mid_i, label = label),
            angle = 0, hjust = 0, vjust = 0.5, size = 1.8, lineheight = 0.85) +
  scale_y_continuous(limits = c(0.5, n_ind + 0.5), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "mm"))

# ----------------------------------------------------------------------------- #
# 3b. Gene header above the genotype heatmap
# ----------------------------------------------------------------------------- #
# Which consecutive SNP columns belong to which gene. EDIT THESE COUNTS if the
# check below fails -- it prints the real number of SNPs and their positions.
gene_blocks <- c("SINCLA"     = 6,
                 "ZGC:109982" = 1,
                 "SEC16B"     = 8,
                 "ZNF648"     = 2,
                 "GLULA"      = 1,
                 "novel gene" = 5)

if (sum(gene_blocks) != n_snp) {
  stop("gene_blocks sums to ", sum(gene_blocks), " but the heatmap has ", n_snp,
       " SNPs. Adjust the counts. Positions in order:\n",
       paste(snp_levels, collapse = ", "))
}

blk_end   <- cumsum(gene_blocks)
blk_start <- blk_end - gene_blocks + 1L
blk_df    <- data.frame(gene  = names(gene_blocks),
                        x1    = blk_start - 0.4,
                        x2    = blk_end   + 0.4,
                        x_mid = (blk_start + blk_end) / 2,
                        stringsAsFactors = FALSE) %>%
  # alternate between two rows so neighbouring genes never sit at the same
  # height -- horizontal text, staggered vertically
  mutate(row   = rep_len(c(0, 1), n()),
         lab_y = ifelse(row == 0, 0.28, 0.72))

# x-domain MUST match the heatmap's: scale_x_discrete(expand = c(0,0)) puts
# n columns in a panel running 0.5 to n+0.5. Mismatch here silently shifts
# every label relative to its block.
gene_header <- ggplot(blk_df) +
  geom_segment(aes(x = x1, xend = x2, y = 0, yend = 0), linewidth = 0.4) +
  geom_segment(aes(x = x_mid, xend = x_mid, y = 0.03, yend = lab_y - 0.05),
               linewidth = 0.15, colour = "grey50") +
  geom_text(aes(x = x_mid, y = lab_y, label = gene),
            angle = 0, hjust = 0.5, vjust = 0, size = 1.8) +
  scale_x_continuous(limits = c(0.5, n_snp + 0.5), expand = c(0, 0)) +
  ylim(-0.05, 1.1) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(t = 1, r = 2, b = 0, l = 1, unit = "mm"))

# Offsets in row/column units, proportional so they scale with the data size.
# These four numbers are the only tuning knobs for how far the strips sit.
HEADER_BASE <- n_ind * 1.01    # gene bracket line, just above the top tile
HEADER_LAB1 <- n_ind * 1.04    # lower row of gene labels
HEADER_LAB2 <- n_ind * 1.11    # upper row
STRIP_X     <- n_snp * 1.04    # species bracket, just right of the last column
STRIP_LAB   <- n_snp * 1.09    # species label text

blk_df <- blk_df %>%
  mutate(lab_y = ifelse(row == 0, HEADER_LAB1, HEADER_LAB2))

plot_genotypes_lab <- geno_dat %>%
  ggplot() +
  rasterise(geom_tile(aes(y = factor(Individuals, levels = geno_levels),
                          x = POS_f, fill = Genotype_f)), dpi = 300) +
  geom_hline(yintercept = head(runs_geno$boundary, -1), colour = "black",
             linewidth = 0.3) +
  
  # --- gene header, drawn above the tiles ---
  geom_segment(data = blk_df, inherit.aes = FALSE,
               aes(x = x1, xend = x2, y = HEADER_BASE, yend = HEADER_BASE),
               linewidth = 0.4) +
  geom_segment(data = blk_df, inherit.aes = FALSE,
               aes(x = x_mid, xend = x_mid, y = HEADER_BASE, yend = lab_y - 0.5),
               linewidth = 0.15, colour = "grey50") +
  geom_text(data = blk_df, inherit.aes = FALSE,
            aes(x = x_mid, y = lab_y, label = gene),
            size = 1.8, hjust = 0.5, vjust = 0) +
  
  # --- species strip, drawn to the right of the tiles ---
  geom_segment(data = runs_sp, inherit.aes = FALSE,
               aes(x = STRIP_X, xend = STRIP_X,
                   y = start_i - 0.4, yend = end_i + 0.4),
               linewidth = 0.8) +
  geom_text(data = runs_sp, inherit.aes = FALSE,
            aes(x = STRIP_LAB, y = mid_i, label = label),
            size = 1.8, hjust = 0, vjust = 0.5, lineheight = 0.85) +
  
  scale_fill_manual(values = genotype_colors, labels = genotype_labels,
                    drop = FALSE, name = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(breaks = runs_geno$mid_row,
                   labels = runs_geno$group,
                   expand = c(0, 0)) +
  # crop the panel to the tiles; clip = "off" lets the strips draw outside it
  coord_cartesian(xlim = c(0.5, n_snp + 0.5),
                  ylim = c(0.5, n_ind + 0.5),
                  clip = "off") +
  ylab(NULL) +
  theme_classic() +
  theme(axis.title.x        = element_blank(),
        axis.text.y         = element_text(size = 5, colour = "black", lineheight = 0.8),
        axis.text.x         = element_text(size = 4, angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.y        = element_line(colour = "grey40", linewidth = 0.2),
        axis.ticks.length.y = unit(0.5, "mm"),
        axis.line.y         = element_blank(),
        legend.position     = "none",
        # room for the strips: t for the gene header, r for the species labels
        plot.margin         = margin(t = 9, r = 16, b = 1, l = 1, unit = "mm"))
# ----------------------------------------------------------------------------- #
# 4. Gene names on the genomic track
# ----------------------------------------------------------------------------- #
# NOTE: in the original script the `yy` stagger vector is computed but never used
# -- geom_text() has a hardcoded y = 2 -- which is why the names collide.
#
# ggrepel measures the ACTUAL rendered text extents and iterates until nothing
# overlaps, which a width estimate can't do. Genes are filtered to the visible
# window first, otherwise repel spends its effort on labels that get cropped.
gene_label_size <- 2.0

genes_in_view <- scan1_match_df %>%
  filter(seqnames == chr) %>%
  filter(as.numeric(end) / 1e6   >= s / 1e6 - padding,
         as.numeric(start) / 1e6 <= e / 1e6 + padding) %>%
  mutate(start_mb = as.numeric(start) / 1e6,
         end_mb   = as.numeric(end) / 1e6)

message("gene track: ", nrow(genes_in_view), " genes in view")

gene_names_lab <- ggplot(genes_in_view) +
  geom_rect(data = candidate_matrix, inherit.aes = FALSE,
            mapping = aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 3),
            fill = "azure3", colour = NA, alpha = 0.5) +
  geom_rect(aes(xmin = start_mb, xmax = end_mb, ymin = 0.05, ymax = 0.35),
            colour = "black", fill = "lightblue", linewidth = 0.25) +
  geom_text_repel(aes(x = start_mb, y = 0.35, label = external_gene_name),
                  size               = gene_label_size,
                  nudge_y            = 0.6,
                  direction          = "both",   # stagger horizontally AND vertically
                  segment.size       = 0.15,
                  segment.colour     = "grey50",
                  min.segment.length = 0,        # always draw the leader line
                  box.padding        = 0.15,
                  point.padding      = 0,
                  max.overlaps       = Inf,
                  force              = 2,
                  seed               = 42) +     # reproducible across runs
  xlim(s / 1e6 - padding, e / 1e6 + padding) +
  ylim(0, 3) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.line  = element_blank(), axis.ticks = element_blank(),
        text = element_text(size = 6), legend.position = "none")

# ----------------------------------------------------------------------------- #
# 5. Re-assemble Figure 3
# ----------------------------------------------------------------------------- #
# Two blocks, because they have incompatible x axes: the top five panels share a
# genomic Mb axis; the genotype heatmap's x axis is discrete SNPs, so aligning
# it to them would be meaningless.
genomic_block <- plot_grid(
  gene_names_lab,
  plot_csq_BS_AS,
  plot_pops_freq_df_lab,
  plot_xpehh_2,
  twisst_plot,
  nrow        = 5,
  align       = "v",
  axis        = "lr",
  rel_heights = c(1.0, 1.0, 2.5, 1.0, 1.0)
)

# The genotype block is a 2x2 grid built in ONE plot_grid call:
#
#     gene_header        |  (blank)
#     plot_genotypes_lab |  species_strip
#
# align = "hv" with axis = "tblr" aligns every panel in both directions at once,
# so the gene brackets land over their SNP columns AND the species brackets land
# beside the right individuals. Doing this as nested calls does NOT work --
# cowplot cannot see through a composite to align against the panel inside it,
# which is what threw the gene header off before.
blank_corner <- ggplot() + theme_void()

genotype_block <- plot_genotypes_lab

# <- plot_grid(
#   gene_header,        blank_corner,
#   plot_genotypes_lab, species_strip,
#   ncol        = 2,
#   align       = "hv",
#   axis        = "tblr",
#   rel_widths  = c(1, 0.5),
#   rel_heights = c(4, 6)
# )

figure3 <- plot_grid(
  genomic_block,
  genotype_block,
  genotype_legend,
  ncol        = 1,
  rel_heights = c(6.5, 3, 0.4)
)

figure3

ggsave(genotype_block, filename = "~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure3/tmp_geno.pdf",
       width = 110, height = 65, units = "mm")

ggsave(figure3,
       filename = "~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure3/sec16b_panel_20260826_tmp.pdf",
       height = 188, width = 110, units = "mm")

ggsave(figure3,
       filename = "~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure3/sec16b_panel_20260826.png",
       height = 210, width = 110, units = "mm", dpi = 400)

