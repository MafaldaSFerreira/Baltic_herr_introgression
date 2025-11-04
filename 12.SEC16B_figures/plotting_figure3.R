# This R code allows to re-plot Figure 3 in the Baltic introgression paper
# All input files can be found either on GitHub or Figshare (if too big)

setwd("Manuscript/Figshare/")

library(ggrastr)
library(tidyverse)
library(cowplot)
library(viridis)

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
intro_reg<-read.table(header=T,"4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txtscan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

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


# PAPER PLOTS #####

# Define regoin to plot
padding <- 0.1
padding_mb <- 1e5

c<-"chr10"
s<-25040001
e<-25220000
chr<-str_remove(c, "chr")

# Subset the input files 
# Positive selection
tmp_xpehh<-xpehh_results_sig %>% filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

# ChiSQ
tmp_csq_BS_AS <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

# Allele Frequencies 
tmp_pops_freq_df <- pops_freq_df %>% 
  filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb) %>%
  pivot_longer(cols = 3:62, names_to = "populations", values_to = "frequencies")

# Define the introgression region:
candidate_matrix <- intro_reg %>% filter(seqnames == c & start >= s-padding_mb & end <= e+ padding_mb)

# Positive selection 
plot_xpehh_2<-ggplot(tmp_xpehh)+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1),
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
plot_csq_BS_AS <- ggplot(tmp_csq_BS_AS)+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_point(aes(x=POS/1e6, y=-log10(chisq_p), color=factor(-log10(chisq_p) > 9.1)), size=0.5)+
  scale_color_manual(values=c("gray", "black"))+
  #geom_line(aes(x=POS/1e6, y=rollmean(dAF_BS_vs_AS, 15, na.pad=TRUE))) +
  theme_classic()+
  xlim(s/1e6-padding, e/1e6 + padding)+
  xlab(paste0(chr," position (Mb)"))+
  ylab(expression(paste(-log[10], " (P value)") ) )+
  theme(legend.position = "none",
        axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_blank())

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
ggsave(all_plots, filename=paste0("sec16b_panel_2025-08-09.pdf"), height = 110, width = 100, units = "mm")

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
sec16B_snps <- moderate_mutations_50kb %>% filter(Chromosome=="chr10") %>% dplyr::select(Position)

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
  filter(CHROM=="chr10") %>% filter(POS %in% sec16B_snps$Position) %>%
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
    axis.text.y = element_text(size=2),
    axis.text.x = element_text(size=6, angle=90, hjust=0.5, vjust=0.5),
    axis.title.y = element_blank(),
    legend.position = "none")

ggsave(plot_genotypes, filename="sec16B_genotypes_big_2025-09-13.pdf", height = 100, width = 190, units = "mm")

