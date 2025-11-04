# This R code allows to re-plot Figure 3 in the Baltic introgression paper
# All input files can be found either on GitHub or Figshare (if too big)

setwd("~/Documents/Postdoc/Project_Herring/Introgression/Manuscript/Figshare/")

# INPUT FILES ####

library(ggrastr)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggsci)
library(gridExtra)

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

annotation_ensembl114 <- read.table("13.THRB_figures/annotations/Ch_v2.0.2v2_ensembl114_chr19-5.88Mb-6.88Mb_annotation_genestart.txt", header=T, sep="\t")

annotation_ensembl114_min_max <- annotation_ensembl114 %>% 
  group_by(Transcript.stable.ID) %>%
  summarise(min=min(Exon.region.start..bp.),
            max=max(Exon.region.end..bp.))

annotation_ensembl114_min_max$n <- 1:nrow(annotation_ensembl114_min_max)

annotation_ensembl114_min_max_to_plot <- left_join(annotation_ensembl114, annotation_ensembl114_min_max, by="Transcript.stable.ID")


annotation_ensembl109 <- read.table("13.THRB_figures/annotations/Ch_v2.0.2_ensembl109_chr19-5.88Mb-6.88Mb_annotation_genestart.txt", header=T, sep="\t")

annotation_ensembl109_min_max <- annotation_ensembl109 %>% 
  group_by(Transcript.stable.ID) %>%
  summarise(min=min(Exon.region.start..bp.),
            max=max(Exon.region.end..bp.))

annotation_ensembl109_min_max$n <- 1:nrow(annotation_ensembl109_min_max)

annotation_ensembl109_min_max_to_plot <- left_join(annotation_ensembl109, annotation_ensembl109_min_max, by="Transcript.stable.ID")

# PAPER PLOTS #####

# Define regoin to plot
# padding <- 0.1
# padding_mb <- 1e5

padding <- 0.25
padding_mb <- 2.5e5

c<-"chr19"
s<-6360001
e<-6380000
chr<-str_remove(c, "chr")

s <- 6.25e6
e <- 6.5e6
chr<-str_remove(c, "chr")

# Subset the input files 
# Positive selection
tmp_xpehh<-xpehh_results_sig %>% filter(CHROM==c & POS >= s & POS <= e)
  
# Chi-square
tmp_csq_BS_AS <- baltic_spring_vs_atlantic_spring_chiseq %>% 
     filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

# Allele Frequency heat map
tmp_pops_freq_df <- pops_freq_df %>% 
    #filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb) %>%
    filter(CHROM==c & POS >= s & POS <= e ) %>%
    pivot_longer(cols = 3:62, names_to = "populations", values_to = "frequencies")
  
# Define the introgression region:
candidate_matrix<-data.frame(matrix(c(6360001/1e6,6380000/1e6),ncol=2,nrow=1))
colnames(candidate_matrix) <- c("start", "end")

# Positive selection 
plot_xpehh_2<-ggplot(tmp_xpehh)+
    geom_rect(data=candidate_matrix,mapping=aes(xmin=start,xmax=end,ymin=-6,ymax=6),
              fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
    geom_point(aes(x=POS/1e6, y=std_XPEHH, color=factor(abs(std_XPEHH) >= 2)), size=0.5)+
    scale_color_manual(values=c("black", "red"))+
    ylim(-6,6)+
    xlim(s/1e6, e/1e6)+
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
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=45),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_point(aes(x=POS/1e6, y=-log10(chisq_p), color=factor(-log10(chisq_p) > 9.1)), size=0.5)+
  scale_color_manual(values=c("gray", "black"))+
  #geom_line(aes(x=POS/1e6, y=rollmean(dAF_BS_vs_AS, 15, na.pad=TRUE))) +
  theme_classic()+
  xlim(s/1e6-padding, e/1e6 + padding)+
  #xlim(s, e)+
  xlab(paste0(chr," position (Mb)"))+
  ylab(expression(paste(-log[10], " (P value)") ) )+
  theme(legend.position = "none",
        axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_blank())
  
# Allele frequency heat map
plot_pops_freq_df <- tmp_pops_freq_df %>% 
  ggplot()+
  rasterise(geom_tile(aes(y=factor(populations, level=rev(pool_order$V1)), x=POS/1e6, color=frequencies)), dpi=300) + 
  scale_color_viridis(direction = -1) +
  ylab("Populations")+
  xlab("Allele frequencies")+
  xlim(s/1e6, e/1e6) + 
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

genes<-scan1_match_df %>%
  filter(seqnames==chr) %>%
  ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=3),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_rect(aes(xmin=as.numeric(start)/1e6,xmax=as.numeric(end)/1e6
                ,ymin=1,ymax=1.25),colour="black",fill="lightblue",linewidth=0.25)+
  geom_text(aes(label=external_gene_name, x=start/1e6, y=2, angle=-45), size=2.5)+
  ylim(0,3)+
  xlim(s/1e6-padding, e/1e6+padding)+
  xlab(paste0("Chr ",chr," position (Mb)"))+
  ylab("gene models")+
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=6),
        legend.position = "none")

# Gene Models With Ensembl 114
#annotation_ensembl114_min_max_to_plot_thrb <- annotation_ensembl114_min_max_to_plot %>% filter(Gene.name=="THRB")
  # genes <- ggplot() +
#   geom_rect(data=annotation_ensembl114_min_max_to_plot_thrb,
#             aes(xmin=min/1e6,
#                 xmax=max/1e6,ymin=n-0.1,ymax=n+0.1, fill=Transcript.stable.ID))+
#   geom_rect(data=annotation_ensembl114_min_max_to_plot_thrb,
#             aes(xmin=as.numeric(Exon.region.start..bp.)/1e6,
#                 xmax=as.numeric(Exon.region.end..bp.)/1e6,
#                 ymin=n-0.1,ymax=n+0.1, fill=Transcript.stable.ID), color="black") +
#   xlim(s/1e6-padding, e/1e6 + padding)+
#   theme_classic()+
#   theme(legend.position = "none")

# Twisst
twisst_plot <- df %>%
  filter(scaffold==c) %>%
  ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=1),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_line(aes(x=mid/1000000,y=topo1_frac),col="#882E72",size=0.5)+
  geom_line(aes(x=mid/1000000,y=topo2_frac),col="#E8601C",size=0.5)+
  geom_line(aes(x=mid/1000000,y=topo3_frac),col="#5289C7",size=0.5)+
  #xlim(s/1e6-padding, e/1e6+padding)+
  xlim(s/1e6, e/1e6)+
  ylab("topology \nsupport (%)") +
  xlab(paste0("Chr ", chr," position (Mb)"))+
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

# Plot all and save
all_plots <- plot_grid(genes, plot_csq_BS_AS, plot_pops_freq_df, plot_xpehh_2, twisst_plot, nrow=5, align = "v")
ggsave(all_plots, filename=paste0("~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure4/thrb_panel_2025-08-09_zoomin_AtlanticSpring.pdf"), height = 110, width = 100, units = "mm")
  
# Plot heatmap independently and bigger
plot_pops_freq_df <- tmp_pops_freq_df %>% 
  ggplot()+
  rasterise(geom_tile(aes(y=factor(populations, level=rev(pool_order$V1)), x=POS/1e6, color=frequencies)), dpi=300) + 
  scale_color_viridis(direction = -1) +
  ylab("Populations")+
  xlab("Allele frequencies")+
  xlim(s/1e6, e/1e6)+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6),
    axis.text.y = element_blank(),
    axis.title.y = element_text(size=6),
    legend.position = "none")
  
ggsave(plot_pops_freq_df, filename=paste0("~/Documents/Postdoc/Project_Herring/Introgression/join_figures/figure4/thrb_AF_heatmap_2025-08-09.pdf"), height = 50, width = 100, units = "mm")

# SUPPLEMENTARY FIGURES ####
## Supplementary Figure 17b Allele frequency heat map ####
thrb_snps <- moderate_mutations_50kb %>% filter(seqnames=="chr19") %>% dplyr::select(start)

thrb_missense_snps_freq <- tmp_pops_freq_df %>% filter(CHROM=="chr19") %>% filter(POS %in% thrb_snps$start)

thrb_missense_snps_freq<-merge(thrb_missense_snps_freq, pool_order, by.x="populations", by.y="V1")

thrb_missense_snps_freq_plot <- thrb_missense_snps_freq %>%
  ggplot()+
  geom_tile(aes(y=factor(populations, level=rev(pool_order$V1)), x=as.character(POS), 
                fill=frequencies)) + 
  scale_fill_viridis(direction = -1) +
  ylab("Populations")+
  xlab("Allele frequencies")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6, angle=90),
    axis.title.y = element_text(size=6),
    axis.text.y = element_text(size=6))

legend <- thrb_missense_snps_freq %>%
  ggplot()+
  geom_tile(aes(y=factor(populations, level=rev(pool_order$V1)), 
                          x=as.character(POS), fill=V2)) + 
  scale_fill_brewer(palette="Paired") +
  ylab("Populations")+
  xlab("Allele frequencies")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6, angle=90),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none")

all_plots <- cowplot::plot_grid(thrb_missense_snps_freq_plot, legend, 
                                ncol=2,
                                rel_widths = c(3, 1))

ggsave(thrb_missense_snps_freq_plot, filename="figure4/thrb_allele_frequencies_big_2025-08-06.pdf", height = 150, width = 70, units = "mm")
ggsave(plot = all_plots, filename="figure4/thrb_allele_frequencies_big_2025-08-13.pdf", 
       height = 150, width = 100, 
       units = "mm")

## Supplementary Figure 17a ####
c<-"chr19"
scan1_match_df %>% filter(external_gene_name=="thrb")
s <-  6344165
e <- 6477279
chr<-str_remove(c, "chr")
padding <- 0
padding_mb <- 0

candidate_matrix<-data.frame(matrix(c(6360001/1e6,6380000/1e6),ncol=2,nrow=1))
colnames(candidate_matrix) <- c("start", "end")

# Annotation:
annotation_ensembl114_min_max_to_plot_thrb <- annotation_ensembl114_min_max_to_plot %>% filter(Gene.name=="THRB")
annotation_ensembl109_min_max_to_plot_thrb <- annotation_ensembl109_min_max_to_plot %>% filter(Gene.name=="thrb")

# CHISQ
tmp_csq_BS_AS <- baltic_spring_vs_atlantic_spring_chiseq %>% 
  filter(CHROM==c & POS >= s-padding_mb & POS <= e + padding_mb)

tmp_csq_BS_AS_missense <- tmp_csq_BS_AS %>% filter(POS==6359267)

#Plots
genes114 <- ggplot() +
  geom_vline(xintercept=6359267/1e6,color="pink")+
  geom_rect(data=annotation_ensembl114_min_max_to_plot_thrb,
            aes(xmin=min/1e6,
                xmax=max/1e6,ymin=n-0.1,ymax=n+0.1),fill="lightgreen")+
  geom_rect(data=annotation_ensembl114_min_max_to_plot_thrb,
            aes(xmin=as.numeric(Exon.region.start..bp.)/1e6,
                xmax=as.numeric(Exon.region.end..bp.)/1e6,
                ymin=n-0.1,ymax=n+0.1), color="black",fill="lightgreen") +
  xlim(s/1e6-padding, e/1e6 + padding)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

genes109 <- ggplot() +
  geom_vline(xintercept=6359267/1e6,color="pink")+
  geom_rect(data=annotation_ensembl109_min_max_to_plot_thrb,
            aes(xmin=min/1e6,
                xmax=max/1e6,ymin=n-0.1,ymax=n+0.1), fill="lightblue")+
  geom_rect(data=annotation_ensembl109_min_max_to_plot_thrb,
            aes(xmin=as.numeric(Exon.region.start..bp.)/1e6,
                xmax=as.numeric(Exon.region.end..bp.)/1e6,
                ymin=n-0.1,ymax=n+0.1),fill="lightblue", color="black") +
  
  xlim(s/1e6-padding, e/1e6 + padding)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())


plot_csq_BS_AS <- ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start,xmax=end,ymin=0,ymax=45),
            fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
  geom_vline(xintercept=6359267/1e6,color="pink")+
  geom_point(data=tmp_csq_BS_AS, aes(x=POS/1e6, y=-log10(chisq_p), color=factor(-log10(chisq_p) > 9.1)), size=0.5)+
  geom_point(data=tmp_csq_BS_AS_missense, aes(x=POS/1e6, y=-log10(chisq_p), color="red"), size=0.8)+
  scale_color_manual(values=c("gray", "black", "red"))+
  
  #geom_line(aes(x=POS/1e6, y=rollmean(dAF_BS_vs_AS, 15, na.pad=TRUE))) +
  theme_classic()+
  xlim(s/1e6-padding, e/1e6 + padding)+
  #xlim(s, e)+
  xlab(paste0(chr," position (Mb)"))+
  ylab(expression(paste(-log[10], " (P value)") ) )+
  theme(legend.position = "none",
        axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_blank())

all_plot<-cowplot::plot_grid(genes114, genes109, plot_csq_BS_AS, nrow=3, align = "v")

ggsave(all_plot, filename="figure4/gene_models_chisq_BS_AS.pdf", 
       units="mm",
       width = 90, height = 90)

## Supplementary Figure 17c ####
genotypes <- read.table("1.mapping_variant_calling/genotypes/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.intro_reg.genotypes.2025-08-06.txt")
indv_v02<-read.table("1.mapping_variant_calling/genotypes/individuals_v01_20231116.txt")
header<-read.table("1.mapping_variant_calling/genotypes/header_vcf_20231116.txt")

# Pacific and Atlantic herring genotypes
colnames(genotypes)<-header[,1]

genotypes[genotypes=="0/0"]<-0
genotypes[genotypes=="0/1"]<-1
genotypes[genotypes=="1/0"]<-1
genotypes[genotypes=="1/1"]<-2
genotypes[genotypes=="./."]<-NA

# European sprat genotype:
sprat.genotype$V4 <- NULL
sprat.genotype[sprat.genotype=="0/0"]<-0
sprat.genotype[sprat.genotype=="0/1"]<-1
sprat.genotype[sprat.genotype=="1/0"]<-1
sprat.genotype[sprat.genotype=="1/1"]<-2
sprat.genotype[sprat.genotype=="./."]<-NA
colnames(sprat.genotype) <- c("CHROM", "POS", "SPRAT")

HerrSprat_genotypes <- left_join(genotypes, sprat.genotype, by=c("CHROM", "POS"))

thrb_missense_Genos_freq <- HerrSprat_genotypes %>% filter(CHROM=="chr19") %>% filter(POS %in% thrb_snps$start) %>%
  pivot_longer(cols=3:128, names_to="Individuals", values_to="Genotype")

thrb_missense_Genos_freq <- merge(thrb_missense_Genos_freq, indv_v02, by.x="Individuals", by.y="V1")

# Remove Balsfjord
thrb_missense_Genos_freq_filter <- thrb_missense_Genos_freq %>% filter(!Individuals %in% c("HWS61_Balsfjord_Atlantic", "HWS62_Balsfjord_Atlantic",
                                                                                               "HWS63_Balsfjord_Atlantic", "HWS64_Balsfjord_Atlantic"))

indv_v03 <- indv_v02 %>% filter(!V1 %in% c("HWS61_Balsfjord_Atlantic", "HWS62_Balsfjord_Atlantic",
                                           "HWS63_Balsfjord_Atlantic", "HWS64_Balsfjord_Atlantic"))


plot_genotypes <- thrb_missense_Genos_freq %>%
  ggplot()+
  rasterise(geom_tile(aes(y=factor(Individuals,level=indv_v02[,1]), x=as.character(POS), fill=Genotype)), dpi=300) + 
  #rasterise(geom_tile(aes(y=Individuals, x=as.character(POS), fill=Genotype)), dpi=300) + 
  #scale_fill_manual(values=c("#1F78B4", "#CCCCCC", "#D95F02", "black"))+
  scale_fill_viridis(discrete = TRUE) +
  ylab("Populations")+
  xlab("Genotypes")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=5, angle=90),
    axis.title.y = element_blank(),
    legend.position = "none")


legend2<- thrb_missense_Genos_freq %>%
  ggplot()+
  geom_tile(aes(y=factor(Individuals,level=indv_v02[,1]), 
                x=as.character(POS), fill=V3)) + 
  scale_fill_simpsons()+
  ylab("Populations")+
  xlab("Genotypes")+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=6, angle=90),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none")

all_plots3 <- cowplot::plot_grid(plot_genotypes, legend2, ncol=2, rel_widths = c(3,1))

ggsave(all_plots3, filename="figure4/thrb_genotypes_big_2025-08-13.pdf", height = 100, width = 70, units = "mm")


# FIGURE 4 Correlation RHO vs THRB ####
# I think the ref and alt in this table is swapped (when comparing with our actual references fasta)
# I want to plot the Baltic versus Baltic allele in our plots. 
# This means I need to plot the REF allele for THRB versus the ALT allele for RHO
# However, Jake's table is swapped, so I actually need to plot ALT for THRB against REF for RHO in his table
table_Jake_RHO_THRB <- read.table("13.THRB_figures/goodall_et_al_2025_data//jake_goodall-THRB/merged_refalt_freq_long_with_metadata.tsv",
                                  header=T, sep="\t")
head(table_Jake_RHO_THRB)

# REF=A and ALT=G
table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("4_11217660")) %>%
  filter(sample_n >=10) %>%
  ggplot()+
  geom_point(aes(y=ALT_FREQ, x=1, color=SEASON))

# REF=T and ALT=A
table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("4_11217516")) %>%
  filter(sample_n >=10) %>%
  ggplot()+
  geom_point(aes(y=ALT_FREQ, x=1, color=SEASON))

# For this SNP, REF=G and ALT=T
table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("19_6359267")) %>%
  filter(sample_n >=10) %>%
  ggplot()+
  geom_point(aes(y=ALT_FREQ, x=1, color=SEASON))

# Plot
p1 <- table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("19_6359267", "4_11217516")) %>%
  pivot_wider(names_from = "SNP", values_from = c("REF_FREQ", "ALT_FREQ", "CHR", "REF", "ALT")) %>%
  filter(sample_n >=10) %>%
  ggplot()+
  geom_point(aes(y=REF_FREQ_4_11217516,x=ALT_FREQ_19_6359267, color=SEASON), size=0.5)+
  geom_smooth(aes(y=REF_FREQ_4_11217516,x=ALT_FREQ_19_6359267), method = "lm", se = TRUE, color = "black", size = 0.6) +  # Single trendline in black
  xlab("THRB Q40H frequency")+
  ylab("RHO F261Y frequency") + 
  theme_classic()+
  theme(legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size=6))


p2 <- table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("19_6359267", "4_11217660")) %>%
  pivot_wider(names_from = "SNP", values_from = c("REF_FREQ", "ALT_FREQ", "CHR", "REF", "ALT")) %>%
  filter(sample_n >=10) %>%
  ggplot()+
  geom_point(aes(y=REF_FREQ_4_11217660,x=ALT_FREQ_19_6359267, color=SEASON), size=0.5)+
  geom_smooth(aes(y=REF_FREQ_4_11217660,x=ALT_FREQ_19_6359267), method = "lm", se = TRUE, color = "black", size = 0.6) +  # Single trendline in black
  xlab("THRB Q40H frequency")+
  ylab("RHO T213I frequency") + 
  theme_classic()+
  theme(legend.position = "none",
        axis.text=element_text(size=6),
        axis.title = element_text(size=6))


all_plots <- grid.arrange(p1, p2, ncol=2)
ggsave(all_plots, filename="figure4/THRB_RHO_all_populations_corr_Jake_data.pdf", units="mm", width=80, height = 40)

SNP4_11217660_vs_19_6359267 <- table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("19_6359267", "4_11217660")) %>%
  pivot_wider(names_from = "SNP", values_from = c("REF_FREQ", "ALT_FREQ", "CHR", "REF", "ALT")) %>%
  filter(sample_n >=10)

SNP4_11217516_vs_19_6359267 <- table_Jake_RHO_THRB %>% 
  filter(SNP %in% c("19_6359267", "4_11217516")) %>%
  pivot_wider(names_from = "SNP", values_from = c("REF_FREQ", "ALT_FREQ", "CHR", "REF", "ALT")) %>%
  filter(sample_n >=10)

lm1 <- summary(lm(ALT_FREQ_19_6359267 ~ REF_FREQ_4_11217660, SNP4_11217660_vs_19_6359267))
lm2 <- summary(lm(ALT_FREQ_19_6359267 ~ REF_FREQ_4_11217516, SNP4_11217516_vs_19_6359267))

lm1$adj.r.squared
lm2$adj.r.squared

lm1
# Call:
#   lm(formula = ALT_FREQ_19_6359267 ~ REF_FREQ_4_11217660, data = SNP4_11217660_vs_19_6359267)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.222702 -0.052117  0.008265  0.052976  0.170954 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.47026    0.02219   21.19   <2e-16 ***
#   REF_FREQ_4_11217660  0.52928    0.04392   12.05   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0806 on 147 degrees of freedom
# Multiple R-squared:  0.497,	Adjusted R-squared:  0.4936 
# F-statistic: 145.3 on 1 and 147 DF,  p-value: < 2.2e-16

lm2
# Call:
#   lm(formula = ALT_FREQ_19_6359267 ~ REF_FREQ_4_11217516, data = SNP4_11217516_vs_19_6359267)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.235117 -0.037362  0.005643  0.055178  0.153267 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.32364    0.02801   11.56   <2e-16 ***
#   REF_FREQ_4_11217516  0.52256    0.03558   14.69   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07236 on 147 degrees of freedom
# Multiple R-squared:  0.5947,	Adjusted R-squared:  0.5919 
# F-statistic: 215.7 on 1 and 147 DF,  p-value: < 2.2e-16
# PLOT Allele frequencies for the missense mutation:

