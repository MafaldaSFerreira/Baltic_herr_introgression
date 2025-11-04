# Plot divergence stats between populations
library(tidyverse)
library(biomaRt)
library(gridExtra)
library(GenomicRanges)

setwd("Manuscript/Figshare/")
# Whole genome plots Supplementary Figure 9 ####
# Overlaps with Introgression Regions.

## Input files ####
# dxy Maf 0.05
table_dxy<- read.table("7.pixy/results_WGS_dxy/2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_dxy.txt", header=T)

# Introgression regions 
intro_reg<-read.table(header=T, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

intro_reg_gr <- GRanges(seqnames=intro_reg$seqnames, IRanges(start=intro_reg$start, end=intro_reg$end))


## Process Dxy tables ####
## Dxy in regions of introgression: ####
table_dxy_Baltic_Spring_Canada_Spring <- table_dxy %>% filter(no_sites>=10000) %>% filter(pop1=="Baltic_Spring" & pop2=="Canada_Spring")
table_dxy_Baltic_Spring_WhiteSea <- table_dxy %>% filter(no_sites>=10000) %>% filter(pop1=="Baltic_Spring" & pop2=="WhiteSea")
table_dxy_Baltic_Spring_BarentSea <- table_dxy %>% filter(no_sites>=10000) %>% filter(pop1=="Baltic_Spring" & pop2=="BarentSea")

table_dxy_Baltic_Spring_Vancouver <- table_dxy %>% filter(no_sites>=10000) %>% filter(pop1=="Baltic_Spring" & pop2=="Vancouver")
table_dxy_Baltic_Spring_Vancouver$type<-"Vancouver"

table_dxy_Baltic_Spring_WhiteSea$mid <- (table_dxy_Baltic_Spring_WhiteSea$window_pos_1 + table_dxy_Baltic_Spring_WhiteSea$window_pos_2)/2
table_dxy_Baltic_Spring_Vancouver$mid <- (table_dxy_Baltic_Spring_Vancouver$window_pos_1 + table_dxy_Baltic_Spring_Vancouver$window_pos_2)/2
table_dxy_Baltic_Spring_BarentSea$mid <- (table_dxy_Baltic_Spring_BarentSea$window_pos_1 + table_dxy_Baltic_Spring_BarentSea$window_pos_2)/2


# convert into range:
table_dxy_Baltic_Spring_WhiteSea_gr<-GRanges(seqnames=table_dxy_Baltic_Spring_WhiteSea$chromosome, 
                                             IRanges(start=table_dxy_Baltic_Spring_WhiteSea$window_pos_1, end=table_dxy_Baltic_Spring_WhiteSea$window_pos_2),
                                             avg_dxy=table_dxy_Baltic_Spring_WhiteSea$avg_dxy,
                                             mid=table_dxy_Baltic_Spring_WhiteSea$mid)

table_dxy_Baltic_Spring_BarentSea_gr<-GRanges(seqnames=table_dxy_Baltic_Spring_BarentSea$chromosome, 
                                             IRanges(start=table_dxy_Baltic_Spring_BarentSea$window_pos_1, end=table_dxy_Baltic_Spring_BarentSea$window_pos_2),
                                             avg_dxy=table_dxy_Baltic_Spring_BarentSea$avg_dxy,
                                             mid=table_dxy_Baltic_Spring_BarentSea$mid)

# Find overlaps between both:
overlaps<-findOverlaps(intro_reg_gr, table_dxy_Baltic_Spring_WhiteSea_gr)

table_dxy_Baltic_Spring_WhiteSea$type<-"NotIntrogressed"

table_dxy_Baltic_Spring_WhiteSea[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

#table_dxy_Baltic_Spring_WhiteSea$type<-ifelse(table_dxy_Baltic_Spring_WhiteSea$type=="Introgressed", paste0(table_dxy_Baltic_Spring_WhiteSea$type,table_dxy_Baltic_Spring_WhiteSea$chromosome), table_dxy_Baltic_Spring_WhiteSea$type)

dxy_to_density<-rbind(table_dxy_Baltic_Spring_WhiteSea,table_dxy_Baltic_Spring_Vancouver)

dxy_boxplot<-ggplot()+
  geom_violin(data=dxy_to_density, aes(y=avg_dxy, x=type, fill=type), )+
  geom_boxplot(data=dxy_to_density, aes(y=avg_dxy, x=type, fill=type), width=0.08, outlier.shape = NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("divergence (dxy)")+
  scale_x_discrete(label=c("introgressed \nregions", "genome without \nintrogressed regions", "Atl herr Baltic Spring \nvs\nPac herr Vancouver"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.y = element_text(size=14, color="black"),
        legend.position = "none")

#ggsave(dxy_boxplot, file="figures/scan1_intro_reg_vs_WG_vs_Vancouver_boxplot_dxy_distb.png", dpi=300, width = 8, height = 8)
#ggsave(dxy_boxplot, file="figures/scan1_intro_reg_vs_WG_vs_Vancouver_boxplot_dxy_distb.pdf", width = 8, height = 8)

dxy_boxplot <- ggplot()+
  geom_violin(data=table_dxy_Baltic_Spring_WhiteSea, aes(y=avg_dxy, x=type, fill=type), )+
  geom_boxplot(data=table_dxy_Baltic_Spring_WhiteSea, aes(y=avg_dxy, x=type, fill=type), width=0.08, outlier.shape = NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("divergence (dxy)")+
  scale_x_discrete(label=c("Introgressed \nregions", "Genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position = "none")

# Using Barents Sea instead of White Sea

# Find overlaps between both:
overlaps<-findOverlaps(intro_reg_gr, table_dxy_Baltic_Spring_BarentSea_gr)

table_dxy_Baltic_Spring_BarentSea$type<-"NotIntrogressed"

table_dxy_Baltic_Spring_BarentSea[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

#table_dxy_Baltic_Spring_BarentSea$type<-ifelse(table_dxy_Baltic_Spring_BarentSea$type=="Introgressed", paste0(table_dxy_Baltic_Spring_BarentSea$type,table_dxy_Baltic_Spring_BarentSea$chromosome), table_dxy_Baltic_Spring_BarentSea$type)

dxy_to_density<-rbind(table_dxy_Baltic_Spring_BarentSea,table_dxy_Baltic_Spring_Vancouver)

ggplot()+
  geom_violin(data=dxy_to_density, aes(y=avg_dxy, x=type, fill=type), )+
  geom_boxplot(data=dxy_to_density, aes(y=avg_dxy, x=type, fill=type), width=0.08, outlier.shape = NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("divergence (dxy)")+
  scale_x_discrete(label=c("introgressed \nregions", "genome without \nintrogressed regions", "Atl herr Baltic Spring \nvs\nPac herr Vancouver"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=14, color="black"),
        axis.text.y = element_text(size=14, color="black"),
        axis.title.y = element_text(size=14, color="black"),
        legend.position = "none")

dxy_boxplot<- ggplot()+
  geom_violin(data=table_dxy_Baltic_Spring_BarentSea, aes(y=avg_dxy, x=type, fill=type), )+
  geom_boxplot(data=table_dxy_Baltic_Spring_BarentSea, aes(y=avg_dxy, x=type, fill=type), width=0.08, outlier.shape = NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("divergence (dxy)")+
  scale_x_discrete(label=c("Introgressed \nregions", "Genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position = "none")

ggplot()+
  geom_density(data=table_dxy_Baltic_Spring_BarentSea, aes(x=avg_dxy, fill=type), alpha=0.6)

# Ok, it makes no sense to not have these plots together... 
table_dxy_Baltic_Spring_BarentSea[table_dxy_Baltic_Spring_BarentSea$type == "Introgressed",]$type <- "IntrogressedBS"
table_dxy_Baltic_Spring_BarentSea[table_dxy_Baltic_Spring_BarentSea$type == "NotIntrogressed",]$type <- "NotIntrogressedBS"

WhiteSeaBarentsSea <- rbind(table_dxy_Baltic_Spring_BarentSea, table_dxy_Baltic_Spring_WhiteSea)
WhiteSeaBarentsSea$type <- factor(WhiteSeaBarentsSea$type, levels=c("Introgressed", "NotIntrogressed", "IntrogressedBS", "NotIntrogressedBS"))

WhiteSeaBarentsSea %>% group_by(type) %>% summarise(mean(avg_dxy)*100)

# FOR MAF5:
# # A tibble: 4 × 2
# type              `mean(avg_dxy) * 100`
# <fct>                             <dbl>
#   1 Introgressed                      0.320
# 2 NotIntrogressed                   0.396
# 3 IntrogressedBS                    0.319
# 4 NotIntrogressedBS                 0.408

dxy_boxplot<-ggplot()+
  geom_violin(data=WhiteSeaBarentsSea, aes(y=avg_dxy, x=type, fill=type), )+
  geom_boxplot(data=WhiteSeaBarentsSea, aes(y=avg_dxy, x=type, fill=type), width=0.08, outlier.shape = NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("divergence (dxy)")+
  scale_x_discrete(label=c("Introgressed \nregions", "Genome without \nintrogressed regions", "Introgressed \nregions", "Genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position = "none")

ggsave(dxy_boxplot, file="figures/scan1_intro_reg_vs_WG_boxplot_dxy_WhiteSea_&_BarentsSea_distb_2025-10-18_maf5.pdf", units="mm", width = 90, height = 50)

signif(t.test(avg_dxy ~ type, data=table_dxy_Baltic_Spring_WhiteSea, alternative="less")$p.value,3)
signif(t.test(avg_dxy ~ type, data=table_dxy_Baltic_Spring_BarentSea, alternative="less")$p.value,3)
