# This R code plots fdM statistics from Dsuite dinvestigate results
# Results in Supplementary Materials of the Baltic herring introgression paper
# Supplementary Figure 10

setwd("Manuscript/Figshare/")

# Libraries ####
library(tidyverse)
library(cowplot)


# Read introgression regions file
# 2025-06-24
# Use introgression regions scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb
intro_reg<-read.table(header=T, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

# fd values
table_BalticSpring<-read.table("3.dsuite/results_dinvestigate/Canada_Spring_Baltic_Spring_WhiteSea_localFstats_run_1004_2023-12-01_w50_s25_50_25.txt", header=T)
table_BalticSpring$mid <- (table_BalticSpring$windowEnd+table_BalticSpring$windowStart)/2

table_BalticSpring
head(table_BalticSpring)

# convert fd and fdM to zero if D is negative
for (x in 1:nrow(table_BalticSpring)){
  table_BalticSpring[x,5] = ifelse(table_BalticSpring[x,4] < 0, 0, table_BalticSpring[x,5])
}

for (x in 1:nrow(table_BalticSpring)){
  table_BalticSpring[x,6] = ifelse(table_BalticSpring[x,4] < 0, 0, table_BalticSpring[x,6])
}

# Zscore 
table_BalticSpring$Zscore<-scale(table_BalticSpring$f_dM,center=T,scale=T)
table_BalticSpring$Zscore_fd<-scale(table_BalticSpring$f_d,center=T,scale=T)
table_BalticSpring$Zscore_df<-scale(table_BalticSpring$d_f,center=T,scale=T)
table_BalticSpring$Zscore_fdM<-scale(table_BalticSpring$f_dM,center=T,scale=T)

# Convert tables to GRanges:
table_BalticSpring_gr<-GRanges(seqnames = table_BalticSpring$chr, 
                               IRanges(start = table_BalticSpring$windowStart, end = table_BalticSpring$windowEnd),
                               f_dM=table_BalticSpring$f_dM,
                               Zscore_fdM=table_BalticSpring$Zscore_fdM,
                               mid=table_BalticSpring$mid)

intro_reg_gr <- GRanges(seqnames=intro_reg$seqnames, IRanges(start=intro_reg$start, end=intro_reg$end))

# Find overlaps between both:
overlaps<-findOverlaps(intro_reg_gr, table_BalticSpring_gr)

table_BalticSpring$type <- "NotIntrogressed"

table_BalticSpring[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

plot_fdM<-ggplot(table_BalticSpring, aes(x=type, y=f_dM, fill=type))+
  geom_boxplot()+
  #geom_violin(fill=NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("fraction of admixture (f_dM)")+
  scale_x_discrete(label=c("introgressed \nregions", "genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position="none")

signif(t.test(f_dM ~ type, data=table_BalticSpring, alternative="greater")$p.value,3)

## Let's read the Barents sea file:

table_BarentSea<-read.table("3.dsuite/results_dinvestigate/Canada_Spring_Baltic_Spring_BarentSea_localFstats_run_1004_2023-12-01_w50_s25_50_25.txt", header=T)
table_BarentSea$mid <- (table_BarentSea$windowEnd+table_BarentSea$windowStart)/2

table_BarentSea
head(table_BarentSea)

# convert fd and fdM to zero if D is negative
for (x in 1:nrow(table_BarentSea)){
  table_BarentSea[x,5] = ifelse(table_BarentSea[x,4] < 0, 0, table_BarentSea[x,5])
}

for (x in 1:nrow(table_BarentSea)){
  table_BarentSea[x,6] = ifelse(table_BarentSea[x,4] < 0, 0, table_BarentSea[x,6])
}

# Zscore
table_BarentSea$Zscore<-scale(table_BarentSea$f_dM,center=T,scale=T)
table_BarentSea$Zscore_fd<-scale(table_BarentSea$f_d,center=T,scale=T)
table_BarentSea$Zscore_df<-scale(table_BarentSea$d_f,center=T,scale=T)
table_BarentSea$Zscore_fdM<-scale(table_BarentSea$f_dM,center=T,scale=T)

# Convert tables to GRanges:
table_BarentSea_gr<-GRanges(seqnames = table_BarentSea$chr, 
                               IRanges(start = table_BarentSea$windowStart, end = table_BarentSea$windowEnd),
                               f_dM=table_BarentSea$f_dM,
                               Zscore_fdM=table_BarentSea$Zscore_fdM,
                               mid=table_BarentSea$mid)

# Find overlaps between both:
overlaps<-findOverlaps(intro_reg_gr, table_BarentSea_gr)

table_BarentSea$type <- "NotIntrogressed"

table_BarentSea[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

ggplot(table_BarentSea, aes(x=type, y=f_dM, fill=type))+
  geom_boxplot()+
  #geom_violin(fill=NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("fraction of admixture (f_dM)")+
  scale_x_discrete(label=c("introgressed \nregions", "genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position="none")

# Let's join both tables
table_BarentSea[table_BarentSea$type == "Introgressed",]$type <- "IntrogressedBS"
table_BarentSea[table_BarentSea$type == "NotIntrogressed",]$type <- "NotIntrogressedBS"

table_BarentSea_WhiteSea <- rbind(table_BarentSea, table_BalticSpring)
table_BarentSea_WhiteSea$type <- factor(table_BarentSea_WhiteSea$type, levels=c("Introgressed", "NotIntrogressed", "IntrogressedBS", "NotIntrogressedBS"))


total_boxplot <- ggplot(table_BarentSea_WhiteSea, aes(x=type, y=f_dM, fill=type))+
  geom_boxplot()+
  #geom_violin(fill=NA)+
  scale_fill_brewer(palette="Dark2")+
  ylab("fraction of admixture (f_dM)")+
  scale_x_discrete(label=c("Introgressed \nregions", "Genome without \nintrogressed regions", "Introgressed \nregions", "Genome without \nintrogressed regions"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position="none")

# PAPER FIGURE 2025-06-25
ggsave(total_boxplot, file="figures/scan1_intro_reg_vs_WG_boxplot_fdM_WhiteSea_&_BarentsSea_distb_2025-06-27.pdf", units="mm", width = 90, height = 50)

signif(t.test(f_dM ~ type, data=table_BalticSpring, alternative="greater")$p.value,3)
#[1] 1.85e-66
signif(t.test(f_dM ~ type, data=table_BarentSea, alternative="greater")$p.value,3)
#[1] 6.29e-57
