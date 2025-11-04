# Introgression Scan ####
# This is code to run the introgression scan in Ferreira et al Baltic herring introgression paper
# We ran several introgression scans as specified in Supplementary Table 2
# This code requires several input files and code which are available in this repository or on Figshare


setwd("Manuscript/Figshare/")

#save.image("4.introgression_scan/intermediate_files/introgression_scan_2023-11-20.RData")
load("4.introgression_scan/intermediate_files/introgression_scan_2023-11-20.RData")

# Necessary packages ####
# install.packages("diptest")
# install.packages("stringdist")
# install.packages("Baltic_herr_introgression/04.introgression_scan/introgression_scan/HaploDistScan_0.2.tar.gz", dependencies = T)
#BiocManager::install("GenomicRanges")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("biomaRt")

# Load libraries ####
require(stringdist)
require(HaploDistScan)
require(GenomicRanges)
require(GenomicFeatures)
require(tidyverse)
require(biomaRt)

# Load the Chromosome sizes ####
Ch_v2.0.2_sizes <- read.table("4.introgression_scan/inputs/Ch_v2.0.2.sizes", sep = "\t", stringsAsFactors = F, header = F)
names(Ch_v2.0.2_sizes) <- c("name", "size")
Ch_v2.0.2_sizes[,"offset"] <- c(0, cumsum(Ch_v2.0.2_sizes[,"size"])[-length(Ch_v2.0.2_sizes[,"size"])])
Ch_v2.0.2_sizes[,"col"] <- c("grey30", "grey70")[(1:length(Ch_v2.0.2_sizes$name) %% 2) + 1]

# Load data:
herring_125 <- generate_geno_df("4.introgression_scan/inputs/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5_GT.txt")
#save(herring_125, file="4.introgression_scan/intermediate_files//herring_125.RData")
load(file="4.introgression_scan/intermediate_files/herring_125.RData")

# Load sample list:
# I modified this file so that we have the correct names based on sample_metadata_20230714
sample_list<-read.table("4.introgression_scan/inputs/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5_sample.txt")
haplotypes_list<-paste0(rep(sample_list$V1,each=2),"_",c(1,2))
herring_125$sample_list<-haplotypes_list

# windows
h125_20k_w <- construct_window_annotation(herring_125$geno, interval = 2e4)

# ---------------------------------------- #
# SCAN1_V01_BALTIC_SPRING: Parent 1: NSSH_Norway and Canada Atlantic Spring; Parent 2: White sea + Pechora; Target: Baltic Spring ####
# ---------------------------------------- #

# This includes the F individuals that are also Baltic spring:

# Confirm we are selecting the correct individuals:
target_group ="B[MF][0-9]|F[0-9]_HastKar"; 
ref_1_group = "NSSH|Canada_Atlantic_Spring"; 
ref_2_group = "HWS[2345]"

length(c(herring_125$sample_list[grep(target_group, herring_125$sample_list)],
         herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)],
         herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]))

scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target= target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
#save(scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df.Rdata")
#load("4.introgression_scan/intermediate_files/scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan1_v01_baltic_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df))*2

# Using min_diff 20 and snp_cutoff 50
scan1_v01_baltic_v_subarctic_alt_ref_intro_filter2 <- introgression_plot_v3(scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df, 
                                                                            sample_list=herring_125$sample_list, 
                                                                            snp_numbers=scan1_v01_baltic_v_subarctic_alt_ref_diff_SNP_count_vec, 
                                                                            pdf_file="scan1_v01_baltic_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", 
                                                                            min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  
                                                                            target_re ="B[MF][0-9]|F[0-9]_HastKar")

scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2 <- complile_intro_lists(intr_obj = scan1_v01_baltic_v_subarctic_alt_ref_intro_filter2, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan1_v01_baltic_alt_ref_summary_filter2 <- summarize_intro(intro_lists = scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2, plot_dir = "outputs/scan1_v01_baltic_summary_regions_filter2/")
save(scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2, file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2.Rdata")
save(scan1_v01_baltic_alt_ref_summary_filter2, file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_alt_ref_summary_filter2.Rdata")

#load(file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2.Rdata")
#load(file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_alt_ref_summary_filter2.Rdata")

## REDUCE RANGES ####
scan1_v01_baltic_alt_ref_summary_filter2$rec_GR # this 99 ranges
GenomicRanges::reduce(scan1_v01_baltic_alt_ref_summary_filter2$rec_GR) #this has 42 ranges
GenomicRanges::reduce(scan1_v01_baltic_alt_ref_summary_filter2$rec_GR, min.gapwidth=100000) # 31 ranges

# If we make coverage a bit higher. So like, 7 30% of total haplotypes:
scan1_v01_baltic_alt_ref_summary_filter2_cov7<-scan1_v01_baltic_alt_ref_summary_filter2$cov[scan1_v01_baltic_alt_ref_summary_filter2$cov[, "cov"] > 7,]
scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr<-makeGRangesFromDataFrame(df=scan1_v01_baltic_alt_ref_summary_filter2_cov7, 
                                                                           seqnames.field = "group_name",
                                                                           keep.extra.columns = T)

sum(data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7)$width)
# 2080000 bp # USE THIS LENGTH IN PUBLICATION 2025-06-19
# Reduce with different min.gapwidth:
scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min10kb<-GenomicRanges::reduce(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr, min.gapwidth=10000)
scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb<-GenomicRanges::reduce(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr, min.gapwidth=50000)

sum(data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min10kb)$width)
# 2080000 bp

sum(data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)$width)
# 2120000 bp 

# 2025-10-24: Used for publication:
2120000/725670187*100 = 0.29 %

# Based on scan1_v01_baltic_alt_ref_summary_filter2_cov7 
2080000/725670187*100 = 0.28 

write.table(data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb), 
            file = "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt",
            col.names = T, row.names = F, sep = "\t", quote=F)

write.table(data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min10kb), 
            file = "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min10kb.txt",
            col.names = T, row.names = F, sep = "\t", quote=F)

# Collapse regions (before this code was in "join_plotting.R":
intro_reg_collapsed<-data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb[c(1,2,3,10,16,18,23,24,25),])
rownames(intro_reg_collapsed)<-1:nrow(intro_reg_collapsed)
intro_reg_collapsed[3,3]<-22060000
intro_reg_collapsed[4,3]<-27560000
intro_reg_collapsed[5,3]<-16240000
intro_reg_collapsed[6,3]<-15360000

write.table(intro_reg_collapsed, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.9regions.collapsed.txt")

# Convert into a bed file and add a buffer (before this code was in "recomb_heatmaps.R")
intro_reg_collapsed_bed<-data.frame(chrom=intro_reg_collapsed$seqnames, start=intro_reg_collapsed$start-1e6-1, end=intro_reg_collapsed$end+1e6)
write.table(intro_reg_collapsed_bed, quote = F, col.names = T, row.names = F, sep = "\t",
            file="4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed")



# ---------------------------------------- #
# SCAN2_V01_BALTIC_AUTUMN: Parent 1: NSSH_Norway and Canada Atlantic Spring; Parent 2: Pechora +  White sea; Target: Baltic Autumn ####
# ---------------------------------------- #

# Confirm we are selecting the correct individuals:
target_group="Fehmarn3_Fehmarn_Baltic_Autumn|Fehmarn44_Fehmarn_Baltic_Autumn|Fehmarn6_Fehmarn_Baltic_Autumn|Gavle100_Gavle_Baltic_Autumn|Gavle54_Gavle_Baltic_Autumn|Gavle98_Gavle_Baltic_Autumn"; 
ref_1_group = "NSSH|Canada_Atlantic_Spring"; 
ref_2_group = "HWS[2345]"

herring_125$sample_list[grep(target_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]

scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target=target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
#save(scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df.Rdata")
#load("4.introgression_scan/intermediate_files/scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan2_v01_baltic_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df))*2

scan2_v01_baltic_v_subarctic_alt_ref_intro <- introgression_plot_v3(scan2_v01_baltic_v_subarctic_dist_alt_ref_20k_df, sample_list=herring_125$sample_list, snp_numbers=scan2_v01_baltic_v_subarctic_alt_ref_diff_SNP_count_vec, pdf_file="outputs/scan2_v01_baltic_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  target_re = target_group)
scan2_v01_baltic_v_subarctic_alt_ref_intro_20k_lists <- complile_intro_lists(intr_obj = scan2_v01_baltic_v_subarctic_alt_ref_intro, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan2_v01_baltic_alt_ref_summary <- summarize_intro(intro_lists = scan2_v01_baltic_v_subarctic_alt_ref_intro_20k_lists, plot_dir = "outputs/scan2_v01_baltic_summary_regions/")
#save(scan2_v01_baltic_alt_ref_summary, file="4.introgression_scan/intermediate_files/scan2_v01_baltic_alt_ref_summary.Rdata")
load(file="4.introgression_scan/intermediate_files/scan2_v01_baltic_alt_ref_summary.Rdata")

## REDUCE REGIONS ####
## 2025-10-18
coverage=7
scan2_v01_baltic_alt_ref_summary_cov7<-scan2_v01_baltic_alt_ref_summary$cov[scan2_v01_baltic_alt_ref_summary$cov[, "cov"] >= coverage,]
scan2_v01_baltic_alt_ref_summary_cov7_gr<-makeGRangesFromDataFrame(df=scan2_v01_baltic_alt_ref_summary_cov7, 
                                                                   seqnames.field = "group_name",
                                                                   keep.extra.columns = T)
# Percentage of introgressioN:
(sum(data.frame(scan2_v01_baltic_alt_ref_summary_cov7_gr)$width)/725670187)*100
0.005512146

scan2_v01_baltic_alt_ref_summary_cov7_gr_min50kb<-GenomicRanges::reduce(scan2_v01_baltic_alt_ref_summary_cov7_gr, min.gapwidth=50000)
(sum(data.frame(scan2_v01_baltic_alt_ref_summary_cov7_gr_min50kb)$width)/725670187)*100
0.008268219


# ---------------------------------------- #
# SCAN3_V01_BALTIC_UK: Parent 1: NSSH_Norway and Canada Atlantic Spring; Parent 2: Pechora +  White sea; Target: UK ####
# ---------------------------------------- #

# Confirm we are selecting the correct individuals:
target_group="Z12_IsleofMan_Atlantic_Autumn|Z14_IsleofMan_Atlantic_Autumn|Z4_IsleofMan_Atlantic_Autumn|CS10_CelticSea_Atlantic_Winter|CS2_CelticSea_Atlantic_Winter|CS4_CelticSea_Atlantic_Winter|CS5_CelticSea_Atlantic_Winter|CS7_CelticSea_Atlantic_Winter|CS8_CelticSea_Atlantic_Winter|AAL1_CelticSea_Atlantic_Winter|AAL2_CelticSea_Atlantic_Winter|AAL3_Celticsea_Atlantic_Winter|AK1_Downs_Atlantic_Winter|AK2_Downs_Atlantic_Winter|AK3_Downs_Atlantic_Winter"; 
ref_1_group = "NSSH|Canada_Atlantic_Spring"; 
ref_2_group = "HWS[2345]"

herring_125$sample_list[grep(target_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]

scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target=target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
save(scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan3_v01_UK_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df))*2

scan3_v01_UK_v_subarctic_alt_ref_intro <- introgression_plot_v3(scan3_v01_UK_v_subarctic_dist_alt_ref_20k_df, sample_list=herring_125$sample_list, snp_numbers=scan3_v01_UK_v_subarctic_alt_ref_diff_SNP_count_vec, pdf_file="outputs/scan3_v01_UK_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  target_re = target_group)
scan3_v01_UK_v_subarctic_alt_ref_intro_20k_lists <- complile_intro_lists(intr_obj = scan3_v01_UK_v_subarctic_alt_ref_intro, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan3_v01_UK_alt_ref_summary <- summarize_intro(intro_lists = scan3_v01_UK_v_subarctic_alt_ref_intro_20k_lists, plot_dir = "outputs/scan3_v01_UK_summary_regions/")
save(scan3_v01_UK_alt_ref_summary, file="4.introgression_scan/intermediate_files/scan3_v01_UK_alt_ref_summary.Rdata")

load("4.introgression_scan/intermediate_files/scan3_v01_UK_alt_ref_summary.Rdata")

## REDUCE REGIONS ####
# 2025-06-19
coverage=7
scan3_v01_UK_alt_ref_summary_cov7<-scan3_v01_UK_alt_ref_summary$cov[scan3_v01_UK_alt_ref_summary$cov[, "cov"] >= coverage,]
scan3_v01_UK_alt_ref_summary_cov7_gr<-makeGRangesFromDataFrame(df=scan3_v01_UK_alt_ref_summary_cov7, 
                                                               seqnames.field = "group_name",
                                                               keep.extra.columns = T)


# Percentage of introgressioN:
(sum(data.frame(scan3_v01_UK_alt_ref_summary_cov7_gr)$width)/725670187)*100
# 0.002756073 %

scan3_v01_UK_alt_ref_summary_cov7_gr_min50kb <- GenomicRanges::reduce(scan3_v01_UK_alt_ref_summary_cov7_gr, min.gapwidth=50000)
(sum(data.frame(scan3_v01_UK_alt_ref_summary_cov7_gr_min50kb)$width)/725670187)*100
# 0.002756073 %

# ---------------------------------------- #
# SCAN4_V01_ATLANTIC_AUTUMN: Parent 1: NSSH_Norway and Canada Atlantic Spring; Parent 2: Pechora +  White sea; Target: Atlantic Autumn ####
# ---------------------------------------- #

# Confirm we are selecting the correct individuals:
target_group="Canada_Atlantic_Autumn"; 
ref_1_group = "NSSH|Canada_Atlantic_Spring"; 
ref_2_group = "HWS[2345]"

herring_125$sample_list[grep(target_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]

scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target=target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
save(scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan4_v01_AtlAut_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df))*2

scan4_v01_AtlAut_v_subarctic_alt_ref_intro <- introgression_plot_v3(scan4_v01_AtlAut_v_subarctic_dist_alt_ref_20k_df, sample_list=herring_125$sample_list, snp_numbers=scan4_v01_AtlAut_v_subarctic_alt_ref_diff_SNP_count_vec, pdf_file="outputs/scan4_v01_AtlAut_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  target_re = target_group)
scan4_v01_AtlAut_v_subarctic_alt_ref_intro_20k_lists <- complile_intro_lists(intr_obj = scan4_v01_AtlAut_v_subarctic_alt_ref_intro, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan4_v01_AtlAut_alt_ref_summary <- summarize_intro(intro_lists = scan4_v01_AtlAut_v_subarctic_alt_ref_intro_20k_lists, plot_dir = "outputs/scan4_v01_AtlAut_summary_regions/")
save(scan4_v01_AtlAut_alt_ref_summary, file="4.introgression_scan/intermediate_files/scan4_v01_AtlAut_alt_ref_summary.Rdata")
load("4.introgression_scan/intermediate_files/scan4_v01_AtlAut_alt_ref_summary.Rdata")

## REDUCE REGIONS ####
# 2025-06-19
coverage=7
scan4_v01_AtlAut_alt_ref_summary_cov7<-scan4_v01_AtlAut_alt_ref_summary$cov[scan4_v01_AtlAut_alt_ref_summary$cov[, "cov"] >= coverage,]
scan4_v01_AtlAut_alt_ref_summary_cov7_gr<-makeGRangesFromDataFrame(df=scan4_v01_AtlAut_alt_ref_summary_cov7, 
                                                              seqnames.field = "group_name",
                                                              keep.extra.columns = T)


# Percentage of introgressioN:
(sum(data.frame(scan4_v01_AtlAut_alt_ref_summary_cov7_gr)$width)/725670187)*100
0 

scan4_v01_AtlAut_alt_ref_summary_cov7_gr_min50kb <- GenomicRanges::reduce(scan4_v01_AtlAut_alt_ref_summary_cov7_gr, min.gapwidth=50000)
(sum(data.frame(scan4_v01_AtlAut_alt_ref_summary_cov7_gr_min50kb)$width)/725670187)*100

# ---------------------------------------- #
# SCAN7_V01_ATLANTIC_SPRING_ALT_REF: Parent 1: Atlantic Autumn; Parent 2: Pechora +  White sea; Target: Atlantic Spring ####
# ---------------------------------------- #

# Confirm we are selecting the correct individuals:
target_group= "NSSH|Canada_Atlantic_Spring"; 
ref_1_group = "Canada_Atlantic_Autumn";  
ref_2_group = "HWS[2345]"

herring_125$sample_list[grep(target_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]

scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target=target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
save(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan7_v01_AtlSpring_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df))*2

scan7_v01_AtlSpring_v_subarctic_alt_ref_intro <- introgression_plot_v3(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df, sample_list=herring_125$sample_list, snp_numbers=scan7_v01_AtlSpring_v_subarctic_alt_ref_diff_SNP_count_vec, pdf_file="outputs/scan7_v01_AtlSpring_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  target_re = target_group)
scan7_v01_AtlSpring_v_subarctic_alt_ref_intro_20k_lists <- complile_intro_lists(intr_obj = scan7_v01_AtlSpring_v_subarctic_alt_ref_intro, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan7_v01_AtlSpring_v_subarctic_alt_ref_summary <- summarize_intro(intro_lists = scan7_v01_AtlSpring_v_subarctic_alt_ref_intro_20k_lists, plot_dir = "outputs/scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_regions/")
save(scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan7_v01_AtlSpring_v_subarctic_dist_alt_ref_20k_df.Rdata")
save(scan7_v01_AtlSpring_v_subarctic_alt_ref_summary, file = "4.introgression_scan/intermediate_files/scan7_v01_AtlSpring_v_subarctic_alt_ref_summary.Rdata")

## REDUCE REGIONS ####
load("4.introgression_scan/intermediate_files/scan7_v01_AtlSpring_v_subarctic_alt_ref_summary.Rdata")
# 2025-06-19
coverage=7
scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7<-scan7_v01_AtlSpring_v_subarctic_alt_ref_summary$cov[scan7_v01_AtlSpring_v_subarctic_alt_ref_summary$cov[, "cov"] >= coverage,]
scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7_gr<-makeGRangesFromDataFrame(df=scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7, 
                                                                                  seqnames.field = "group_name",
                                                                                  keep.extra.columns = T)

# Percentage of introgression:
(sum(data.frame(scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7_gr)$width)/725670187)*100
# 0.05512146 %


scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7_gr_min50kb <- GenomicRanges::reduce(scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7_gr, min.gapwidth=50000)
(sum(data.frame(scan7_v01_AtlSpring_v_subarctic_alt_ref_summary_cov7_gr_min50kb)$width)/725670187)*100
# 0.05512146


# ---------------------------------------- #
# SCAN10_V01_NorthSea: Parent 1: Atlantic Spring; Parent 2: Pechora +  White sea; Target: NorthSea ####
# ---------------------------------------- #

# Confirm we are selecting the correct individuals:
target_group= "NorthSea"
ref_1_group= "NSSH|Canada_Atlantic_Spring"
ref_2_group = "HWS[2345]"

herring_125$sample_list[grep(target_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_1_group, herring_125$sample_list)]
herring_125$sample_list[grep(ref_2_group, herring_125$sample_list)]

scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df <- introgression_contrast(target=target_group, ref_1 = ref_1_group, ref_2 = ref_2_group, w_df=h125_20k_w, sample_list=herring_125$sample_list, geno=herring_125$geno)
save(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df.Rdata")
scan10_v01_NorthSea_v_subarctic_alt_ref_diff_SNP_count_vec <- rowSums(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df[,6:ncol(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df)])/length(6:ncol(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df))*2

scan10_v01_NorthSea_v_subarctic_alt_ref_intro <- introgression_plot_v3(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df, sample_list=herring_125$sample_list, snp_numbers=scan10_v01_NorthSea_v_subarctic_alt_ref_diff_SNP_count_vec, pdf_file="outputs/scan10_v01_NorthSea_v_subarctic_alt_ref_min_diff20_snp_cutoff50.pdf", min_diff = 20, snp_cutoff = 50, assoc_tresh = 8, assoc_dir = "up",  target_re = target_group)
scan10_v01_NorthSea_v_subarctic_alt_ref_intro_20k_lists <- complile_intro_lists(intr_obj = scan10_v01_NorthSea_v_subarctic_alt_ref_intro, win_df = h125_20k_w, fuse_thresh = 5e4, assoc_t = 8)
scan10_v01_NorthSea_alt_ref_summary <- summarize_intro(intro_lists = scan10_v01_NorthSea_v_subarctic_alt_ref_intro_20k_lists, plot_dir = "outputs/scan10_v01_NorthSea_summary_regions/")
save(scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df, file = "4.introgression_scan/intermediate_files/scan10_v01_NorthSea_v_subarctic_dist_alt_ref_20k_df.Rdata")
save(scan10_v01_NorthSea_alt_ref_summary, file = "4.introgression_scan/intermediate_files/scan10_v01_NorthSea_alt_ref_summary.Rdata")
load("4.introgression_scan/intermediate_files/scan10_v01_NorthSea_alt_ref_summary.Rdata")

## REDUCE REGIONS ####
scan10_v01_NorthSea_alt_ref_summary$rec_GR # has Nothing genome range, so let's simply consider this one

coverage=7
scan10_v01_NorthSea_alt_ref_summary_cov7<-scan10_v01_NorthSea_alt_ref_summary$cov[scan10_v01_NorthSea_alt_ref_summary$cov[, "cov"] >= coverage,]
scan10_v01_NorthSea_alt_ref_summary_cov7_gr<-makeGRangesFromDataFrame(df=scan10_v01_NorthSea_alt_ref_summary_cov7, 
                                                                            seqnames.field = "group_name",
                                                                            keep.extra.columns = T)

# Percentage of introgression:
(sum(data.frame(scan10_v01_NorthSea_alt_ref_summary_cov7_gr)$width)/725670187)*100
# 0 %

scan10_v01_NorthSea_alt_ref_summary_cov7_gr_min50kb <- GenomicRanges::reduce(scan10_v01_NorthSea_alt_ref_summary_cov7_gr, min.gapwidth=50000)



#--------------------------------------------- #
# FINDING OVERLAPS ####
#--------------------------------------------- #

makeGRanges

## Scan2 vs Scan1 ####
# Overlap between introgressed regions in Baltic Spring and Baltic Autumn:

scan2_nrow<-nrow(data.frame(scan2_v01_baltic_alt_ref_summary_cov5_gr))

scan2_vs_scan1_nrow<-nrow(data.frame(subsetByOverlaps(scan2_v01_baltic_alt_ref_summary_cov5_gr,scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)))
# GRanges object with 4 ranges and 5 metadata columns:
#   seqnames            ranges strand |     group       cov global_start global_end         col
# <Rle>         <IRanges>  <Rle> | <integer> <integer>    <numeric>  <numeric> <character>
#   92    chr17 25300001-25320000      * |         3         8    517448151  517468150      grey70
# 93    chr17 25320001-25340000      * |         3         6    517468151  517488150      grey70
# 94    chr17 25340001-25360000      * |         3         7    517488151  517508150      grey70
# 146    chr16 14120001-14140000      * |         9         6    478494329  478514328      grey30
# -------
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

scan2_vs_scan1_nrow/ scan2_nrow
# 1

# Overlap between introgressed regions in Baltic Spring and UK:
scan3_nrow<-nrow(data.frame(scan3_v01_UK_alt_ref_summary_cov5_gr))

scan3_vs_scan1_nrow<-nrow(data.frame(subsetByOverlaps(scan3_v01_UK_alt_ref_summary_cov5_gr,scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)))

scan3_vs_scan1_nrow/ scan3_nrow
# 0

# Overlap between introgressed regions in Baltic Spring and Atlantic Autumn:
#. No regions in scan 4
#scan4_nrow<-nrow(data.frame(scan4_v01_AtlAut_alt_ref_summary_gr))

#scan4_vs_scan1_nrow<-nrow(data.frame(subsetByOverlaps(scan4_v01_AtlAut_alt_ref_summary_gr,scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)))

#scan4_vs_scan1_nrow/ scan4_nrow
#0 

# Overlap between introgressed regions in Baltic Spring and Atlantic Spring:
scan5_nrow<-nrow(data.frame(scan5_v01_AtlSpring_alt_ref_summary_cov5_gr))

scan5_vs_scan1_nrow<-nrow(data.frame(subsetByOverlaps(scan5_v01_AtlSpring_alt_ref_summary_cov5_gr,scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)))

scan5_vs_scan1_nrow/ scan5_nrow
#0 

# ------------------------------ #
# PLOT REGIONS ####
# ------------------------------ #

herring_txdb<-makeTxDbFromEnsembl(organism="Clupea harengus",
                                  release=NA,
                                  circ_seqs=NULL,
                                  server="ensembldb.ensembl.org",
                                  username="anonymous", password=NULL, port=0L,
                                  tx_attrib=NA)

# retrieve the annotation from biomart
bm <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="charengus_gene_ensembl")
ensembl = useEnsembl(biomart="ensembl",dataset="charengus_gene_ensembl")
attributes<-listAttributes(ensembl)
attr = c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id","wikigene_id","entrezgene_id","chromosome_name","start_position",
         "end_position",
         "description")

regions <- getBM(attributes = attr, mart = bm)

# you can also get a TxDB from biomart
herring_txdb_biomart<-makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "charengus_gene_ensembl")

## Scan 1 ####
scan1_v01_baltic_alt_ref_summary_filter2_cov_df<-scan1_v01_baltic_alt_ref_summary_filter2$cov
scan1_v01_baltic_alt_ref_summary_filter2_cov_df$chromosome<-str_remove(scan1_v01_baltic_alt_ref_summary_filter2_cov_df$group_name, "chr")
scan1_v01_baltic_alt_ref_summary_filter2_cov_df$chromosome<-as.numeric(scan1_v01_baltic_alt_ref_summary_filter2_cov_df$chromosome)
scan1_v01_baltic_alt_ref_summary_filter2_cov_df<-scan1_v01_baltic_alt_ref_summary_filter2_cov_df[order(scan1_v01_baltic_alt_ref_summary_filter2_cov_df$chromosome, scan1_v01_baltic_alt_ref_summary_filter2_cov_df$start, scan1_v01_baltic_alt_ref_summary_filter2_cov_df$end),]
scan1_v01_baltic_alt_ref_summary_filter2_cov_df$n <- 1:nrow(scan1_v01_baltic_alt_ref_summary_filter2_cov_df)

axisdf <- scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% group_by(chromosome) %>% summarize(center=(max(n) + min(n) ) / 2 )
chroms <- data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% group_by(chromosome) %>% summarize(end=max(n), start= min(n)))
chroms$chromosome<-as.character(chroms$chromosome)

ggplot(scan1_v01_baltic_alt_ref_summary_filter2_cov_df)+
  #geom_rect(data = chroms, mapping=aes(xmin = start, xmax = end, ymin = -10, ymax = 10) , color="gray", fill=NA, alpha=0.2)+
  geom_point(aes(x=n, y=cov, color=as.factor(chromosome)))+
  scale_color_manual(values = rep(c("black","gray"), 23 )) +
  scale_x_continuous(label = axisdf$chromosome, breaks= axisdf$center , expand = c(0.01, 0.01))+
  geom_hline(yintercept=7, color="purple")+
  xlab("Chromosomes")+
  ylab("Number of haplotypes")+
  theme_classic()+
  theme(axis.text.x=element_text())

axisdf <- scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% group_by(group_name) %>% summarize(center=(max(global_start) + min(global_start) ) / 2 )
chroms <- data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% group_by(group_name) %>% summarize(end=max(global_start), start= min(global_start)))
chroms$chromosome<-as.numeric(str_remove(chroms$group_name, "chr"))
chroms<-chroms[order(chroms$chromosome),]

plot_coverage<-ggplot(scan1_v01_baltic_alt_ref_summary_filter2_cov_df)+
  #geom_rect(data = chroms, mapping=aes(xmin = start, xmax = end, ymin = 0, ymax = 25, fill=group_name) , alpha=0.2)+
  #scale_fill_manual(values = rep(c("black","gray"), 23 ))+
  geom_point(aes(x=(global_start+global_end)/2, y=cov,color=group_name))+
  #scale_color_manual(values = rep(c("black","gray"), 23 ))+
  scale_x_continuous(label = str_remove(axisdf$group_name, "chr"), breaks= axisdf$center , expand = c(0.01, 0.01))+
  geom_hline(yintercept=7, color="black", linetype="dashed")+
  xlab("Chromosomes")+
  ylab("Number of haplotypes")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=0, hjust=0), legend.position = "none")


ggsave(plot_coverage, file="figures/scan1_v01_baltic_alt_ref_summary_filter2_cov_df.coverage.png", height=3, width=10, dpi=300)

# Plot chromosome by chromosome

chr17_coverage<-scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% 
  filter(group_name=="chr17") %>%
  ggplot()+
  geom_point(aes(x=(start+end)/2e6, y=cov), color="#97A900")+
  theme_classic()+
  geom_hline(yintercept=7, color="black", linetype="dashed")+
  xlab("Chr17 Position (Mb)")+
  ylab("No. of haplotypes")

chr10_coverage<-scan1_v01_baltic_alt_ref_summary_filter2_cov_df %>% 
  filter(group_name=="chr10") %>%
  ggplot()+
  geom_point(aes(x=(start+end)/2e6, y=cov), color="#DD8D00")+
  theme_classic()+
  geom_hline(yintercept=7, color="black", linetype="dashed")+
  xlab("Chr10 Position (Mp)")+
  ylab("No. of haplotypes")

chr10_chr17_coverage<-gridExtra::grid.arrange(chr10_coverage, chr17_coverage, nrow=1)
ggsave(chr10_chr17_coverage, file="figures/scan1_v01_baltic_alt_ref_summary_filter2_cov_df_chr10_chr17.png", width=10, height = 3, dpi=300)

write.table(scan1_v01_baltic_alt_ref_summary_filter2_cov_df, file="introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2.coverage.txt", col.names = T, row.names = F, quote =F, sep="\t")

intro_chromosomes<-unique(scan1_v01_baltic_alt_ref_summary_filter2_cov_df$group_name)
colors_chromosomes<-hue_pal()(25)

## Annotation ####
seqlevelsStyle(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb) <-"NCBI"

# Exact overlap:
scan1_genes<-subsetByOverlaps(genes(herring_txdb_biomart), scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb)
scan1_match<-regions[match(scan1_genes$gene_id, regions$ensembl_gene_id),c(1,2, 5, 6, 7,8,9)]
scan1_match_gr<-GRanges(seqnames=scan1_match$chromosome_name, 
                        ranges = IRanges(start = scan1_match$start_position, end = scan1_match$end_position),
                        ensembl_gene_id = scan1_match$ensembl_gene_id, external_gene_name = scan1_match$external_gene_name,
                        description = scan1_match$description)
scan1_match_df <- data.frame(scan1_match_gr)
write.table(scan1_match_df, file="introgression_regions_annotations/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.nogap.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# Overlap + 20 kb
overlaps <- findOverlaps(query = genes(herring_txdb_biomart), scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb, maxgap = 20000)

scan1_genes<-subsetByOverlaps(genes(herring_txdb_biomart), scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb, maxgap = 20000)
scan1_match<-regions[match(scan1_genes$gene_id, regions$ensembl_gene_id),c(1,2, 5, 6, 7,8,9)]
scan1_match_gr<-GRanges(seqnames=scan1_match$chromosome_name, 
                        ranges = IRanges(start = scan1_match$start_position, end = scan1_match$end_position),
                        ensembl_gene_id = scan1_match$ensembl_gene_id, external_gene_name = scan1_match$external_gene_name,
                        description = scan1_match$description)
scan1_match_df <- data.frame(scan1_match_gr)
write.table(scan1_match_df, file="introgression_regions_annotations/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.txt", quote = F, col.names = T, row.names = F, sep = "\t")

scan1_match_df<-read.table("4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.modified.txt", sep="\t", header=T, row.names=NULL)
scan1_match_gr <- makeGRangesFromDataFrame(scan1_match_df)

scan1_anno_overl <- findOverlaps(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb, scan1_match_gr, maxgap = 2e4)  

scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df <- data.frame(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb[scan1_anno_overl@from,])
scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$genes <- NA
scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$genes <- scan1_match_df[scan1_anno_overl@to,]$external_gene_name

scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$names <- paste0(scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$seqnames, "_", scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$start, "_", scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df$end)

supplementary_table3 <- scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb_df %>% 
  group_by(names) %>% 
  summarise(Annotations = paste(unique(genes), collapse =", ")) %>%
  print(n=100)

# Results present in Supplementary Table 4
write.table(supplementary_table3, "introgression_regions_annotations/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.asTableS3.txt", col.names = T, row.names = F, quote = F, sep = "\t")


### FUNCTIONS ####
calc_cov <- function(GR_list, size_df = Ch_v2.0.2_sizes){
  cov_obj <- coverage(unlist(GR_list))
  cov_df <- as.data.frame(ranges(cov_obj))
  cov_df[,"cov"] <- unlist(runValue(cov_obj))
  cov_df[,"global_start"] <- cov_df[,"start"] + size_df[match(cov_df[, "group_name"], size_df[,"name"]), "offset"]
  cov_df[,"global_end"] <-  cov_df[,"end"] + size_df[match(cov_df[, "group_name"], size_df[,"name"]), "offset"]
  cov_df[,"col"] <-  size_df[match(cov_df[, "group_name"], size_df[,"name"]), "col"]
  return(cov_df)
}

#adjusted for more general group names
introgression_contrast <- function(target = "HWS-6", ref_1 = "A[MF]", ref_2 = "HWS-1", w_df, sample_list, geno){
  #target_haps <- grep("HWS-6", h67_clean_samples)
  #ref_atl_haps <- grep("A[MF]", h67_clean_samples)
  #ref_pac_haps <- grep("Paci", h67_clean_samples)
  #ref_pac_haps <- grep("HWS-1", h67_clean_samples)
  target_haps <- grep(target, sample_list)
  ref_1_haps <- grep(ref_1, sample_list)
  ref_2_haps <- grep(ref_2, sample_list)
  
  #hws6_dist_df <- h67_w
  dist_df <- w_df
  
  for (w_idx in 1:dim(dist_df)[1]){
    win <- dist_df[w_idx,]
    scaff_geno <- geno[geno[,1] ==  win[1,1],]
    target_markers <- scaff_geno[,2] >= win[1,2] & scaff_geno[,2] < win[1,3]
    target_geno <- scaff_geno[target_markers,]
    if(dim(target_geno)[1] > 0){
      n_haps <- nchar(target_geno[1,3])
      hap_vec <- array()
      for(i in 1:n_haps){
        hap_vec[i] <- paste(substr(target_geno[,3], i, i), collapse = "")
      }
      hap_dist_mat <- outer(hap_vec, hap_vec, FUN = "stringdist", method = "hamming")
      diag(hap_dist_mat) <- NA
    } else{
      print(paste("Win", w_idx, "empty; skipping."))
    }
    for(th in target_haps){
      #hws6_dist_df[w_idx, paste(h67_clean_samples[th],"_atl", sep = "")] <- mean(hap_dist_mat[th, ref_atl_haps])
      #hws6_dist_df[w_idx, paste(h67_clean_samples[th],"_pac", sep = "")] <- mean(hap_dist_mat[th, ref_pac_haps])
      dist_df[w_idx, paste(sample_list[th],"_1", sep = "")] <- min(hap_dist_mat[th, ref_1_haps])
      dist_df[w_idx, paste(sample_list[th],"_2", sep = "")] <- min(hap_dist_mat[th, ref_2_haps])
    }
  }
  return(dist_df)
}
introgression_plot_v3 <- function(dist_df, sample_list,  snp_numbers, assoc_tresh = 100, snp_cutoff = 200, min_diff = 100, pdf_file = "", eps = 1e-3, assoc_dir = "both", target_re = "[0-9][._].+"){
  #target <- sub(target_re, "", names(dist_df)[6])
  if (pdf_file != "") pdf(file = pdf_file, width = 15, height = 5)
  #target_haps <- grep(target, sample_list)
  target_haps <- grep(target_re, sample_list)
  plot(x= 1, y=1, xlim = c(0,max(dist_df[,4])), ylim = c(log10(eps)-3, -log10(eps)+3), type = "n", main = paste(target_re," vs Subarctic C. Pallasi", sep = ""), xlab = "Cumulative position", ylab = "Log10 of distance ratio")
  hap_ratio_df <- data.frame(hap = sample_list[target_haps], stringsAsFactors= F)
  #assoc_tresh <- 100
  #snp_cutoff <- 200
  #min_diff <- 100
  
  for(bf_hap in 1:dim(hap_ratio_df)[1]){
    bf_ind <- ceiling(bf_hap/2)
    hap_no <- 2 - bf_hap%%2
    hap_ratio_df[bf_hap, "ref_1_mean"] <- mean(dist_df[,2*(bf_hap-1)+6])
    hap_ratio_df[bf_hap, "ref_2_mean"] <- mean(dist_df[,2*(bf_hap-1)+7])
    mean_ratio <- hap_ratio_df[bf_hap, "ref_1_mean"]/hap_ratio_df[bf_hap, "ref_2_mean"]
    
    rv <- (dist_df[,2*(bf_hap-1)+6]+eps)/(dist_df[,2*(bf_hap-1)+7]+eps)
    rv[snp_numbers <= snp_cutoff|(dist_df[,2*(bf_hap-1)+6] < min_diff & dist_df[,2*(bf_hap-1)+7] < min_diff)] <- NA #1
    
    dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")] <- rv
    if(assoc_dir == "both"){
      assoc_filter <- which(rv > assoc_tresh * mean_ratio | rv < mean_ratio/assoc_tresh) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    if(assoc_dir == "up"){
      assoc_filter <- which(rv > assoc_tresh * mean_ratio ) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    if(assoc_dir == "down") {
      assoc_filter <- which(rv < mean_ratio/assoc_tresh) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
    }
    
    
    hap_ratio_df[bf_hap, "pac"] <- sum(rv > mean_ratio*assoc_tresh, na.rm = T)
    hap_ratio_df[bf_hap, "atl"] <- sum(rv < mean_ratio/assoc_tresh, na.rm = T)
    bg_col_vec <- c("grey20", "grey80")[match(as.integer(as.factor(dist_df[, 1])), unique(as.integer(as.factor(dist_df[, 1]))))%%2 + 1]
    if(assoc_dir != "none") {
      points(y = log10(dist_df[-assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[-assoc_filter,4], col = bg_col_vec[-assoc_filter], pch = 16, cex=0.2)
      points(y = log10(dist_df[assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[assoc_filter,4], col = rgb(t(col2rgb(bf_ind+1)), alpha = 120, maxColorValue = 255), pch = 15+hap_no, cex = 0.8)
    }
    if(assoc_dir == "none") {
      points(y = log10(dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")]), x = dist_df[,4], col = bg_col_vec, pch = 16, cex=0.2)
    }
  }
  if(assoc_dir == "none") {
    abline(h = log10(c(mean_ratio/assoc_tresh, assoc_tresh * mean_ratio)), lty = "dashed", col = "red", lwd = 2)
  }
  if (pdf_file != "")dev.off()
  return(list(dist = dist_df, hap= hap_ratio_df))
}
complile_intro_lists <- function(intr_obj, win_df,  fuse_thresh = 2e5, assoc_t = 5){
  ratio_cols <- grep("ratio", names(intr_obj$dist))
  #GRanges-based quantifications
  a_GR_list <- GRangesList()
  a_GR_red_list <- GRangesList()
  p_GR_list <- GRangesList()
  p_GR_red_list <- GRangesList()
  
  for (i in 1:length(ratio_cols)) {
    hap_name <- sub("_ratio", "", names(intr_obj$dist)[ratio_cols[i]])
    dist_vec <- intr_obj$dist[, ratio_cols[i]]
    a_win <- which(dist_vec <= median(dist_vec, na.rm = T) / assoc_t)
    a_GR <-
      GRanges(
        seqnames = win_df$scaffold[a_win],
        ranges =  IRanges(start = win_df$start[a_win], end = win_df$stop[a_win])
      )
    p_win <- which(dist_vec >= median(dist_vec, na.rm = T) * assoc_t)
    p_GR <-
      GRanges(
        seqnames = win_df$scaffold[p_win],
        ranges =  IRanges(start = win_df$start[p_win], end = win_df$stop[p_win])
      )
    a_GR_list[[hap_name]] <- a_GR
    a_GR_red_list[[hap_name]] <- GenomicRanges::reduce(a_GR, min.gapwidth = fuse_thresh)
    p_GR_list[[hap_name]] <- p_GR
    p_GR_red_list[[hap_name]] <- GenomicRanges::reduce(p_GR, min.gapwidth = fuse_thresh)
  }
  return(list(atl = a_GR_list, atl_red = a_GR_red_list, pac = p_GR_list, pac_red = p_GR_red_list))
}

plot_top_regions_2 <- function(top_pos, win_df, geno, sample_list = NULL, plot_dir = "", snp_col = "orange", snp_lwd = 0.1, snp_alpha = 50){
  require(stringdist)
  require(diptest)
  hclust_list <- list()
  for (pos in top_pos){
    win <- win_df[pos,]
    scaff_geno <- geno[geno[,1] ==  win[1,1],]
    target_markers <- scaff_geno[,2] >= win[1,2] & scaff_geno[,2] < win[1,3]
    target_geno <- scaff_geno[target_markers,]
    n_haps <- nchar(target_geno[1,3])
    if(is.na(n_haps)) next # added this...
    hap_vec <- array()
    for(i in 1:n_haps){
      hap_vec[i] <- paste(substr(target_geno[,3], i, i), collapse = "")
    }
    hap_dist_mat <- outer(hap_vec, hap_vec, FUN = "stringdist", method = "hamming")
    diag(hap_dist_mat) <- NA
    hap_dist_data <- array(data = hap_dist_mat[!is.na(hap_dist_mat)])
    
    if (is.null(sample_list)) sample_list <- 1:n_haps
    rownames(hap_dist_mat) <- sample_list
    colnames(hap_dist_mat) <- sample_list
    win_loc <- paste(win[1], "_", round(win[2]/1000), "_to_", round(win[3]/1000), "_kb", sep = "")
    pdf(paste(plot_dir, win_loc, ".pdf", sep = ""), width = 14, height = 10, paper = "a4r")
    
    #Distance distribution histogram
    hist(hap_dist_data, main = win_loc, xlab = "Edit distance")
    
    #Dendrogram
    hap_hc <- hclust(as.dist(hap_dist_mat))
    plot(hap_hc, xlab ="", ann = F, cex = min(c(1, 60/length(sample_list))))
    title(main = win_loc)
    
    #Heatmap
    block_hm <- heatmap(x= hap_dist_mat, scale = "none", margins = c(4,4), Rowv = as.dendrogram(hap_hc), Colv = as.dendrogram(hap_hc), cexRow = min(c(1, 35/length(sample_list))), labCol = "")
    
    #Haplotype visualisation
    ext_start <- max(pos-5, min(which(win_df[,1] == win[1,1])))
    ext_stop <- min(pos+5, max(which(win_df[,1] == win[1,1])))
    ext_target_m <- scaff_geno[,2] >= win_df[ext_start,2] & scaff_geno[,2] < win_df[ext_stop,3]
    ext_target_g <- scaff_geno[ext_target_m,]
    for(i in 1:n_haps){
      hap_vec[i] <- paste(substr(ext_target_g[,3], i, i), collapse = "")
    }
    plot(1:n_haps, type = "n", xlim = range(ext_target_g[,2]), axes = F, ann = F)
    axis(2, labels = sample_list[block_hm$rowInd], at = 1:n_haps, las = 2, cex.axis= min(c(0.6, 20/n_haps)), col=1:2, pos = range(ext_target_g[,2])[1] - 0.0175*diff(range(ext_target_g[,2])))
    axis(1)
    title(main = win_loc)
    rect(range(ext_target_g[,2])[1]- 0.0025*diff(range(ext_target_g[,2])),0,range(ext_target_g[,2])[2] + 0.0025*diff(range(ext_target_g[,2])), n_haps+1, col = "grey90")
    rect(range(target_geno[,2])[1],0.1,range(target_geno[,2])[2], n_haps+0.9, col = "grey80", border = NA)
    #rect(range(ext_target_g[,2])[1] - 0.015*diff(range(ext_target_g[,2])), which(block_hm$rowInd == 1)-0.4, range(ext_target_g[,2])[1] - 0.005*diff(range(ext_target_g[,2])), which(block_hm$rowInd == 1) +0.4, col = "orange")
    #ref_col <- factor(unlist(strsplit(hap_vec[1], "")), levels = c("G", "C", "T", "A"))
    allele_df <- data.frame("G" = 1:dim(ext_target_g)[1], "C" = 0, "T" = 0, "A" = 0, "major" = NA, stringsAsFactors = F)
    allele_df[,"G"] <- nchar(gsub("[^G]", "",ext_target_g[,3]))
    allele_df[,"C"] <- nchar(gsub("[^C]", "",ext_target_g[,3]))
    allele_df[,"T"] <- nchar(gsub("[^T]", "",ext_target_g[,3]))
    allele_df[,"A"] <- nchar(gsub("[^A]", "",ext_target_g[,3]))			
    allele_max_vec <- pmax(allele_df[,"G"], allele_df[,"C"], allele_df[,"T"], allele_df[,"A"])
    allele_df[allele_df[,"G"] == allele_max_vec,"major"] <- "G"
    allele_df[allele_df[,"C"] == allele_max_vec,"major"] <- "C"
    allele_df[allele_df[,"T"] == allele_max_vec,"major"] <- "T"
    allele_df[allele_df[,"A"] == allele_max_vec,"major"] <- "A"
    ref_col <- factor(allele_df[,"major"], levels = c("G", "C", "T", "A"))
    for(i in 1:n_haps){
      tmp_col <- factor(unlist(strsplit(hap_vec[i], "")), levels = c("G", "C", "T", "A"))
      y_vec <- rep(which(block_hm$rowInd == i), length(tmp_col))[tmp_col != ref_col]
      x_vec <- ext_target_g[,2][tmp_col != ref_col]
      #rect(x_vec, y_vec-0.4 , x_vec + 0.001*diff(range(ext_target_g[,2])), y_vec+0.4, col = rgb(t(col2rgb("orange")), alpha = 100, max= 255), border = NA)
      segments(x0 = x_vec, y0 = y_vec-0.4, y1 = y_vec+0.4, col = rgb(t(col2rgb(snp_col)),alpha = min(c(snp_alpha,100*(10000/length(tmp_col)))) , max= 255), lwd = snp_lwd)
      #alpha = min(c(50,100*(1000/length(x_vec))))
    }
    segments(x0 = range(target_geno[,2]), y0 = 0.05, y1 = n_haps + 0.95)
    dev.off()
    hclust_list[[win_loc]] <- hap_hc
  }
  return(hclust_list)
}

summarize_intro <- function(intro_lists, target_list = "pac_red", plot_dir = "~/Projects/Herring/doc/Baltic_v_Pac/HapDist/", win_df = h125_20k_w, HapDistObj = herring_125, rec_co = 5){
  tmp_list <- intro_lists[[target_list]]
  tmp_intro_cov <- calc_cov(tmp_list)
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  if(!dir.exists(paste0(plot_dir, "/any_intro"))) dir.create(paste0(plot_dir, "/any_intro/"))
  
  tmp_recurring_intro_cov_GR <- NULL
  if(any(tmp_intro_cov$cov >= rec_co)){
    tmp_recurring_intro_cov_df <- tmp_intro_cov[tmp_intro_cov$cov >= rec_co,c(2:4)]
    names(tmp_recurring_intro_cov_df )[1] <- "seqnames"
    tmp_recurring_intro_cov_GR <- GRanges(tmp_recurring_intro_cov_df)
  }
  
  tmp_all_intro_cov_GR <- NULL
  if(any(tmp_intro_cov$cov >= 1)){
    tmp_all_intro_cov_df <- tmp_intro_cov[tmp_intro_cov$cov >= 1,c(2:4)]
    names(tmp_all_intro_cov_df )[1] <- "seqnames"
    tmp_all_intro_cov_GR <- GRanges(tmp_all_intro_cov_df)
  }
  tmp_win_GR <- GRanges(seqnames = win_df$scaffold, ranges = IRanges(start = win_df$start, end =  win_df$stop))
  
  if(!is.null(tmp_recurring_intro_cov_GR)){
    tmp_rec_win_hits <- findOverlaps(tmp_win_GR, tmp_recurring_intro_cov_GR) #Using reduced version
    plot_top_regions_2(top_pos = tmp_rec_win_hits@from, win_df = win_df, geno = HapDistObj$geno, sample_list = HapDistObj$sample_list, plot_dir = plot_dir, snp_col = "blue")
  }
  
  #all hits
  if(!is.null(tmp_all_intro_cov_GR)){
    tmp_all_win_hits <- findOverlaps(tmp_win_GR, tmp_all_intro_cov_GR) #Using reduced version
    plot_top_regions_2(top_pos = tmp_all_win_hits@from, win_df = win_df, geno = HapDistObj$geno, sample_list = HapDistObj$sample_list, plot_dir = paste0(plot_dir, "/any_intro/"), snp_col = "blue")
  }
  return(list(cov = tmp_intro_cov, rec_GR = tmp_recurring_intro_cov_GR, all_GR = tmp_all_intro_cov_GR))
}

