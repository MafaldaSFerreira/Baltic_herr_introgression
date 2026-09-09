install.packages("poolfstat")
library(poolfstat)
library(tidyverse)
#library(separate) # or use base R as below

# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
setwd("../Figshare/5.processing_pool_seq_data/")

# Recalculate DAFs with Han et al Populations and AD ####
## Populations From Han Et al #### 
Baltic_Spring_col <- c("A_Kalix_Baltic_Spring", 
                       "B_Vaxholm_Baltic_Spring", 
                       "G_Gamleby_Baltic_Spring", 
                       "HGS1_Riga_Baltic_Spring", 
                       "HGS2_Riga_Baltic_Spring", 
                       "PB11_Kalmar_Baltic_Spring",
                       "PB12_Karlskrona_Baltic_Spring", 
                       "PB1_HastKar_Baltic_Spring", 
                       "PB4_Hudiksvall_Baltic_Spring",
                       "PB5_Galve_Baltic_Spring",
                       "PN3_CentralBaltic_Baltic_Spring",
                       "PB6_Galve_Baltic_Summer")

Atlantic_Spring_col <- c("HGS15_NSSH_Atlantic_Spring", 
                         "PB2_Iceland_Atlantic_Spring", 
                         "Q_Norway_Atlantic_Atlantic_Spring",
                         "PB9_Kattegat_Atlantic_Spring",
                         "PB10_Skagerrak_Atlantic_Spring",
                         "O_Hamburgsund_Atlantic_Spring",
                         "HGS8_KattegatNorth_Atlantic_Spring")

Baltic_Autumn_col <- c("HGS12_BornholmBasin_Baltic_Autumn",
                       "HGS3_Riga_Baltic_Autumn",
                       "HGS4_Riga_Baltic_Autumn",
                       "H_Fehmarn_Baltic_Autumn",
                       "PB7_Galve_Baltic_Autumn")

Atlantic_Autumn_col <- c("HGS16_Orkney_NorthSea_Autumn",
                         "HGS17_IsleOfMan_IrishSea_Autumn",
                         "HGS22_CapeWrath_Atlantic_Autumn",
                         "N_NorthSea_Atlantic_Autumn")


# READ ALLELE COUNTS ####
AD<-read.table("inputs/60.Neff.AD", header=T)
names <- names(AD)


# Combine all your population vectors into one target list
target_pops <- c(Baltic_Spring_col, Atlantic_Spring_col, 
                 Baltic_Autumn_col, Atlantic_Autumn_col)

# Subset AD: keep Chrom (col 1), Pos (col 2), and any column name in our target list
# Note: Check if your AD file actually has "Chrom" and "Pos" as the first two column names
AD_filtered <- AD[, c(names(AD)[1:2], intersect(names(AD), target_pops))]


# Identify population columns in the new subsetted data
pop_idx <- 3:ncol(AD_filtered)

# Process row-wise filtering
keep_indices <- apply(AD_filtered[, pop_idx], 1, function(row) {
  # Split "Ref,Alt" into a numeric vector
  counts <- as.numeric(unlist(strsplit(as.character(row), ",")))
  
  # Sum all Ref alleles (odd positions) and Alt alleles (even positions)
  total_ref <- sum(counts[seq(1, length(counts), by = 2)], na.rm = TRUE)
  total_alt <- sum(counts[seq(2, length(counts), by = 2)], na.rm = TRUE)
  global_total <- total_ref + total_alt
  
  # 1. Check if biallelic (both alleles must be present at least once)
  if (total_ref == 0 | total_alt == 0) return(FALSE)
  
  # 2. Calculate MAF
  maf <- min(total_ref, total_alt) / global_total
  
  # Return TRUE if MAF is at least 0.01
  return(maf >= 0.01)
})

# Apply the row filter
AD_final <- AD_filtered[keep_indices, ]

# A) The SNP position matrix
snp_pos_matrix <- as.matrix(AD_final[, 1:2])

# B) The TreeMix allele count file (Populations only)
treemix_data <- AD_final[, 3:ncol(AD_final)]

# Write to disk
# save(snp_pos_matrix, file="results/filtered_han_et_al.snp_pos_matrix")
# 
# output_file <- gzfile("results/filtered_han_et_al.treemix.gz", "w")
# write.table(treemix_data, 
#             file = output_file, 
#             quote = FALSE, 
#             row.names = FALSE, 
#             sep = " ")
# close(output_file)

# C) Read into poolfstat
pooldata <- genotreemix2countdata(
  genotreemix.file = "results/filtered_han_et_al.treemix.gz",
  snp.pos = snp_pos_matrix
)

save(pooldata, file="fst_results/pooldata_for_poolfstat.RData")
load("fst_results/pooldata_for_poolfstat.RData")

# # 1. Get the population names in the order they exist in your pooldata object
# ordered_pops <- pooldata@popnames # or pooldata$popnames depending on object type
# 
# # 2. Create the struct vector by checking which group each population belongs to
# struct_vector <- ifelse(ordered_pops %in% Baltic_Spring_col, "Baltic_Spring",
#                         ifelse(ordered_pops %in% Atlantic_Spring_col, "Atlantic_Spring",
#                                ifelse(ordered_pops %in% Baltic_Autumn_col, "Baltic_Autumn",
#                                       ifelse(ordered_pops %in% Atlantic_Autumn_col, "Atlantic_Autumn", "Unknown"))))
# 
# 
# # 3. Check if any populations were missed
# if(any(struct_vector == "Unknown")) {
#   warning("Some populations were not assigned to a group!")
# }
# 

# # Estimate Fst 
# # Lets test with a big window size
# fst_results_1e5 <- computeFST(
#   pooldata, 
#   method = "Anova", 
#   struct = struct_vector,
#   sliding.window.size = 1e5,
# )
# 
# fst_results_1 <- computeFST(
#   pooldata, 
#   method = "Anova", 
#   struct = struct_vector,
#   sliding.window.size = 1,
# )
# 
# # Extract per-SNP estimates
# snp_fst_df <- fst_results_1$snp.Fstats
# 
# # Combine with SNP info (Chr and Pos) for plotting
# final_snp_results <- cbind(pooldata@snp.info[, 1:2], snp_fst_df)
# 
# 
# head(final_snp_results)
# 
# final_snp_results %>% filter(Chromosome=="chr1") %>%
#   ggplot() + 
#   geom_point(aes(x=Position, y=Fst))
# 

# I need to create a subset of the pooldata so I can contrast spring and autumn. 

# Find indices for Spring populations
spring_idx <- which(pooldata@popnames %in% c(Baltic_Spring_col, Atlantic_Spring_col))

# Create the subset with a MAF filter of 0.01
countdata_spring_maf0.01 <- countdata.subset(
  pooldata, 
  pop.index = spring_idx, 
  min.maf = 0.01
)

# Create the subset with a MAF filter of 0.05
countdata_spring_maf0.05 <- countdata.subset(
  pooldata, 
  pop.index = spring_idx, 
  min.maf = 0.05
)

# Define the structure for the two groups
struct_spring <- ifelse(countdata_spring_maf0.05@popnames %in% Baltic_Spring_col, "Baltic", "Atlantic")

# Run FST
fst_spring_1SNP_maf0.01 <- computeFST(countdata_spring_maf0.01, 
                                      method = "Anova", 
                                      struct = struct_spring,
                                      sliding.window.size = 1)

fst_spring_1SNP_maf0.05 <- computeFST(countdata_spring_maf0.05, 
                         method = "Anova", 
                         struct = struct_spring,
                         sliding.window.size = 1)



# Find indices for Autumn populations
autumn_idx <- which(pooldata@popnames %in% c(Baltic_Autumn_col, Atlantic_Autumn_col))

# Create the subset with a MAF filter of 0.01
countdata_autumn_maf0.01 <- countdata.subset(
  pooldata, 
  pop.index = autumn_idx, 
  min.maf = 0.01
)

# Create the subset with a MAF filter of 0.05
countdata_autumn_maf0.05 <- countdata.subset(
  pooldata, 
  pop.index = autumn_idx, 
  min.maf = 0.05
)


# Define the structure for the two groups
struct_autumn <- ifelse(countdata_autumn_maf0.05@popnames %in% Baltic_Autumn_col, "Baltic", "Atlantic")

# Run FST
fst_autumn_1SNP_maf0.01 <- computeFST(countdata_autumn_maf0.01, 
                                      method = "Anova", 
                                      struct = struct_autumn,
                                      sliding.window.size = 1)


fst_autumn_1SNP_maf0.05 <- computeFST(countdata_autumn_maf0.05, 
                           method = "Anova", 
                           struct = struct_autumn,
                           sliding.window.size = 1)

# Maf 0.01:
# Per-SNP Fst for Spring
df_spring_maf0.01 <- data.frame(
  Chr = countdata_spring_maf0.01@snp.info[,1],
  Pos = as.numeric(countdata_spring_maf0.01@snp.info[,2]),
  FST = fst_spring_1SNP_maf0.01$snp.Fstats$Fst
)

# Per-SNP Fst for Autumn
df_autumn_maf0.01 <- data.frame(
  Chr = countdata_autumn_maf0.01@snp.info[,1],
  Pos = as.numeric(countdata_autumn_maf0.01@snp.info[,2]),
  FST = fst_autumn_1SNP_maf0.01$snp.Fstats$Fst
)


# Maf 0.05:
# Per-SNP Fst for Spring
df_spring_maf0.05 <- data.frame(
  Chr = countdata_spring_maf0.05@snp.info[,1],
  Pos = as.numeric(countdata_spring_maf0.05@snp.info[,2]),
  FST = fst_spring_1SNP_maf0.05$snp.Fstats$Fst
)

# Per-SNP Fst for Autumn
df_autumn_maf0.05 <- data.frame(
  Chr = countdata_autumn_maf0.05@snp.info[,1],
  Pos = as.numeric(countdata_autumn_maf0.05@snp.info[,2]),
  FST = fst_autumn_1SNP_maf0.05$snp.Fstats$Fst
)


# Plotting:
plotmaf0.05 <- df_spring_maf0.05 %>% 
  filter(Chr=="chr10") %>%
  filter(FST > 0) %>% 
  ggplot()+
  geom_point(aes(x=Pos, y=FST))

plotmaf0.01 <- df_spring_maf0.01 %>% 
  filter(Chr=="chr10") %>%
  filter(FST > 0) %>% 
  ggplot()+
  geom_point(aes(x=Pos, y=FST))

df_spring_maf0.01 %>% filter(FST >=0 & !is.na(FST)) %>%
  ggplot() + geom_histogram(aes(x=FST), bins=50)

df_spring_maf0.05 %>% filter(FST >=0 & !is.na(FST)) %>%
  ggplot() + geom_histogram(aes(x=FST), bins=50)

# Subset:
df_spring_maf0.05_filtered <- df_spring_maf0.05 %>% filter(FST>=0 & !is.na(FST))
df_spring_maf0.05_filt_500k <- df_spring_maf0.05_filtered[sample(nrow(df_spring_maf0.05_filtered), 5e5),]
ggplot(df_spring_maf0.05_filt_500k) + geom_histogram(aes(x=FST), bins=100)

df_autumn_maf0.05_filtered <- df_autumn_maf0.05 %>% filter(FST>=0 & !is.na(FST))
df_autumn_maf0.05_filt_500k <- df_autumn_maf0.05_filtered[sample(nrow(df_autumn_maf0.05_filtered), 5e5),]
ggplot(df_autumn_maf0.05_filt_500k) + geom_histogram(aes(x=FST), bins=100)

# 50 K SNPS
df_spring_maf0.05_filt_50K <- df_spring_maf0.05_filtered[sample(nrow(df_spring_maf0.05_filtered), 5e4),]
ggplot(df_spring_maf0.05_filt_50K) + geom_histogram(aes(x=FST), bins=20, color="black", fill="white")

df_autumn_maf0.05_filt_50K <- df_autumn_maf0.05_filtered[sample(nrow(df_autumn_maf0.05_filtered), 5e4),]
ggplot(df_autumn_maf0.05_filt_50K) + geom_histogram(aes(x=FST), bins=20, color="black", fill="white")


save(df_spring_maf0.05, file="fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData")

library(GenomicRanges)
library(dplyr)


# 2. LOAD INTROGRESSION REGIONS ####
intro_regions <- read.table("../4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt", header=T)
intro_reg_gr <- makeGRangesFromDataFrame(intro_regions)
intro_regions$name <- paste0(intro_regions$seqnames, "_", intro_regions$start, "_", intro_regions$end)

# 3. OVERLAP SPRING FST ####
# Convert Spring SNPs to Genomic Ranges
df_spring_gr <- GRanges(seqnames = df_spring$Chr,
                        IRanges(start = df_spring$Pos, end = df_spring$Pos))

overlaps_spring <- findOverlaps(df_spring_gr, intro_reg_gr)

# Extract SNPs that fall within regions
overlaps_spring_df <- df_spring[overlaps_spring@from, ]
overlaps_spring_df$intro_reg <- intro_regions[overlaps_spring@to, ]$name

# Summarize top FST SNP per region for Spring
results_fst_spring <- overlaps_spring_df %>%
  group_by(intro_reg) %>%
  summarise(
    Position_Spring = Pos[which.max(FST)],
    Max_FST_Spring = max(FST, na.rm = TRUE),
    Mean_FST_Spring = mean(FST, na.rm = TRUE),
    n_SNPs = n()
  )

# 4. OVERLAP AUTUMN FST ####
# Convert Autumn SNPs to Genomic Ranges
df_autumn_gr <- GRanges(seqnames = df_autumn$Chr,
                        IRanges(start = df_autumn$Pos, end = df_autumn$Pos))

overlaps_autumn <- findOverlaps(df_autumn_gr, intro_reg_gr)

# Extract SNPs that fall within regions
overlaps_autumn_df <- df_autumn[overlaps_autumn@from, ]
overlaps_autumn_df$intro_reg <- intro_regions[overlaps_autumn@to, ]$name

# Summarize top FST SNP per region for Autumn
results_fst_autumn <- overlaps_autumn_df %>%
  group_by(intro_reg) %>%
  summarise(
    Position_Autumn = Pos[which.max(FST)],
    Max_FST_Autumn = max(FST, na.rm = TRUE),
    Mean_FST_Autumn = mean(FST, na.rm = TRUE),
    n_SNPs = n()
  )

# 5. FINAL TABLE GENERATION ####
# Combine Spring and Autumn results into one table for comparison
final_comparison_table <- full_join(results_fst_spring, results_fst_autumn, by = "intro_reg")

# Format columns for easy viewing (Position / Fst)
final_comparison_table$Spring_Summary <- paste0(
  formatC(final_comparison_table$Position_Spring, big.mark=","), " / Fst: ", 
  sprintf("%.3f", final_comparison_table$Max_FST_Spring)
)

final_comparison_table$Autumn_Summary <- paste0(
  formatC(final_comparison_table$Position_Autumn, big.mark=","), " / Fst: ", 
  sprintf("%.3f", final_comparison_table$Max_FST_Autumn)
)

# 6. OUTPUT ####
View(final_comparison_table[, c("intro_reg", "Spring_Summary", "Autumn_Summary")])

write.table(final_comparison_table, "results/FST_overlap_introgression_table.txt", 
            col.names = T, row.names = F, quote = F, sep = "\t")

# Write Fst outputs
save(df_spring, file = "fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst")
save(df_autumn, file = "fst_results/baltic_autumn_vs_atlantic_autumn_han_pops_Fst")



#### Calculate average Fst in introgression regions ####
load("fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05.RData")
load("fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst")

df_spring_maf0.05 <- df_spring_maf0.05 %>% filter(FST >=0 & !is.na(FST))
# 3. OVERLAP SPRING FST ##### 3. OVERLAP SPRING FST ####FST
# Convert Spring SNPs to Genomic Ranges
df_spring_gr <- GRanges(seqnames = df_spring_maf0.05$Chr,
                        IRanges(start = df_spring_maf0.05$Pos, end = df_spring_maf0.05$Pos))

overlaps_spring <- findOverlaps(df_spring_gr, intro_reg_gr)

# Extract SNPs that fall within regions
overlaps_spring_df <- df_spring_maf0.05[overlaps_spring@from, ]
overlaps_spring_df$intro_reg <- intro_regions[overlaps_spring@to, ]$name

mean(overlaps_spring_df$FST, na.rm=T)
0.2125821
median(overlaps_spring_df$FST, na.rm=T)
0.1676966

mean(df_spring_maf0.05$FST, na.rm=T)
0.0391534
median(df_spring_maf0.05$FST, na.rm=T)
0.02199862

save(overlaps_spring_df, file="fst_results/baltic_spring_vs_atlantic_spring_han_pops_Fst_maf0.05_introgression_regions.RData")

# table <- read.table("../7.pixy/results_WGS_dxy/2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_fst.txt", header = T)
# #table <- read.table("../../../../Inversion_Project/Dxy/pixy_miss20maf0.01_28032023/datafiles/all_herring_populations.wg.20kb.popgenpixy.out_fst.txt", header=T)
# table %>% filter(pop1=="Baltic_Spring" & pop2=="Canada_Spring") %>% filter(no_snps > 100) %>% summarise(mean(avg_wc_fst))
# table %>% filter(pop1=="Baltic_Spring" ) %>% filter(no_snps > 100 & avg_wc_fst >= 0) %>% group_by(pop2) %>% summarise(mean(avg_wc_fst))
# 
# table %>% filter(pop1=="Canada_Spring" )  %>% group_by(pop2) %>% summarise(median(avg_wc_fst, na.rm=T))
# 
# table_fst <- table
# chr6<-table_fst[table_fst$chromosome=="chr6" & table_fst$window_pos_1>=22282765 & table_fst$window_pos_2 <=24868581,]
# chr12<-table_fst[table_fst$chromosome=="chr12" & table_fst$window_pos_1>=17826318 & table_fst$window_pos_2 <=25603093,]
# chr17<-table_fst[table_fst$chromosome=="chr17" & table_fst$window_pos_1>=25805445 & table_fst$window_pos_2 <=27568511,]
# chr23<-table_fst[table_fst$chromosome=="chr23" & table_fst$window_pos_1>=16226443 & table_fst$window_pos_2 <=17604273,]
# 
# inversions<-rbind(chr6, chr12, chr17, chr23)
# to_exclude<-rownames(inversions)
# all_rows<-rownames(table_fst)
# all_rows %in% to_exclude
# 
# table_fst_without_inversions<-table_fst[!(all_rows %in% to_exclude),]
# 
# table_fst_without_inversions %>% 
#   group_by(pop1, pop2) %>%
#   filter(avg_wc_fst >= 0) %>%
#   summarise(median(avg_wc_fst, na.rm=T))
# 
