# LIBRARIES ####
library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(data.table)

# INPUTS ####
# Set the working directory
# Let's load the snpEFF annotations of SNPs based on the high coverage data:
setwd("Manuscript/Figshare/")

# Read in introgressed regions coordinates. Let's use non-collapsed for now
intro_regions<-read.table(header=T, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

# Annotations #
scan1_match_df<-read.table("4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_summary_filter2_cov7_gr_min50kb.maxgap20K.modified.txt", sep="\t", header=T, row.names=NULL)
scan1_match_df_gr <-GRanges(seqnames=scan1_match_df$seqnames, 
                            IRanges(start=scan1_match_df$start, end=scan1_match_df$end),
                            ensemble_gene_id=scan1_match_df$ensembl_gene_id, 
                            external_gene_name=scan1_match_df$external_gene_name,
                            description=scan1_match_df$description_type)

# ChiSquare
load("5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_arctic_pacific_han_pops_chiseq.output.RData")
load("5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData")
load("5.processing_pool_seq_data/chisquare_results/baltic_autumn_vs_atlantic_autumn_han_pops_chiseq.output.RData")

# Calculate MAF for each contrast
baltic_autumn_vs_atlantic_autumn_chiseq$maf<- with(baltic_autumn_vs_atlantic_autumn_chiseq, pmin(pop1_A + pop2_A, pop1_D + pop2_D) / 
                                                               (pop1_A + pop1_D + pop2_A + pop2_D))

baltic_spring_vs_arctic_pacific_chiseq$maf<- with(baltic_spring_vs_arctic_pacific_chiseq, pmin(pop1_A + pop2_A, pop1_D + pop2_D) / 
                                                     (pop1_A + pop1_D + pop2_A + pop2_D))


baltic_spring_vs_atlantic_spring_chiseq$maf<- with(baltic_spring_vs_atlantic_spring_chiseq, pmin(pop1_A + pop2_A, pop1_D + pop2_D) / 
                                                     (pop1_A + pop1_D + pop2_A + pop2_D))


# baltic_spring_vs_arctic_pacific_chiseq_maf005 <- baltic_spring_vs_arctic_pacific_chiseq %>% filter(maf > 0.05)
# baltic_spring_vs_atlantic_spring_chiseq_maf005 <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(maf > 0.05)
# baltic_autumn_vs_atlantic_autumn_chiseq_maf005 <- baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(maf > 0.05)


baltic_spring_vs_arctic_pacific_chiseq_maf001 <- baltic_spring_vs_arctic_pacific_chiseq %>% filter(maf > 0.01)
baltic_spring_vs_atlantic_spring_chiseq_maf001 <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(maf > 0.01)
baltic_autumn_vs_atlantic_autumn_chiseq_maf001 <- baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(maf > 0.01)

# ANALYSIS ####
# Find mutations in introgressed regions ####
# Convert introgression regions to Granges 
intro_regions_df<-makeGRangesFromDataFrame(intro_regions)

# SNPEff mutations
annotations_list<-list.files(path=".", pattern="SNP_annotations")

annotations<-lapply(annotations_list, FUN=read.table)

annotations_df<-rbindlist(annotations, use.names=TRUE)
colnames(annotations_df)<-c("CHR", "POS", "REF", "ALT", "GENE", "TRANSCRIPT", "TYPE", "TYPE_2")

annotations_df_gr<-GRanges(seqnames=annotations_df$CHR, IRanges(start=annotations_df$POS, end=annotations_df$POS), 
                           REF=annotations_df$REF,
                           ALT=annotations_df$ALT,
                           GENE=annotations_df$GENE,
                           TRANSCRIPT=annotations_df$TRANSCRIPT,
                           TYPE=annotations_df$TYPE,
                           TYPE_2=annotations_df$TYPE_2)

# Mutations of interest are the ones in and around the introgression regions:
mutations_of_interest <-subsetByOverlaps(annotations_df_gr,intro_regions_df, maxgap = 2e4)
#mutations_of_interest$GENE == scan1_match_df$ensembl_gene_id

# Thank you chatGPT. I created a dictionary like thing so I can add gene names to
# the Granges above (so freaking easy......)
annotation_dict <- setNames(scan1_match_df[[7]], scan1_match_df[[6]])
mutations_of_interest$gene_name <- annotation_dict[c(mutations_of_interest$GENE)]

# So, now convert the chisquare tables into Granges so that I can find the overlap:

## 2. Function to turn a contrast data.frame into GRanges
contrast_to_gr <- function(df) {
  
  if (!all(c("CHROM", "POS") %in% names(df))) {
    stop("df must contain CHROM and POS columns")
  }
  df <- df[!is.na(df$CHROM) & !is.na(df$POS), , drop = FALSE]
  
  # Ensure columns exist
  if (!"dAF"       %in% names(df)) df$dAF       <- NA_real_
  if (!"chisq_p"   %in% names(df)) df$chisq_p   <- NA_real_
  if (!"freqpop1"  %in% names(df)) df$freqpop1  <- NA_real_
  if (!"freqpop2"  %in% names(df)) df$freqpop2  <- NA_real_
  
  GRanges(
    seqnames  = df$CHROM,
    ranges    = IRanges(start = as.integer(df$POS), end = as.integer(df$POS)),
    dAF       = as.numeric(df$dAF),
    pval      = as.numeric(df$chisq_p),
    freqpop1  = as.numeric(df$freqpop1),
    freqpop2  = as.numeric(df$freqpop2)
  )
}

## 3. Make GRanges for each contrast
# gr_contrast_BS_AS <- contrast_to_gr(baltic_spring_vs_atlantic_spring_chiseq)
# gr_contrast_BA_AA <- contrast_to_gr(baltic_autumn_vs_atlantic_autumn_chiseq)
# gr_contrast_BS_AP <- contrast_to_gr(baltic_spring_vs_arctic_pacific_chiseq)

# gr_contrast_BS_AS_maf5 <- contrast_to_gr(baltic_spring_vs_atlantic_spring_chiseq_maf005)
# gr_contrast_BA_AA_maf5 <- contrast_to_gr(baltic_autumn_vs_atlantic_autumn_chiseq_maf005)
# gr_contrast_BS_AP_maf5 <- contrast_to_gr(baltic_spring_vs_arctic_pacific_chiseq_maf005)

gr_contrast_BS_AS_maf1 <- contrast_to_gr(baltic_spring_vs_atlantic_spring_chiseq_maf001)
gr_contrast_BA_AA_maf1 <- contrast_to_gr(baltic_autumn_vs_atlantic_autumn_chiseq_maf001)
gr_contrast_BS_AP_maf1 <- contrast_to_gr(baltic_spring_vs_arctic_pacific_chiseq_maf001)

## 4. Function to match and add columns
add_contrast_info <- function(mut_gr, contrast_gr, stat_prefix, pop1_name, pop2_name) {
  
  #mut_gr <- mutations_of_interest
  #contrast_gr <- gr_contrast_BS_AS
  #prefix <- "BS_AS"
  
  # Column names for stats
  nm_daf  <- paste0(stat_prefix, "_dAF")
  nm_pval <- paste0(stat_prefix, "_pval")
  # Column names for frequencies
  nm_freq1 <- paste0(pop1_name, "_freq")
  nm_freq2 <- paste0(pop2_name, "_freq")
  
  # Preallocate columns
  mcols(mut_gr)[[nm_daf]]   <- rep(NA_real_, length(mut_gr))
  mcols(mut_gr)[[nm_pval]]  <- rep(NA_real_, length(mut_gr))
  mcols(mut_gr)[[nm_freq1]] <- rep(NA_real_, length(mut_gr))
  mcols(mut_gr)[[nm_freq2]] <- rep(NA_real_, length(mut_gr))
  
  # Match GRanges
  hits <- findOverlaps(mut_gr, contrast_gr)
  if (length(hits) == 0) return(mut_gr)
  
  q <- queryHits(hits)
  s <- subjectHits(hits)
  
  # If multiple subject hits map to the same query, keep the first one by default.
  # I don't think this will happen! There won't be more than one position in the contrast 
  # data frame
  
  keep <- !duplicated(q)
  q_unique <- q[keep]
  s_first  <- s[keep]
  
  mcols(mut_gr)[[nm_daf]][q_unique]  <- mcols(contrast_gr)$dAF[s_first]
  mcols(mut_gr)[[nm_pval]][q_unique] <- mcols(contrast_gr)$pval[s_first]
  mcols(mut_gr)[[nm_freq1]][q_unique] <- mcols(contrast_gr)$freqpop1[s_first]
  mcols(mut_gr)[[nm_freq2]][q_unique] <- mcols(contrast_gr)$freqpop2[s_first]
  
  mut_gr
  
}


## 5. Add data from all three contrasts
# Let's make a copy so we retain the initial output:
# mut_gr <- mutations_of_interest

# mut_gr <- add_contrast_info(mut_gr, gr_contrast_BS_AS, "BS_AS", "BS", "AS")
# mut_gr <- add_contrast_info(mut_gr, gr_contrast_BA_AA, "BA_AA", "BA", "AA")
# mut_gr <- add_contrast_info(mut_gr, gr_contrast_BS_AP, "BS_AP", "BS", "AP")

# mut_gr_df <- data.frame(mut_gr)

# mut_gr_maf5 <- mutations_of_interest
# mut_gr_maf5 <- add_contrast_info(mut_gr_maf5, gr_contrast_BS_AS_maf5, "BS_AS", "BS", "AS")
# mut_gr_maf5 <- add_contrast_info(mut_gr_maf5, gr_contrast_BA_AA_maf5, "BA_AA", "BA", "AA")
# mut_gr_maf5 <- add_contrast_info(mut_gr_maf5, gr_contrast_BS_AP_maf5, "BS_AP", "BS", "AP")

# mut_gr_maf5_df <- data.frame(mut_gr_maf5)

mut_gr_maf1 <- mutations_of_interest

mut_gr_maf1 <- add_contrast_info(mut_gr_maf1, gr_contrast_BS_AS_maf1, "BS_AS", "BS", "AS")
mut_gr_maf1 <- add_contrast_info(mut_gr_maf1, gr_contrast_BA_AA_maf1, "BA_AA", "BA", "AA")
mut_gr_maf1 <- add_contrast_info(mut_gr_maf1, gr_contrast_BS_AP_maf1, "BS_AP", "BS", "AP")

mut_gr_maf1_df <- data.frame(mut_gr_maf1)


# Let's double check our code:
mut_gr_df %>% filter(seqnames=="chr19" & start==6359267)

baltic_spring_vs_arctic_pacific_chiseq %>% filter(CHROM=="chr19" & POS==6359267)
baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM=="chr19" & POS==6362829)
baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(CHROM=="chr19" & POS==6362829)

# mut_gr_maf5_df %>% filter(seqnames=="chr19" & start==6359267)
mut_gr_maf1_df %>% filter(seqnames=="chr19" & start==6359267)

# Create a new data.frame and 

# Define priority order for TYPE_2 when TYPE is not "HIGH"
type2_priority <- c(
  "missense_variant",
  "synonymous_variant",
  "intron_variant",
  "5_prime_UTR_variant",
  "3_prime_UTR_variant",
  "downstream_gene_variant",
  "upstream_gene_variant",
  "intergenic_region"
)

# Assign a numeric rank based on TYPE / TYPE_2
# mut_gr_df_2 <- mut_gr_df %>%
#   mutate(
#     priority = case_when(
#       TYPE == "HIGH" ~ 0L,  # all equally top priority if HIGH
#       TYPE_2 %in% type2_priority ~ match(TYPE_2, type2_priority),
#       TRUE ~ length(type2_priority) + 1L  # lowest priority if not in list
#     )
#   )


# mut_gr_maf5_df_2 <- mut_gr_maf5_df %>%
#   mutate(
#     priority = case_when(
#       TYPE == "HIGH" ~ 0L,  # all equally top priority if HIGH
#       TYPE_2 %in% type2_priority ~ match(TYPE_2, type2_priority),
#       TRUE ~ length(type2_priority) + 1L  # lowest priority if not in list
#     )
#   )

mut_gr_maf1_df_2 <- mut_gr_maf1_df %>%
  mutate(
    priority = case_when(
      TYPE == "HIGH" ~ 0L,  # all equally top priority if HIGH
      TYPE_2 %in% type2_priority ~ match(TYPE_2, type2_priority),
      TRUE ~ length(type2_priority) + 1L  # lowest priority if not in list
    )
  )

# Group by genomic coordinates (and allele if needed), keep best priority
# dedup_df <- mut_gr_df_2 %>%
#   group_by(seqnames, start, end, REF, ALT) %>%
#   slice_min(order_by = priority, with_ties = FALSE) %>%
#   ungroup() %>%
#   dplyr::select(-priority)  # remove helper column

# dedup_df %>% filter(seqnames=="chr19" & start==6359267)


# dedup_df_maf5 <- mut_gr_maf5_df_2 %>%
#   group_by(seqnames, start, end, REF, ALT) %>%
#   slice_min(order_by = priority, with_ties = FALSE) %>%
#   ungroup() %>%
#   dplyr::select(-priority)  # remove helper column

dedup_df_maf1 <- mut_gr_maf1_df_2 %>%
  group_by(seqnames, start, end, REF, ALT) %>%
  slice_min(order_by = priority, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-priority)  # remove helper column

# Ask questions about the data:



# # So, I chatted with Leif and we agreed to use a "normal" Pvalue. I will still use 0.01 
# # dedup_df has 142216 unique SNPs. So we can apply a bonferroni correction of 6.5

# # 91 SNPs
# MODERATE_SNPs <- dedup_df %>% 
#   filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
#   filter(TYPE=="MODERATE") %>% 
#   dplyr::select(seqnames, start, REF, ALT, GENE, gene_name, TYPE, TYPE_2,  
#                 BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# # 4 SNPs
# HIGH_SNPS <- dedup_df %>% 
#   filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
#   filter(TYPE=="HIGH") %>% 
#   dplyr::select(seqnames, start, REF, ALT, GENE, gene_name, TYPE, TYPE_2,  
#                 BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# # 99 SNPs
# LOW_SNPs <- dedup_df %>% 
#   filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
#   filter(TYPE=="LOW") %>% 
#   dplyr::select(seqnames, start, REF, ALT, GENE, gene_name, TYPE, TYPE_2,  
#                 BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# # 3673 SNPs
# MODIFIER_SNPs <- dedup_df %>% 
#   filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
#   filter(TYPE=="MODIFIER") %>% 
#   dplyr::select(seqnames, start, REF, ALT, GENE, gene_name, TYPE, TYPE_2,  
#                 BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# colnames(MODERATE_SNPs)<-c("Chromosome", "Position", "Reference", "Alternative", "Gene_ID", "Gene_name",
#                            "Type", "Function", "BS_AS_dAF", "BS_AS_pval", "BS_AP_dAF", "BS_AP_pval", "BA_AA_dAF", "BA_AA_pval")

# colnames(HIGH_SNPS)<-c("Chromosome", "Position", "Reference", "Alternative", "Gene_ID", "Gene_name",
#                            "Type", "Function", "BS_AS_dAF", "BS_AS_pval", "BS_AP_dAF", "BS_AP_pval", "BA_AA_dAF", "BA_AA_pval")

# colnames(LOW_SNPs)<-c("Chromosome", "Position", "Reference", "Alternative", "Gene_ID", "Gene_name",
#                        "Type", "Function", "BS_AS_dAF", "BS_AS_pval", "BS_AP_dAF", "BS_AP_pval", "BA_AA_dAF", "BA_AA_pval")

# colnames(MODIFIER_SNPs)<-c("Chromosome", "Position", "Reference", "Alternative", "Gene_ID", "Gene_name",
#                       "Type", "Function", "BS_AS_dAF", "BS_AS_pval", "BS_AP_dAF", "BS_AP_pval", "BA_AA_dAF", "BA_AA_pval")


# write.table(MODERATE_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "resuls_candidate_mutations_2025-08-31/moderate_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
# write.table(HIGH_SNPS, quote=F, col.names = T, row.names=F, sep="\t", file = "resuls_candidate_mutations_2025-08-31/high_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
# write.table(LOW_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "resuls_candidate_mutations_2025-08-31/low_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
# write.table(MODIFIER_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "resuls_candidate_mutations_2025-08-31/modifier_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")


# Supplementary Table 6 ####
# Apply a maf filter of 0.01
# Filter also biallelic sites:
dedup_df_maf1_bi<-dedup_df_maf1 %>% filter(nchar(ALT) == 1)

# -log10(0.05/137098)
# 6.438061

# 84 SNPs
MODERATE_SNPs <- dedup_df_maf1_bi %>% 
  filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
  filter(TYPE=="MODERATE") %>% 
  dplyr::select(seqnames, start, REF, ALT, GENE,TYPE, TYPE_2,   gene_name, 
                BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval) 

# 3 SNPs
HIGH_SNPs <- dedup_df_maf1_bi %>% 
  filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
  filter(TYPE=="HIGH") %>% 
  dplyr::select(seqnames, start, REF, ALT, GENE, TYPE, TYPE_2,   gene_name, 
                BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# 96 SNPs
LOW_SNPs <- dedup_df_maf1_bi %>% 
  filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
  filter(TYPE=="LOW") %>% 
  dplyr::select(seqnames, start, REF, ALT, GENE, TYPE, TYPE_2,  gene_name, 
                BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval)

# 3352 SNPs
MODIFIER_SNPs <- dedup_df_maf1_bi  %>% 
  filter(BS_AS_pval < 1e-7 & BS_AP_pval > 1e-7) %>%
  filter(TYPE=="MODIFIER") %>% 
  dplyr::select(seqnames, start, REF, ALT, GENE, TYPE, TYPE_2, gene_name,  
                BS_AS_dAF, BS_AS_pval, BS_AP_dAF, BS_AP_pval, BA_AA_dAF, BA_AA_pval) 

write.table(MODERATE_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "results/moderate_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
write.table(HIGH_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "results/high_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
write.table(LOW_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "results/low_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")
write.table(MODIFIER_SNPs, quote=F, col.names = T, row.names=F, sep="\t", file = "results/modifier_mutations_20kb_BS_AS_sig_BS_AP_nonsig.txt")

