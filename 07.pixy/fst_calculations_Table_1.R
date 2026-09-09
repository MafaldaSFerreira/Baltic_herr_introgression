# Read packages
library(tidyverse)
library(GenomicRanges)

# Session Info
sessionInfo()

# Set working directory
# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
setwd("../Figshare/7.pixy/results_WGS_dxy/")

# Read in Fst results 
fst_table <- read.table("2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_fst.txt", header=T)

# We should compare Atlantic spring vs Baltic spring; Atlantic Autumn vs Baltic Autumn
# Based on Table S1, these will Baltic Spring, Canada Spring and Norway Pelagic individuals
# Based on other analysis, Norway Costal may not follow the genetic patterns of Norway Pelagic individuals,
# so I will not include them here.
fst_table_spring <- fst_table %>% 
  filter(pop1 %in% c("Baltic_Spring")) %>% 
  filter(pop2 %in% c("Canada_Spring", "Norway_Pelagic")) %>%
  filter(no_snps >200 & avg_wc_fst >0) # let's try this

# For Autumn, I will include also Canada Autumn and Baltic Autumn. 
# I could include from UK, but for now let's see these as I need not separate them all
# in the analysis
fst_table_autumn <- fst_table %>% 
  filter(pop1 %in% c("Baltic_Autumn")) %>% 
  filter( pop2 %in% c("Canada_Autumn")) %>%
  filter(no_snps >200 & avg_wc_fst >0)

# convert to GRanges
fst_table_spring_gr <- GRanges(seqnames=fst_table_spring$chromosome, 
                               IRanges(start=fst_table_spring$window_pos_1,
                                       end=fst_table_spring$window_pos_2), 
                              fst=fst_table_spring$avg_wc_fst)

fst_table_autumn_gr <- GRanges(seqnames=fst_table_autumn$chromosome, 
                               IRanges(start=fst_table_autumn$window_pos_1,
                                       end=fst_table_autumn$window_pos_2), 
                               fst=fst_table_autumn$avg_wc_fst)

# Table 1 ####
# Find overlaps with the introgression regions:
# Read in introgressed regions coordinates. Let's use non-collapsed for now
intro_regions<-read.table(header=T, "../../4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")
intro_reg_gr <- makeGRangesFromDataFrame(intro_regions)
intro_regions$name <- paste0(intro_regions$seqnames, "_", intro_regions$start, "_", intro_regions$end)

# Spring
overlaps_spring_intro_reg <- findOverlaps(fst_table_spring_gr, intro_reg_gr)

overlaps_spring_intro_reg_df <- fst_table_spring[overlaps_spring_intro_reg@from,]
overlaps_spring_intro_reg_df$intro_reg <- NA
overlaps_spring_intro_reg_df$intro_reg <- intro_regions[overlaps_spring_intro_reg@to,]$name



# --------------------- #
# R version 4.5.1 (2025-06-13)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Europe/Stockholm
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lubridate_1.9.4        forcats_1.0.1          stringr_1.5.2          dplyr_1.1.4           
# [5] purrr_1.1.0            readr_2.1.5            tidyr_1.3.1            tibble_3.3.0          
# [9] ggplot2_4.0.0          tidyverse_2.0.0        GenomicFeatures_1.60.0 AnnotationDbi_1.70.0  
# [13] Biobase_2.68.0         GenomicRanges_1.60.0   GenomeInfoDb_1.44.3    IRanges_2.42.0        
# [17] S4Vectors_0.46.0       BiocGenerics_0.54.1    generics_0.1.4        
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1            farver_2.1.2                blob_1.2.4                 
# [4] filelock_1.0.3              Biostrings_2.76.0           S7_0.2.0                   
# [7] bitops_1.0-9                fastmap_1.2.0               RCurl_1.98-1.17            
# [10] BiocFileCache_2.16.2        GenomicAlignments_1.44.0    XML_3.99-0.19              
# [13] digest_0.6.37               timechange_0.3.0            lifecycle_1.0.4            
# [16] KEGGREST_1.48.1             RSQLite_2.4.3               magrittr_2.0.4             
# [19] compiler_4.5.1              rlang_1.1.6                 progress_1.2.3             
# [22] tools_4.5.1                 yaml_2.3.10                 rtracklayer_1.68.0         
# [25] prettyunits_1.2.0           S4Arrays_1.8.1              bit_4.6.0                  
# [28] curl_7.0.0                  DelayedArray_0.34.1         xml2_1.4.0                 
# [31] RColorBrewer_1.1-3          abind_1.4-8                 BiocParallel_1.42.2        
# [34] txdbmaker_1.4.2             withr_3.0.2                 grid_4.5.1                 
# [37] scales_1.4.0                biomaRt_2.64.0              SummarizedExperiment_1.38.1
# [40] cli_3.6.5                   crayon_1.5.3                rstudioapi_0.17.1          
# [43] httr_1.4.7                  tzdb_0.5.0                  rjson_0.2.23               
# [46] DBI_1.2.3                   cachem_1.1.0                parallel_4.5.1             
# [49] XVector_0.48.0              restfulr_0.0.16             matrixStats_1.5.0          
# [52] vctrs_0.6.5                 Matrix_1.7-4                jsonlite_2.0.0             
# [55] hms_1.1.4                   bit64_4.6.0-1               glue_1.8.0                 
# [58] codetools_0.2-20            stringi_1.8.7               gtable_0.3.6               
# [61] BiocIO_1.18.0               UCSC.utils_1.4.0            pillar_1.11.1              
# [64] rappdirs_0.3.3              GenomeInfoDbData_1.2.14     R6_2.6.1                   
# [67] dbplyr_2.5.1                httr2_1.2.1                 lattice_0.22-7             
# [70] png_0.1-8                   Rsamtools_2.24.1            memoise_2.0.1              
# [73] SparseArray_1.8.1           MatrixGenerics_1.20.0       pkgconfig_2.0.3    
