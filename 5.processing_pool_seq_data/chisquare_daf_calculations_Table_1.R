# This R code will take the Neff Corrected Allele Counts file from Han et al 2020 and perform Chi-Square tests
# to identify SNPs with significant allele frequency differences between several population constrasts.

# Input files for this can be found in FigShare


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
AD<-read.table("60.Neff.AD", header=T)
names <- names(AD)

# Calculate allele frequencies
# I generated this function to perform the chisq calculations between populations

pop_allele_chisq <- function(count_df, pop1, pop2){
  out_df <- count_df[,c("CHROM", "POS")]
  out_df[,"pop1_A"] <- 0
  out_df[,"pop1_D"] <- 0
  for(p1_col in pop1){
    out_df[,"pop1_A"] <- out_df[,"pop1_A"] + as.integer(sub("([0-9]+),([0-9]+)", "\\1", count_df[,p1_col]))
    out_df[,"pop1_D"] <- out_df[,"pop1_D"] + as.integer(sub("([0-9]+),([0-9]+)", "\\2", count_df[,p1_col]))
  }
  
  out_df[,"pop2_A"] <- 0
  out_df[,"pop2_D"] <- 0
  for(p2_col in pop2){
    out_df[,"pop2_A"] <- out_df[,"pop2_A"] + as.integer(sub("([0-9]+),([0-9]+)", "\\1", count_df[,p2_col]))
    out_df[,"pop2_D"] <- out_df[,"pop2_D"] + as.integer(sub("([0-9]+),([0-9]+)", "\\2", count_df[,p2_col]))
  }
  out_df[,"total"] <- out_df[,"pop1_A"] + out_df[,"pop1_D"] + out_df[,"pop2_A"] + out_df[,"pop2_D"]
  out_df[,"pop1_A_exp"] <- ((out_df[,"pop1_A"] + out_df[,"pop2_A"])/out_df[,"total"]) * (out_df[,"pop1_A"] + out_df[,"pop1_D"])
  out_df[,"pop1_D_exp"] <- ((out_df[,"pop1_D"] + out_df[,"pop2_D"])/out_df[,"total"]) * (out_df[,"pop1_A"] + out_df[,"pop1_D"])
  out_df[,"pop2_A_exp"] <- ((out_df[,"pop1_A"] + out_df[,"pop2_A"])/out_df[,"total"]) * (out_df[,"pop2_A"] + out_df[,"pop2_D"])
  out_df[,"pop2_D_exp"] <- ((out_df[,"pop1_D"] + out_df[,"pop2_D"])/out_df[,"total"]) * (out_df[,"pop2_A"] + out_df[,"pop2_D"])
  
  chisq_site_filter <- (out_df[,"pop1_A_exp"] > 0 | out_df[,"pop2_A_exp"] > 0) & 
    (out_df[,"pop1_D_exp"] > 0 | out_df[,"pop2_D_exp"] > 0) & 
    (out_df[,"pop1_A_exp"] + out_df[,"pop1_D_exp"]) > 10*length(pop1) & 
    (out_df[,"pop2_A_exp"] + out_df[,"pop2_D_exp"]) > 10*length(pop2)   
  
  
  out_df <-  out_df[chisq_site_filter,]
  
  #Continuity correction
  pop1_A_diff <- abs(out_df[,"pop1_A"] - out_df[,"pop1_A_exp"]) - 0.5
  pop1_A_diff[pop1_A_diff < 0] <- 0
  pop1_D_diff <- abs(out_df[,"pop1_D"] - out_df[,"pop1_D_exp"]) - 0.5
  pop1_D_diff[pop1_D_diff < 0] <- 0
  pop2_A_diff <- abs(out_df[,"pop2_A"] - out_df[,"pop2_A_exp"]) - 0.5
  pop2_A_diff[pop2_A_diff < 0] <- 0
  pop2_D_diff <- abs(out_df[,"pop2_D"] - out_df[,"pop2_D_exp"]) - 0.5
  pop2_D_diff[pop2_D_diff < 0] <- 0
  
  out_df[,"chisq_stat"] <-  ((pop1_A_diff^2)/out_df[,"pop1_A_exp"] + 
                               (pop1_D_diff^2)/out_df[,"pop1_D_exp"] +
                               (pop2_A_diff^2)/out_df[,"pop2_A_exp"] + 
                               (pop2_D_diff^2)/out_df[,"pop2_D_exp"]) 
  #Uncorrectred version
  #out_df[,"chisq_stat"] <-  (((out_df[,"pop1_A"] - out_df[,"pop1_A_exp"])^2)/out_df[,"pop1_A_exp"] + 
  #                             ((out_df[,"pop1_D"] - out_df[,"pop1_D_exp"])^2)/out_df[,"pop1_D_exp"] +
  #                             ((out_df[,"pop2_A"] - out_df[,"pop2_A_exp"])^2)/out_df[,"pop2_A_exp"] + 
  #                             ((out_df[,"pop2_D"] - out_df[,"pop2_D_exp"])^2)/out_df[,"pop2_D_exp"]) 
  
  out_df[,"chisq_p"] <- pchisq(q = out_df[,"chisq_stat"], df = 1, lower.tail = F)
  return(out_df)
}



## Baltic spring vs Atlantic spring ####
baltic_spring_vs_atlantic_spring_chiseq<-pop_allele_chisq(AD, Baltic_Spring_col, Atlantic_Spring_col)
baltic_spring_vs_atlantic_spring_chiseq$freqpop1 <- baltic_spring_vs_atlantic_spring_chiseq$pop1_A / (baltic_spring_vs_atlantic_spring_chiseq$pop1_A + baltic_spring_vs_atlantic_spring_chiseq$pop1_D)
baltic_spring_vs_atlantic_spring_chiseq$freqpop2 <- baltic_spring_vs_atlantic_spring_chiseq$pop2_A / (baltic_spring_vs_atlantic_spring_chiseq$pop2_A + baltic_spring_vs_atlantic_spring_chiseq$pop2_D)
baltic_spring_vs_atlantic_spring_chiseq$dAF <- abs(baltic_spring_vs_atlantic_spring_chiseq$freqpop1 - baltic_spring_vs_atlantic_spring_chiseq$freqpop2)
baltic_spring_vs_atlantic_spring_chiseq$SNP<-paste0(baltic_spring_vs_atlantic_spring_chiseq$CHROM,"_", baltic_spring_vs_atlantic_spring_chiseq$POS)
save(baltic_spring_vs_atlantic_spring_chiseq, file = "baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData")
load("baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData")

## Baltic Autumn vs Atlantic Autumn ####
baltic_autumn_vs_atlantic_autumn_chiseq<-pop_allele_chisq(AD, Baltic_Autumn_col, Atlantic_Autumn_col)
baltic_autumn_vs_atlantic_autumn_chiseq$freqpop1 <- baltic_autumn_vs_atlantic_autumn_chiseq$pop1_A / (baltic_autumn_vs_atlantic_autumn_chiseq$pop1_A + baltic_autumn_vs_atlantic_autumn_chiseq$pop1_D)
baltic_autumn_vs_atlantic_autumn_chiseq$freqpop2 <- baltic_autumn_vs_atlantic_autumn_chiseq$pop2_A / (baltic_autumn_vs_atlantic_autumn_chiseq$pop2_A + baltic_autumn_vs_atlantic_autumn_chiseq$pop2_D)
baltic_autumn_vs_atlantic_autumn_chiseq$dAF <- abs(baltic_autumn_vs_atlantic_autumn_chiseq$freqpop1 - baltic_autumn_vs_atlantic_autumn_chiseq$freqpop2)
baltic_autumn_vs_atlantic_autumn_chiseq$SNP<-paste0(baltic_autumn_vs_atlantic_autumn_chiseq$CHROM,"_", baltic_autumn_vs_atlantic_autumn_chiseq$POS)
save(baltic_autumn_vs_atlantic_autumn_chiseq, file = "baltic_autumn_vs_atlantic_autumn_han_pops_chiseq.output.RData")
load("baltic_autumn_vs_atlantic_autumn_han_pops_chiseq.output.RData")

## Baltic Spring vs White + Pechora Sea ####
WhitePechora_sea_col <- names[grep("HWS[2345]", names)]
baltic_spring_vs_arctic_pacific_chiseq<-pop_allele_chisq(AD, Baltic_Spring_col, WhitePechora_sea_col)
baltic_spring_vs_arctic_pacific_chiseq$freqpop1 <- baltic_spring_vs_arctic_pacific_chiseq$pop1_A / (baltic_spring_vs_arctic_pacific_chiseq$pop1_A + baltic_spring_vs_arctic_pacific_chiseq$pop1_D)
baltic_spring_vs_arctic_pacific_chiseq$freqpop2 <- baltic_spring_vs_arctic_pacific_chiseq$pop2_A / (baltic_spring_vs_arctic_pacific_chiseq$pop2_A + baltic_spring_vs_arctic_pacific_chiseq$pop2_D)
baltic_spring_vs_arctic_pacific_chiseq$dAF <- abs(baltic_spring_vs_arctic_pacific_chiseq$freqpop1 - baltic_spring_vs_arctic_pacific_chiseq$freqpop2)
baltic_spring_vs_arctic_pacific_chiseq$SNP<-paste0(baltic_spring_vs_arctic_pacific_chiseq$CHROM,"_", baltic_spring_vs_arctic_pacific_chiseq$POS)
save(baltic_spring_vs_arctic_pacific_chiseq, file = "baltic_spring_vs_arctic_pacific_han_pops_chiseq.output.RData")
load("baltic_spring_vs_arctic_pacific_han_pops_chiseq.output.RData")


# Table 1 ####
# Find overlaps with the introgression regions:
# Read in introgressed regions coordinates. Let's use non-collapsed for now
intro_regions<-read.table(header=T, "~/Documents/Postdoc/Project_Herring/Introgression/introgression_scan/2023-11-16_analysis/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")
intro_reg_gr <- makeGRangesFromDataFrame(intro_regions)
intro_regions$name <- paste0(intro_regions$seqnames, "_", intro_regions$start, "_", intro_regions$end)

# Convert to GRanges:
baltic_spring_vs_atlantic_spring_chiseq_gr <- GRanges(seqnames=baltic_spring_vs_atlantic_spring_chiseq$CHROM,
                                                      IRanges(start=baltic_spring_vs_atlantic_spring_chiseq$POS,
                                                              end=baltic_spring_vs_atlantic_spring_chiseq$POS))


overlaps_bs_as_intro_reg <- findOverlaps(baltic_spring_vs_atlantic_spring_chiseq_gr, intro_reg_gr)

overlaps_bs_as_intro_reg_df <- baltic_spring_vs_atlantic_spring_chiseq[overlaps_bs_as_intro_reg@from,]
overlaps_bs_as_intro_reg_df$intro_reg <- NA
overlaps_bs_as_intro_reg_df$intro_reg <- intro_regions[overlaps_bs_as_intro_reg@to,]$name

data.frame(overlaps_bs_as_intro_reg_df) %>%
  group_by(intro_reg) %>%
  summarise(
    Position_Spring = POS[which.min(chisq_p)],
    dAF_Spring = dAF[which.min(chisq_p)],
    chi_Spring_top = min(chisq_p)) %>% 
  print(n=100) 


# For autumn
baltic_autumn_vs_atlantic_autumn_chiseq_gr <- GRanges(seqnames=baltic_autumn_vs_atlantic_autumn_chiseq$CHROM,
                                                      IRanges(start=baltic_autumn_vs_atlantic_autumn_chiseq$POS,
                                                              end=baltic_autumn_vs_atlantic_autumn_chiseq$POS))


overlaps_ba_aa_intro_reg <- findOverlaps(baltic_autumn_vs_atlantic_autumn_chiseq_gr, intro_reg_gr)

overlaps_ba_aa_intro_reg_df <- baltic_autumn_vs_atlantic_autumn_chiseq[overlaps_ba_aa_intro_reg@from,]
overlaps_ba_aa_intro_reg_df$intro_reg <- NA
overlaps_ba_aa_intro_reg_df$intro_reg <- intro_regions[overlaps_ba_aa_intro_reg@to,]$name


results_ba_aa <- data.frame(overlaps_ba_aa_intro_reg_df) %>%
  group_by(intro_reg) %>%
  summarise(
    Position_Autumn = POS[which.min(chisq_p)],
    dAF_Autumn = dAF[which.min(chisq_p)],
    AlleleFreq1 = freqpop1[which.min(chisq_p)],
    AlleleFreq2 = freqpop2[which.min(chisq_p)],
    chi_Autumn_top = min(chisq_p)) %>% 
  print(n=100) 

results_bs_as <- data.frame(overlaps_bs_as_intro_reg_df) %>%
  group_by(intro_reg) %>%
  summarise(
    Position_Spring = POS[which.min(chisq_p)],
    dAF_Spring = dAF[which.min(chisq_p)],
    AlleleFreq1 = freqpop1[which.min(chisq_p)],
    AlleleFreq2 = freqpop2[which.min(chisq_p)],
    chi_Spring_top = min(chisq_p)) %>% 
  print(n=100) 

results_bs_as$results <- paste0(formatC(results_bs_as$Position_Spring, big.mark=","), " / ", sprintf("%.1f", -log10(results_bs_as$chi_Spring_top))  ," / ", sprintf("%.2f", results_bs_as$dAF_Spring))

results_ba_aa$results <- paste0(formatC(results_ba_aa$Position_Autumn, big.mark=","), " / ", sprintf("%.1f", -log10(results_ba_aa$chi_Autumn_top))  ," / ", sprintf("%.2f", results_ba_aa$dAF_Autumn))

View(cbind(results_bs_as$results, results_ba_aa$results))

write.table(results_bs_as, "~/Documents/Postdoc/Repositories/Baltic_herr_introgression/5.processing_pool_seq_data/results/results_bs_as_table_1.txt", col.names = T, row.names = F, quote = F, sep = "\t")
write.table(results_ba_aa, "~/Documents/Postdoc/Repositories/Baltic_herr_introgression/5.processing_pool_seq_data/results/results_ba_aa_table_1.txt", col.names = T, row.names = F, quote = F, sep = "\t")



