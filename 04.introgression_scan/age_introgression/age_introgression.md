# Estimating the age of introgression

We run this analysis in Uppmax using R code written by Mats Pettersson for the Balsfjord paper (Pettersson et al 2023).

We need to create an input file that is generated during the introgression scan. 

This object is `scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2` and is stored in the following path of our repository:

`4.introgression_scan/intermediate_files/scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2.Rdata`

In the cluster, I simplified the name of this file to `intro_20k_lists.RData`.

We run `age_sim.sh`

~~~bash
#! /bin/bash -l 
#SBATCH -A naiss2024-22-819
#SBATCH -p core -n 1
#SBATCH -t 4-00:00:00
#SBATCH -J Age_sim
#SBATCH -e age_sim_%A_%a.err
#SBATCH -o age_sim_%A_%a.out
#SBATCH --array=2,14
#SBATCH --mail-user=mafalda.ferreira@scilifelab.se
#SBATCH --mail-type=ALL


module load bioinfo-tools R/4.3.1 R_packages/4.3.1

CHR=chr$SLURM_ARRAY_TASK_ID
echo $CHR
R --vanilla --quiet --args ${CHR} < Balsfjord_Age_Sim_UPPMAX.R
~~~

This calls `Balsfjord_Age_Sim_UPPMAX.R`.

~~~R
#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se
#Balsfjord_Age_Sim_UPPMAX.R

argv <- commandArgs(trailingOnly = TRUE)
chr <- argv[1]

print(chr)

require(GenomicRanges)

#save(Ch_v2_linkage_map, Ch_v2_recombination_profile, Ch_v2.0.2_sizes, intro_20k_lists, file = "~/Projects/Herring/data/Balsfjord/age_sim_ref_data.RData")
source("./Balsfjord_Age_Sim_UPPMAX_fun.R")
load("./age_sim_ref_data_baltic.RData")
load("./intro_20k_lists.RData")

#intro_20k_lists <- scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2

# we need to change the chromosome nomenclature
#GenomeInfoDb::seqlevelsStyle(intro_20k_lists$atl) <- "NCBI"
#GenomeInfoDb::seqlevelsStyle(intro_20k_lists$atl_red) <- "NCBI"
#GenomeInfoDb::seqlevelsStyle(intro_20k_lists$pac) <- "NCBI"
#GenomeInfoDb::seqlevelsStyle(intro_20k_lists$pac_red) <- "NCBI"

target_intro_GR <- unlist(intro_20k_lists$atl_red)
pac_target_intro_GR <- unlist(intro_20k_lists$pac_red)

print(paste0("Starting chr: ", chr))
flush.console()
chr_GR <- target_intro_GR[target_intro_GR@seqnames == chr]
pac_chr_GR <- pac_target_intro_GR[pac_target_intro_GR@seqnames == chr]

head(chr_GR)
head(pac_chr_GR)

atl_tmp <- chr_age_estimation(chr = chr, chr_target_GR = chr_GR, gen_map = Ch_v2_linkage_map, gen_prof = Ch_v2_recombination_profile, size_df = Ch_v2.0.2_sizes, n_runs = 500, n_generations = 1e6)
assign(paste("atl_age", chr, sep = "_"),atl_tmp )

pac_tmp <-  chr_age_estimation(chr = chr, chr_target_GR = pac_chr_GR, gen_map = Ch_v2_linkage_map, gen_prof = Ch_v2_recombination_profile, size_df = Ch_v2.0.2_sizes, n_runs = 500, n_generations = 1e6)
assign(paste("pac_age", chr, sep = "_"),pac_tmp )

save(list=c(paste("pac_age", chr, sep = "_"), paste("atl_age", chr, sep = "_")), file = paste0("./", chr, "_age_est_UPPMAX.RData"))
~~~

This script calls the functions in `/Balsfjord_Age_Sim_UPPMAX_fun.R`.

~~~R
#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

simulate_linkage_breakdown_GW <- function(sim_linkage_prof, n_gen = 100){
  f1 <- Vectorize(FUN = function(x,y){min(which(y > x))}, vectorize.args = "x")
  gen_map_max <- max(sim_linkage_prof$interval_dist)
  rec_vec <- runif(n=n_gen, max = max(c(gen_map_max, 100)))
  #Core simulation block
  rec_df <- data.frame(rec_pos = 1:n_gen, lower_bp_tmp = NA, upper_bp_tmp = NA, size = NA, stringsAsFactors = F)
  rec_df$upper_bp_tmp[1] <- sum(sim_linkage_prof$interval_size)
  rec_df$lower_bp_tmp[1] <- 0
  rec_df$rec_pos <- NA 
  rec_df$rec_event <- (rec_vec <= gen_map_max)
  rec_int_vec <- f1(x = rec_vec[rec_df$rec_event], y = sim_linkage_prof$interval_dist) 
  rec_df$rec_interval[rec_df$rec_event] <- rec_int_vec 
  rec_df$rec_event_shift[rec_df$rec_event] <- sample.int(1e6, sum(rec_df$rec_event), replace = T)/1e6
  rec_df$int_start[rec_df$rec_event] <- sim_linkage_prof[rec_int_vec,"BIN_START"]
  rec_df$int_size[rec_df$rec_event] <- sim_linkage_prof[rec_int_vec,"interval_size"]
  rec_df$rec_pos <- rec_df$int_start + round(rec_df$int_size*rec_df$rec_event_shift)
  return(invisible(rec_df))
}

age_per_pos <- function(focus_pos, focus_size, rec_df){
  rec_df$before_focus <- rec_df$rec_pos <= focus_pos
  rec_df$after_focus <- rec_df$rec_pos > focus_pos
  
  rec_df$upper_bp_tmp[which(rec_df$after_focus & rec_df$rec_event)] <- cummin(rec_df$rec_pos[which(rec_df$after_focus & rec_df$rec_event)])
  rec_df$upper_bp <- rep(rec_df$upper_bp_tmp[!is.na(rec_df$upper_bp_tmp)], times = diff(c(which(!is.na(rec_df$upper_bp_tmp)),(length(rec_df$upper_bp_tmp)+1))))
  rec_df$lower_bp_tmp[which(rec_df$before_focus & rec_df$rec_event)] <- cummax(rec_df$rec_pos[which(rec_df$before_focus & rec_df$rec_event)])
  rec_df$lower_bp <- rep(rec_df$lower_bp_tmp[!is.na(rec_df$lower_bp_tmp)], times = diff(c(which(!is.na(rec_df$lower_bp_tmp)),(length(rec_df$lower_bp_tmp)+1))))
  rec_df$size <- rec_df$upper_bp - rec_df$lower_bp
  age_est <- which.min(abs(focus_size - rec_df$size))
  if(all(rec_df$size > focus_size)) age_est <- 2*length(rec_df$size)
  return(age_est)
}

age_pos_vec <- Vectorize("age_per_pos",vectorize.args = c("focus_pos", "focus_size"))

chr_age_estimation <- function(chr, chr_target_GR, gen_map, gen_prof, size_df, n_generations = 2000, n_runs = 100){
  reg_linkage_max <- max(gen_map[gen_map$chr == sub("chr", "", chr),]$all_map)
  phys_size <- size_df[size_df$name == chr,"size"]
  #Evening out the rates across each "step" of the map
  #Adding refinement by using recombination profile scaled by length of map per chromosome
  reg_prof <- gen_prof[gen_prof$CHR ==  sub("chr", "", chr),]
  reg_prof[,"interval_dist"] <- cumsum((reg_prof$RR/sum(reg_prof$RR, na.rm = T)) * reg_linkage_max)
  reg_prof$interval_size <- diff(c(reg_prof$BIN_START, phys_size))
  
  chr_pos_vec <- chr_target_GR@ranges@start + chr_target_GR@ranges@width/2
  chr_size_vec <- chr_target_GR@ranges@width
  age_est_df <- data.frame(hap = names(chr_target_GR), gr = paste(chr_target_GR), stringsAsFactors = F)
  for (i in 1:n_runs){
    chr_rec_df <- simulate_linkage_breakdown_GW(sim_linkage_prof = reg_prof, n_gen = n_generations)
    chr_age_vec <-  age_pos_vec(chr_pos_vec, chr_size_vec, chr_rec_df)
    age_est_df[,paste("age", i, sep = "_")] <- chr_age_vec
  }
  age_est_df$median_age <- apply(age_est_df[, grepl("age_", names(age_est_df))],1, "median")
  return(invisible(age_est_df))
~~~

We also load `age_sim_ref_data_baltic.RData` which contains the Ch_v2.0.2 linkage map, chromosome sizes and recombination profile. This is important for the simulation itself.

The code will output a file with simulated data for all chromosomes, finishing in "_age_est_UPPMAX.RData".

We finally plot data we use the code in `age_introgression.R`

All input and output files described here are available in the repository in the following path:

`4.introgression_scan/age_introgression`