#Longer run, from UPPMAX
atl_age_vec <- numeric()
pac_age_vec <- numeric()
setwd("~/Documents/Postdoc/Project_Herring/Introgression/age_introgression")

n<-c(1,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26)

for(chr in paste0("chr",n)){ 
  #chr <- "chr10"
  load(paste0("./",chr, "_age_est_UPPMAX.RData"))
  #load(paste0("chr10_age_est_UPPMAX.RData"))
  tmp_atl <- get(paste0("atl_age_", chr))
  tmp_atl_age <- tmp_atl$median_age
  names(tmp_atl_age) <- tmp_atl$gr
  atl_age_vec <- c(atl_age_vec, tmp_atl_age)
  
  tmp_pac <- get(paste0("pac_age_", chr))
  tmp_pac_age <- tmp_pac$median_age
  names(tmp_pac_age) <- tmp_pac$gr
  pac_age_vec <- c(pac_age_vec, tmp_pac_age)
}

#pdf(file = "~/Projects/Herring/doc/Balsfjord/age_estimate/GW_age_hist_UPPMAX.pdf")
age_co <- 1e5
atl_col <- col2rgb("red")
atl_col <- rgb(t(atl_col),alpha = 150, maxColorValue = 255)
hist(atl_age_vec[atl_age_vec < age_co], breaks = 100, main = "Atlantic introgressions", xlab = "Age (generations)", col = atl_col)
#hist(atl_age_vec, breaks = 100, main = "Atlantic introgressions", xlab = "Age (generations)", col = atl_col)

sum(atl_age_vec < age_co)
sum(!atl_age_vec < age_co)

pac_col <- col2rgb("blue",)
pac_col <- rgb(t(pac_col), alpha = 150, maxColorValue = 255)

hist(pac_age_vec[pac_age_vec < age_co], breaks = 100, main = "Pacific introgressions", xlab = "Age (generations)", col = pac_col, add=T)
sum(!pac_age_vec < age_co)
sum(pac_age_vec < age_co)

text(x = 6.0e4, y = 380, labels = paste("\"Atlantic\" loci below 1e5:",sum(atl_age_vec < age_co)), pos = 2, cex = 1)
text(x = 6.0e4, y = 380, labels = paste("\"Atlantic\" loci above 1e5:",sum(!atl_age_vec < age_co)), pos = 4, cex = 1)
text(x = 6.0e4, y = 400, labels = paste("\"Pacific\" loci below 1e5:",sum(pac_age_vec < age_co)), pos = 2, cex = 1)
text(x = 6.0e4, y = 400, labels = paste("\"Pacific\" loci above 1e5:",sum(!pac_age_vec < age_co)), pos = 4, cex = 1)
text(x = 6.0e4, y = 360, labels = paste0("P(wilcox) = ", signif(wilcox.test(x = pac_age_vec, y = atl_age_vec)$p.value, digits = 2)))

abline(v = median(atl_age_vec), col = "red", lwd = 2)
abline(v = median(pac_age_vec), col = "blue", lwd = 2)
legend(x = "right", legend = c("Atlantic", "Pacific", "Atlantic median", "Pacific median"), fill = c(atl_col, pac_col, NA, NA), border = c(1,1,NA, NA), col = c(NA, NA, "red", "blue"), lwd = c(NA, NA, 2, 2), merge = T, seg.len = 1, cex = 1.1)
dev.off()

# 2025-01-19
median(pac_age_vec)
# 6938.5
median(atl_age_vec)
# 8008.25

library(tidyverse)

## FOR PAPER ###
histogram <- data.frame(pac_age_vec) %>%
  filter(pac_age_vec < age_co) %>%
  ggplot()+
  geom_histogram(aes(x=pac_age_vec), bins = 50, color="white")+ 
  theme_classic()+
  xlab("Age (generations)")+
  ylab("Frequency")+
  theme(axis.title.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6))
  
ggsave(histogram, filename="figures/histogram_age_simulation.pdf", width=90, height=45, units = "mm")

#extracting the recurring introgressions
atl_rec_age <- numeric()
for(i in 1:length(atl_reg_region_results)){
  tmp_atl_rec_GR <- atl_reg_region_results[[i]]$reg_intros
  target_chr <- as.character(seqnames(tmp_atl_rec_GR)[1])
  tmp_atl_age_est <- get(paste0("atl_age_", target_chr))
  atl_age_est_hits <- findOverlaps(tmp_atl_rec_GR, GRanges(tmp_atl_age_est$gr), type = "equal")
  tmp_atl_rec_age <- tmp_atl_age_est$median_age[unique(atl_age_est_hits@to)]
  atl_rec_age <- c(atl_rec_age, tmp_atl_rec_age)
}
median(atl_rec_age)

pac_rec_age <- numeric()
for(i in 1:length(pac_reg_region_results)){
  tmp_pac_rec_GR <- pac_reg_region_results[[i]]$reg_intros
  target_chr <- as.character(seqnames(tmp_pac_rec_GR)[1])
  tmp_pac_age_est <- get(paste0("pac_age_", target_chr))
  pac_age_est_hits <- findOverlaps(tmp_pac_rec_GR, GRanges(tmp_pac_age_est$gr), type = "equal")
  tmp_pac_rec_age <- tmp_pac_age_est$median_age[unique(pac_age_est_hits@to)]
  pac_rec_age <- c(pac_rec_age, tmp_pac_rec_age)
}
median(pac_rec_age)

#plotting some old regions
head(atl_age_vec[atl_age_vec > 2e5])
atl_age_vec[atl_age_vec > 2e5][40:50]
old_regs <- sub("(chr[0-9]{1,2})[:]([0-9]+).+", "\\1_\\2", names(atl_age_vec[atl_age_vec > 2e5]))
old_windows <- which(paste(h79_20k_w$scaffold, h79_20k_w$start, sep = "_") %in% old_regs )
which(h79_20k_w$scaffold == "chr1" & h79_20k_w$start == 18360001)
#[1] 919
which(h79_20k_w$scaffold == "chr7" & h79_20k_w$start == 13080001)
#[1] 10356
which(h79_20k_w$scaffold == "chr8" & h79_20k_w$start == 11260001)
#[1] 11815

plot_top_regions_2(top_pos = old_windows, win_df = h79_20k_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/Balsfjord/old_regions/", snp_col = "blue")



#Support functions
region_age_estimation <- function(region_GR, intro_list = intro_20k_lists$atl_red, gen_map = Ch_v2_linkage_map, gen_prof = Ch_v2_recombination_profile, size_df = Ch_v2.0.2_sizes,  n_generations = 2000, n_runs = 100, pdf_file = NULL){
  introgression_pos <- region_GR@ranges@start + region_GR@ranges@width/2
  focal_GR <- GRanges(seqnames = region_GR@seqnames, ranges = IRanges(start = introgression_pos, end = introgression_pos))
  reg_linkage_max <- max(gen_map[gen_map$chr == sub("chr", "", as.character(focal_GR@seqnames)),]$all_map)
  phys_size <- size_df[size_df$name == as.character(focal_GR@seqnames),"size"]
  
  
  #Evening out the rates across each "step" of the map
  #reg_linkage_map <- cbind(reg_linkage_map[which(diff(reg_linkage_map[,"all_map"]) > 0) +1 ,], diff(reg_linkage_map[,"all_map"])[which(diff(reg_linkage_map[,"all_map"]) > 0)])
  #reg_linkage_map$HiC_pos[dim(reg_linkage_map)[1]] <- phys_size #Stretching the final interval to the end of the chromosome
  #reg_linkage_map <- reg_linkage_map[,c(2,5,7)]
  #names(reg_linkage_map)[3] <- "interval_dist"
  #reg_linkage_map$interval_size <- diff(c(0,reg_linkage_map$HiC_pos))
  #Adding refinement by using recombination profile scaled by length of map per chromosome
  reg_prof <- gen_prof[gen_prof$CHR ==  sub("chr", "", as.character(focal_GR@seqnames)),]
  reg_prof[,"interval_dist"] <- cumsum((reg_prof$RR/sum(reg_prof$RR, na.rm = T)) * reg_linkage_max)
  reg_prof$interval_size <- diff(c(reg_prof$BIN_START, phys_size))
  
  block_size_df <- data.frame(dummy = 1:n_generations)
  rec_pos_df <- data.frame(dummy = 1:n_generations)
  for (i in 1:n_runs){
    tm_sim_df <- simulate_linkage_breakdown(sim_linkage_prof = reg_prof, focus_pos = introgression_pos,n_gen = n_generations, draw_plot = F, phys_size = phys_size)
    block_size_df[,i] <- tm_sim_df$size
    rec_pos_df[,i] <- tm_sim_df$rec_pos
  }
  target_intros <- GRanges()
  for(i in 1:length(intro_list)){
    target_intros[names(intro_list)[i]] <- intro_list[[i]][findOverlaps(focal_GR, intro_list[[i]])@to]
  }
  median_vec <- (apply(block_size_df,1, "median"))
  age_est_vec <- array()
  for(i in 1:length(target_intros)){
    age_est_vec[i] <- which.min(abs(target_intros@ranges@width[i] - median_vec))
  }
  
  if(!is.null(pdf_file)) pdf(pdf_file)
  simulate_linkage_breakdown(sim_linkage_prof = reg_prof, focus_pos = introgression_pos,n_gen = n_generations, draw_plot = T, phys_size = phys_size)
  plot(x = 1:n_generations, y = rowSums(block_size_df)/n_runs/1e6, xlab = "Generations after introgression", ylab = "Block size", pch = 16, col = "firebrick", main = "Results from 100 simulations", type  ="n", ylim = c(0,1.5))
  points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[2,], pch = 16, col = "grey50")
  points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[3,], pch = 17, col = "darkorchid")
  points(x = 1:n_generations, y = (apply(block_size_df,1, "summary")/1e6)[5,], pch = 18, col = "grey50")
  if(!is.null(pdf_file)) dev.off()
  
  return(list(block_size = block_size_df, age_est_vec =  age_est_vec, reg_intros = target_intros))
}

#Version that uses the recombiation profile instead, calibrated by the lenght of the linkage map 
simulate_linkage_breakdown <- function(sim_linkage_prof,focus_pos, n_gen = 100, phys_size = NULL, draw_plot = F){
  rec_df <- data.frame(rec_pos = 1:n_gen, lower_bp = NA, upper_bp = NA, size = NA)
  lower_bp <- 0
  upper_bp <- sum(sim_linkage_prof$interval_size)
  gen_map_max <- max(sim_linkage_prof$interval_dist)
  rec_vec <- runif(n=n_gen, max = max(c(gen_map_max, 100)))
  for(i in 1:n_gen){
    if(rec_vec[i] <= max(gen_map_max)){
      rec_interval <- min(which(sim_linkage_prof$interval_dist > rec_vec[i]))
      rec_pos <- sim_linkage_prof$BIN_START[rec_interval] + sample(0:sim_linkage_prof$interval_size[rec_interval], 1)
      if(rec_pos < focus_pos & rec_pos > lower_bp) lower_bp <- rec_pos
      if(rec_pos > focus_pos & rec_pos < upper_bp) upper_bp <- rec_pos
    } else{
      rec_pos <- NA
    }
    rec_df[i,] <- c(rec_pos, lower_bp, upper_bp, upper_bp-lower_bp)
  }
  if(draw_plot){
    plot(x = c(0, phys_size), y = c(0, n_gen), xlab = "Position", ylab = "Generation", type = "n", axes = F)
    axis(1)
    axis(2, at = seq(from = 0, to = n_gen, by = 10), labels = rev(seq(from = 0, to = n_gen, by = 10)))
    segments(x0 = rec_df[,"lower_bp"], x1 = rec_df[,"upper_bp"], y0 = rev(1:n_gen), col = "blue", lwd = 2)
    points(x = rec_df[,"rec_pos"], y = rev(1:n_gen), pch = 4, col = "red")
  }
  return(invisible(rec_df))
}

simulate_linkage_breakdown_original <- function(sim_linkage_map, focus_pos, n_gen = 100, phys_size = NULL, draw_plot = F){
  rec_df <- data.frame(rec_pos = 1:n_gen, lower_bp = NA, upper_bp = NA, size = NA)
  lower_bp <- 0
  upper_bp <- max(sim_linkage_map$HiC_pos)
  
  rec_vec <- runif(n=n_gen, max = max(c(sim_linkage_map$all_map, 100)))
  for(i in 1:n_gen){
    if(rec_vec[i] <= max(sim_linkage_map$all_map)){
      rec_interval <- min(which(sim_linkage_map$all_map > rec_vec[i]))
      rec_pos <- sim_linkage_map$HiC_pos[rec_interval] - sample(0:sim_linkage_map$interval_size[rec_interval], 1)
      if(rec_pos < focus_pos & rec_pos > lower_bp) lower_bp <- rec_pos
      if(rec_pos > focus_pos & rec_pos < upper_bp) upper_bp <- rec_pos
    } else{
      rec_pos <- NA
    }
    rec_df[i,] <- c(rec_pos, lower_bp, upper_bp, upper_bp-lower_bp)
  }
  if(draw_plot){
    plot(x = c(0, phys_size), y = c(0, n_gen), xlab = "Position", ylab = "Generation", type = "n", axes = F)
    axis(1)
    axis(2, at = seq(from = 0, to = n_gen, by = 10), labels = rev(seq(from = 0, to = n_gen, by = 10)))
    segments(x0 = rec_df[,"lower_bp"], x1 = rec_df[,"upper_bp"], y0 = rev(1:n_gen), col = "blue", lwd = 2)
    points(x = rec_df[,"rec_pos"], y = rev(1:n_gen), pch = 4, col = "red")
  }
  return(invisible(rec_df))
}
