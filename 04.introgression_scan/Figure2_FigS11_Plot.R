# Libraries ####
library(tidyverse)
library(GenomicRanges)

setwd("Manuscript/Figshare/")

# Load Input files ####
load(file="4.introgression_scan/intermediate_files/herring_125.RData")
load("4.introgression_scan/intermediate_files/introgression_scan_2023-11-20.RData")

# Introgression regNULL# Introgression regions: 
intro_reg<-read.table(header=T,"4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")
head(intro_reg)

intro_reg_collapsed<-read.table(header=T,"4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.9regions.collapsed.txt")
head(intro_reg_collapsed)

# Compare with coverage #
load("4.introgression_scan/intermediate_files/scan1_v01_baltic_alt_ref_summary_filter2.Rdata")


# Modified code from HaploDistScan so that we can plot the 
# Introgression scan for Figure 2

dist_df <- scan1_v01_baltic_v_subarctic_dist_alt_ref_20k_df
sample_list=herring_125$sample_list
snp_numbers=scan1_v01_baltic_v_subarctic_alt_ref_diff_SNP_count_vec
min_diff = 20 #20
snp_cutoff = 50 #50
assoc_tresh = 8 #8
assoc_dir = "up"
target_re ="B[MF][0-9]|F[0-9]_HastKar"
eps = 1

target <- sub(target_re, "", names(dist_df)[6])
target_haps <- grep(target_re, sample_list)

hap_ratio_df <- data.frame(hap = sample_list[target_haps], stringsAsFactors= F)

for(bf_hap in 1:dim(hap_ratio_df)[1]){
  #bf_hap <- 1
  bf_ind <- ceiling(bf_hap/2)
  hap_no <- 2 - bf_hap%%2
  
  # this calculates the mean distance of each haplotype to the reference
  hap_ratio_df[bf_hap, "ref_1_mean"] <- mean(dist_df[,2*(bf_hap-1)+6])
  hap_ratio_df[bf_hap, "ref_2_mean"] <- mean(dist_df[,2*(bf_hap-1)+7])
  mean_ratio <- hap_ratio_df[bf_hap, "ref_1_mean"]/hap_ratio_df[bf_hap, "ref_2_mean"]
  
  # ratio between the two distances.
  rv <- (dist_df[,2*(bf_hap-1)+6]+eps)/(dist_df[,2*(bf_hap-1)+7]+eps)
  
  # how many snps are in the window. if bellow 50 SNPs or 
  # there is less than 20 differences between the sequences (in 20kb, that is 0.1%)
  # convert the ratio no NA:
  rv[snp_numbers <= snp_cutoff | (dist_df[,2*(bf_hap-1)+6] < min_diff & dist_df[,2*(bf_hap-1)+7] < min_diff)] <- NA #1
  
  dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")] <- rv
  
  # we only consider regions as associated if the ratio of distance to each target
  # is 8 times the average (expected). 
  # but why 8 times? 
  assoc_filter <- which(rv > assoc_tresh * mean_ratio ) #which(rv > assoc_tresh | rv < 1/assoc_tresh)
  
  # these regions will be regions where the individual is very similar to one of the
  # the references. for this particular individual we see 72  20 kb windows, 
  # and for atl is 1199 20 kb windows 
  hap_ratio_df[bf_hap, "pac"] <- sum(rv > mean_ratio*assoc_tresh, na.rm = T)
  hap_ratio_df[bf_hap, "atl"] <- sum(rv < mean_ratio/assoc_tresh, na.rm = T)
  
  dist_df[assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_filter", sep = "")]<-"assoc"
  dist_df[-assoc_filter,paste(hap_ratio_df[bf_hap, "hap"], "_filter", sep = "")]<-"not_assoc"
  
}

# Now, let's convert into -log10:
for(bf_hap in 1:dim(hap_ratio_df)[1]){
  
  dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_log10", sep = "")] <- log10(dist_df[,paste(hap_ratio_df[bf_hap, "hap"], "_ratio", sep = "")])
  
}

# Modify coverage file ####
coverage<-scan1_v01_baltic_alt_ref_summary_filter2$cov
scan1_v01_baltic_alt_ref_summary_filter2_cov7<-scan1_v01_baltic_alt_ref_summary_filter2$cov[scan1_v01_baltic_alt_ref_summary_filter2$cov[, "cov"] > 7,]
coverage_gr <- GRanges(seqnames=coverage$group_name, IRanges(start=coverage$start,
                                                             end=coverage$end),
                       cov=coverage$cov)

# let's find overlaps...
names(dist_df)[1] <- "seqnames"
names(dist_df)[3] <- "end"
dist_df_gr <- makeGRangesFromDataFrame(dist_df)

overlaps_cov_dist <- findOverlaps(coverage_gr, dist_df_gr )

dist_df$cov <- NA

dist_df[overlaps_cov_dist@to,]$cov <- coverage[overlaps_cov_dist@from,]$cov

dist_df_log <- dist_df %>% 
  select(seqnames, start, end, matches("log10"), cov) %>%
  pivot_longer(cols=c(4:31), names_to="Individual", values_to="Log10")

dist_df_log$Individual <- str_remove(dist_df_log$Individual, pattern="_log10")

dist_df_filter <- dist_df %>% 
  select(seqnames, start, end, matches("filter")) %>%
  pivot_longer(cols=c(4:31), names_to="Individual", values_to="Filter")

dist_df_filter$Individual <- str_remove(dist_df_filter$Individual, pattern="_filter")

# join both
dist_df_to_plot <- left_join(dist_df_log, 
                             dist_df_filter, 
                             by = c("seqnames"="seqnames", 
                                    "start"="start",
                                    "end"="end",
                                    "Individual"="Individual"))

# there should be no "NA" in filter
unique(dist_df_to_plot$Filter)

# remove missing data
dist_df_to_plot_nomiss <- dist_df_to_plot %>% filter(!is.na(Log10))

dist_df_to_plot_nomiss$FinalFilter <- "black"
rows_both <- which(dist_df_to_plot_nomiss$cov > 7 & dist_df_to_plot_nomiss$Filter=="assoc")
rows_assoc <- which(dist_df_to_plot_nomiss$Filter=="assoc" & dist_df_to_plot_nomiss$cov <= 7)
rows_notassoc <- which(dist_df_to_plot_nomiss$Filter=="not_assoc" )

dist_df_to_plot_nomiss[rows_both, ]$FinalFilter <- "red"
dist_df_to_plot_nomiss[rows_assoc, ]$FinalFilter <- "pink"
dist_df_to_plot_nomiss[rows_notassoc, ]$FinalFilter <- "black"

# Manhattan plot in Figure 2 ####
# Manhattan plot ####

dist_df_to_plot_nomiss$Chrom<-as.numeric(str_remove(dist_df_to_plot_nomiss$seqnames,"chr"))
# Order table by chromosome, start and end coordinate of each gene.
# If you only have a single position, such as when you're plotting p-values for 
# SNPs across the genome, just substitute Start or End by the single column
# containing the $Position
dist_df_to_plot_nomiss <- dist_df_to_plot_nomiss[order(dist_df_to_plot_nomiss$Chrom, dist_df_to_plot_nomiss$start, dist_df_to_plot_nomiss$end),]
# Create a column that will have the order in which to plot each data point
dist_df_to_plot_nomiss$n<-1:nrow(dist_df_to_plot_nomiss)
# Key point: Find the center of each chromosome to place axis marks
axisdf <- dist_df_to_plot_nomiss %>% 
  group_by(Chrom) %>% 
  summarize(center=(max(n) + min(n) ) / 2 )

gray <- dist_df_to_plot_nomiss %>% filter(FinalFilter=="pink")
red <- dist_df_to_plot_nomiss %>% filter(FinalFilter=="red")
black <- dist_df_to_plot_nomiss %>% filter(FinalFilter=="black")


# Manhattan plot:
introgression_scan_baltic <- black %>%
  # Here, you give "n" as x, and not the coordinates inside each chromosome.
  # Color points according to chromosome
  # y is your stat of interest.
  ggplot(aes(x=n, y=Log10, color = as.factor(Chrom)))+
  geom_point(size=0.2) +
  scale_color_manual(values = rep(c("black","gray"), 26 )) +
  
  # Here I color outlier genes in red:
  geom_point(data = gray, aes(x=n, y=Log10), colour = "pink", size=0.3) +
  geom_point(data = red, aes(x=n, y=Log10), colour = "red", size=0.3) +
  
  # Modify the axis to contain chromosome names in order.
  scale_x_continuous(label = axisdf$Chrom, breaks= axisdf$center , expand = c(0.01, 0.01)) +
  labs(x="Chromosome", y=expression(paste(log[10], " distance ratio")))+
  theme_classic() +
  theme(legend.position="none", panel.border = element_blank(),
        axis.title = element_text(size=6),
        axis.text = element_text(size=6),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(introgression_scan_baltic, 
       filename="figures/scan1_v01_baltic_alt_ref_summary_filter2_cov7_red_pink_fixed.png",
       units="mm", 
       height = 38, 
       width = 170)

ggsave(introgression_scan_baltic, 
       filename="figures/scan1_v01_baltic_alt_ref_summary_filter2_cov7_red_pink_fixed.pdf",
       units="mm", 
       height = 38, 
       width = 170)


# Plots in Supplementary Figure 11 ####
# Load ChiSquare Test Results #
load("5.processing_pool_seq_data/chisquare_results/baltic_spring_vs_atlantic_spring_han_pops_chiseq.output.RData")
baltic_spring_vs_atlantic_spring_chiseq$Log10 <- -log10(baltic_spring_vs_atlantic_spring_chiseq$chisq_p)

# Autumn:
load("5.processing_pool_seq_data/chisquare_results/baltic_autumn_vs_atlantic_autumn_han_pops_chiseq.output.RData")
baltic_autumn_vs_atlantic_autumn_chiseq$Log10 <- -log10(baltic_autumn_vs_atlantic_autumn_chiseq$chisq_p)

padding <- 1
padding_mb <- 1e6

list_plots <- list()

for(i in 1:nrow(intro_reg_collapsed)){

  if(i==2){
    
    # Define regions
    c<-intro_reg_collapsed[i, "seqnames"]
    s<-intro_reg_collapsed[i, "start"]
    e<-intro_reg_collapsed[i, "end"]
    chr<-paste0("chr", c)
    
    candidate_matrix <- intro_reg %>% filter(seqnames == chr & start >= s & end <= e)
    
    tmpBSAS <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    tmpBAAA <- baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    
    tmp <- dist_df_to_plot_nomiss %>% filter(seqnames==chr & 
                                               start >= s-padding_mb & 
                                               end <= e + padding_mb)
    
    tmp$mid <- (tmp$start + tmp$end)/2
    
    tmpBSAS_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBSAS, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      xlim(s/1e6 - padding, 32.6)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    tmpBAAA_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBAAA, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      scale_y_reverse() + 
      scale_x_continuous(position = "top", limits = c(s/1e6 - padding, 32.6))+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size=6),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=6))
    
    
    
    dist_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=-4,ymax=4),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmp, aes(x=mid/1e6, y=Log10, color=FinalFilter), size=0.5)+
      scale_color_identity()+
      #geom_point(data=tmp, aes(x=mid/1e6, y=Log10, shape=FinalFilter, fill), color=NA, size=1)+
      #scale_fill_manual(values=c("purple", "black"))+
      xlab(paste0(c," position (Mb)"))+
      xlim(s/1e6 - padding, 32.6)+
      ylab(expression(paste(log[10], " (distance ratio)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    
    region_plot <- cowplot::plot_grid(dist_plot, tmpBSAS_plot, tmpBAAA_plot, ncol=1)
    

    list_plots[[i]] <- region_plot 
    
    
  }else if(i==8){
    
    # Define regions
    c<-intro_reg_collapsed[i, "seqnames"]
    s<-intro_reg_collapsed[i, "start"]
    e<-intro_reg_collapsed[i, "end"]
    chr<-paste0("chr", c)
    
    candidate_matrix <- intro_reg %>% filter(seqnames == chr & start >= s & end <= e)
    
    tmpBSAS <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    tmpBAAA <- baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    
    tmp <- dist_df_to_plot_nomiss %>% filter(seqnames==chr & 
                                               start >= s-padding_mb & 
                                               end <= e + padding_mb)
    
    tmp$mid <- (tmp$start + tmp$end)/2
    
    tmpBSAS_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=150),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBSAS, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      xlim(s/1e6 - padding, e/1e6 + padding)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    tmpBAAA_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=150),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBAAA, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      scale_y_reverse() + 
      scale_x_continuous(position = "top", limits = c(s/1e6 - padding, e/1e6 + padding))+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size=6),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=6))
    
    
    
    dist_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=-4,ymax=4),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmp, aes(x=mid/1e6, y=Log10, color=FinalFilter), size=0.5)+
      scale_color_identity()+
      #geom_point(data=tmp, aes(x=mid/1e6, y=Log10, shape=FinalFilter, fill), color=NA, size=1)+
      #scale_fill_manual(values=c("purple", "black"))+
      xlab(paste0(c," position (Mb)"))+
      xlim(s/1e6 - padding, e/1e6 + padding)+
      ylab(expression(paste(log[10], " (distance ratio)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    
    region_plot <- cowplot::plot_grid(dist_plot, tmpBSAS_plot, tmpBAAA_plot, ncol=1)
    
    #directory="~/Documents/Postdoc/Project_Herring/Introgression/Pool_data/"
    #ggsave(region_plot, file=paste0(directory,"figures/combined_intro_scan_chiseq_",c,"_",s,"_",e,".png"),
    #       units="mm",
    #       height=72.5,
    #       width=72.5)
    
    list_plots[[i]] <- region_plot
    
  }else{
    
    
    # Define regions
    c<-intro_reg_collapsed[i, "seqnames"]
    s<-intro_reg_collapsed[i, "start"]
    e<-intro_reg_collapsed[i, "end"]
    chr<-paste0("chr", c)
    
    candidate_matrix <- intro_reg %>% filter(seqnames == chr & start >= s & end <= e)
    
    tmpBSAS <- baltic_spring_vs_atlantic_spring_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    tmpBAAA <- baltic_autumn_vs_atlantic_autumn_chiseq %>% filter(CHROM==chr & 
                                                                    POS >= s-padding_mb & 
                                                                    POS <= e + padding_mb)
    
    
    tmp <- dist_df_to_plot_nomiss %>% filter(seqnames==chr & 
                                               start >= s-padding_mb & 
                                               end <= e + padding_mb)
    
    tmp$mid <- (tmp$start + tmp$end)/2
    
    tmpBSAS_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBSAS, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      xlim(s/1e6 - padding, e/1e6 + padding)+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    tmpBAAA_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=80),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmpBAAA, aes(x=POS/1e6, y=Log10, color="black"), size=0.5)+
      scale_color_identity()+
      scale_y_reverse() + 
      scale_x_continuous(position = "top", limits = c(s/1e6 - padding, e/1e6 + padding))+
      xlab(paste0(c," position (Mb)"))+
      ylab(expression(paste(log[10], "(Pvalue)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size=6),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size=6))
    
    
    
    dist_plot <- ggplot()+     
      geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=-4,ymax=4),
                fill="azure3",color=NA,alpha=1,inherit.aes = FALSE)+
      geom_point(data=tmp, aes(x=mid/1e6, y=Log10, color=FinalFilter), size=0.5)+
      scale_color_identity()+
      #geom_point(data=tmp, aes(x=mid/1e6, y=Log10, shape=FinalFilter, fill), color=NA, size=1)+
      #scale_fill_manual(values=c("purple", "black"))+
      xlab(paste0(c," position (Mb)"))+
      xlim(s/1e6 - padding, e/1e6 + padding)+
      ylab(expression(paste(log[10], " (distance ratio)")))+
      theme_classic()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(size=6),
            axis.title.y = element_text(size=6),
            axis.text.x = element_text(size=6),
            axis.text.y = element_text(size=6))
    
    
    region_plot <- cowplot::plot_grid(dist_plot, tmpBSAS_plot, tmpBAAA_plot, ncol=1)
    
    #directory="~/Documents/Postdoc/Project_Herring/Introgression/Pool_data/"
    #ggsave(region_plot, file=paste0(directory,"figures/combined_intro_scan_chiseq_",c,"_",s,"_",e,".png"),
    #       units="mm",
    #       height=72.5,
    #       width=72.5)
    
    list_plots[[i]] <- region_plot
    
  }
  
}



margin = theme(plot.margin = unit(rep(0.5,4), "cm"))

plots1.6 <- gridExtra::grid.arrange(grobs = lapply(list_plots[c(1:6)], "+", margin), ncol=2)
plots7.9 <- gridExtra::grid.arrange(grobs = lapply(list_plots[c(7:9)], "+", margin), ncol=2)


directory="directory/"
ggsave(plots1.6, file=paste0(directory,"figures/combined_daf_introgression_scan_intro_reg_1_to_6.pdf"),
       units="mm",
       height=200,
       width=180)

ggsave(plots7.9, file=paste0(directory,"figures/combined_daf_introgression_scan_intro_reg_7_to_9.pdf"),
       units="mm",
       height=200,
       width=180)


ggsave(plots1.6, file=paste0(directory,"figures/combined_daf_introgression_scan_intro_reg_1_to_6.png"),
       units="mm",
       height=200,
       width=180)

ggsave(plots7.9, file=paste0(directory,"figures/combined_daf_introgression_scan_intro_reg_7_to_9.png"),
       units="mm",
       height=200,
       width=180)

