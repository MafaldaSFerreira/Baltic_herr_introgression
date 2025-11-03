# Code to plot Twisst results in Figures 2 and Supplementary Figures 6, 7 and 8

# Introgression regions ####
intro_reg_collapsed<-read.table(header=T,"scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.9regions.collapsed.txt")

intro_reg<-read.table(header=T,"scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")


# TWISST: ####
# Vancouver as outgroup
weights <- read.table("/Users/mafaldaferreira/Library/CloudStorage/Dropbox/Mac (2)/Documents/Postdoc/Project_Herring/Introgression/Twisst/2023-12-04_twisst_new/output.weights_100_soft.csv", skip = 3, header=T)
windows <- read.table("~/Dropbox/Mac (2)/Documents/Postdoc/Project_Herring/Introgression/Twisst/2023-12-04_twisst_new/output.phyml_bionj.w100_soft.data.tsv", header=T)

# Use sprat output:
weights <- read.table("/Users/mafaldaferreira/Library/CloudStorage/Dropbox/Mac (2)/Documents/Postdoc/Project_Herring/Introgression/Twisst/2025-02-24_TwisstSprat/results/output.wg.phyml.w100.Mi10.weights.csv.gz", skip = 3, header=T)
windows <- read.table("/Users/mafaldaferreira/Library/CloudStorage/Dropbox/Mac (2)/Documents/Postdoc/Project_Herring/Introgression/Twisst/2025-02-24_TwisstSprat/results/output.wg.phyml.w100.Mi10.data.tsv", header=T)

weights = weights / apply(weights, 1, sum)
good_rows = which(is.na(apply(weights, 1, sum)) == F)
weights <- weights[good_rows,]
windows <- windows[good_rows,]

twisst_df <-cbind(windows,weights)

# Suplementary Figure 6 - Summary Plot ####
head(twisst_df)

twisst_df %>%
  pivot_longer(cols=c(7:9)) %>%
  ggplot()+
  geom_bar(aes(x=value, color=name))


twisst_summary_plot <- twisst_df %>% pivot_longer(cols=7:9) %>%
  group_by(name) %>%
  summarise(mean_value=mean(value)) %>%
  ggplot(aes(y=mean_value, x=name, fill=name))+
  geom_col()+
  theme_classic()+
  scale_fill_manual(values=c("#882E72", "#E8601C","#5289C7" ))+
  ylab("Average weighting")+
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position="none")
  
  
intro_reg_gr <- GRanges(seqnames=intro_reg$seqnames, IRanges(start=intro_reg$start, end=intro_reg$end))

# Convert tables to GRanges:
twisst_df_gr<-GRanges(seqnames = twisst_df$scaffold, 
                            IRanges(start = twisst_df$start, end = twisst_df$end),
                            topo1=twisst_df$topo1_frac,
                      topo2=twisst_df$topo2_frac,
                      topo3=twisst_df$topo3_frac)

# Find overlaps between both:
overlaps<-findOverlaps(intro_reg_gr, twisst_df_gr)

twisst_df$type <- "NotIntrogressed"

twisst_df[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

twisst_weights_by_topo_plot <- twisst_df %>% pivot_longer(cols=7:9) %>%
  ggplot()+
  geom_violin(aes(y=value, x=type), width=0.2, )+
  geom_boxplot(aes(y=value, x=type), width=0.08, outlier.shape = NA)+
  theme_classic()+
  facet_wrap(facets = ~name)+
  scale_x_discrete(label=c("Introgressed \nregions", "Genome without \nintrogressed regions", "Introgressed \nregions", "Genome without \nintrogressed regions"))+
  ylab("Topology weight")+
  theme(axis.title.x = element_blank(),
        strip.text  = element_text(size=6, color="black"),
        axis.text.x = element_text(size=6, color="black", angle=45, hjust=1),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"),
        legend.position = "none")

signif(t.test(topo2 ~ type, data=twisst_df, alternative="greater")$p.value, 2)
# 9.8e-29

all_plot <- grid.arrange(twisst_summary_plot, twisst_weights_by_topo_plot, nrow=2, ncol=2, layout_matrix=cbind(c(0,1),c(2,1)))

ggsave(all_plot, filename="figures/scan1_intro_reg_vs_WG_boxplot_Twisst_2025-06-26.pdf", units="mm", width = 90, height = 90)


# Supplementary Figures 7 and 8 ####

## Comparing sprat or vancouver as outgroup ####

# Let's plot the introgre collapsed, but highlight the introgression blocks
# intro_reg and intro_reg_collapsed.


padding <- 0.25
padding_mb <- 2.5e5

# If you want to make a list of plots for A4 plotting:
list_plots <- list()

for(i in 1:nrow(intro_reg_collapsed)){

  # Define regions
  c<-intro_reg_collapsed[i, "seqnames"]
  s<-intro_reg_collapsed[i, "start"]
  e<-intro_reg_collapsed[i, "end"]
  chr<-paste0("chr", c)
  
  candidate_matrix <- intro_reg %>% filter(seqnames == chr & start >= s & end <= e)
  
  tmp_df <- twisst_df %>% filter(scaffold == chr & start >= s - padding_mb & end <= e + padding_mb)
  
  twisst_plot <- tmp_df %>%
    ggplot()+
    geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1),
              fill="azure3",color=NA,alpha=0.5,inherit.aes = FALSE)+
    geom_line(aes(x=mid/1000000,y=topo1_frac),col="#882E72",linewidth=1)+
    geom_line(aes(x=mid/1000000,y=topo2_frac),col="#E8601C",linewidth=1)+
    geom_line(aes(x=mid/1000000,y=topo3_frac),col="#5289C7",linewidth=1)+
    xlim(s/1e6-padding, e/1e6+padding)+
    ylab("topology support (%)") +
    xlab(paste0("Chr ",c," position (Mb)"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_text(colour="black",size=6),
          axis.text.x = element_text(colour="black",size=6),
          axis.title.y = element_text(colour="black",size=6),
          axis.title.x = element_text(colour="black",size=6))
  
  list_plots[[i]] <- twisst_plot
  
}

# Outgroup sprat:
all_sprat_plots <- grid.arrange(grobs = list_plots, ncol=2)
ggsave(all_sprat_plots, filename="~/Documents/Postdoc/Project_Herring/Introgression/Twisst/2025-02-24_TwisstSprat/figures/all_intro_reg.twisst_OutSprat_Mi10_introg_reg_collapsed.pdf", units="mm", height = 150, width = 210)

# Outgroup pacific herring:
all_pacific_plots <- grid.arrange(grobs = list_plots, ncol=2)
ggsave(all_pacific_plots, filename="~/Documents/Postdoc/Project_Herring/Introgression/Twisst/2023-12-04_twisst_new/figures_2025-03-03/all_intro_reg.twisst_OutVancouver_introg_reg_collapsed.pdf", units="mm", height = 150, width = 210)


# Figure 2 Zoom In into Chromosome 10 Regions ####

# 2025-07-14: Figure 2 Twisst plot with entire chromosome 10 ####
#weights <- "results/output.wg.phyml.w100.Mi10.weights.csv"
#windows <- "results/output.wg.phyml.w100.Mi10.data.tsv"
#twisst_data <- import.twisst(weights, windows)

df<-cbind(window,weights)
twisst_df<-mutate(twisst_df,
           topo1_frac=topo1/(topo1+topo2+topo3),
           topo2_frac=topo2/(topo1+topo2+topo3),
           topo3_frac=topo3/(topo1+topo2+topo3))

xlim_min=s-50000
xlim_max=e+50000
plot.weights(weights_dataframe=twisst_data$weights[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
             line_cols=topo_cols, fill_cols=F, stacked=F)

c=10
chr="chr10"

  candidate_matrix <- intro_reg %>% filter(seqnames == chr & start==26300001)
  
  twisst_df %>%
    filter(scaffold==chr) %>%
    ggplot()+
    geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1.2),
              fill="azure3",color=NA,inherit.aes = FALSE)+
    geom_line(aes(x=mid/1000000,y=topo3),col="#5289C7",size=0.5)+
    geom_line(aes(x=mid/1000000,y=topo1),col="#882E72",size=0.5)+
    geom_line(aes(x=mid/1000000,y=topo2),col="#E8601C",size=0.5)+
    #xlim(s/1e6-padding, e/1e6+padding)+
    ylab("topology support (%)") +
    xlab(paste0(chr," position (Mb)"))+
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.y = element_text(colour="black",size=6),
          axis.text.x = element_text(colour="black",size=6),
          axis.title.y = element_text(colour="black",size=6),
          axis.title.x = element_text(colour="black",size=6))
  
 
ggsave(twisst_plot, filename="~/Documents/Postdoc/Project_Herring/Introgression/Twisst/2025-02-24_TwisstSprat/figures/chr10_overall.pdf", height = 40, width = 180, units = "mm")
 
twisst_df_long <- twisst_df %>% pivot_longer(col=c(7:9), names_to = "topology", values_to="weights")

twisst_df_long$topology <- factor(twisst_df_long$topology, levels=c("topo3", "topo2", "topo1"))

twisst_plot_area <-twisst_df_long %>% 
  filter(scaffold=="chr10") %>%
  ggplot(aes(x=mid/1e6, y=weights, fill=topology)) + 
  ylab("topology support (%)") +
  xlab(paste0(chr," position (Mb)"))+
  geom_area()+
  scale_fill_manual(values=c("#5289C7","#E8601C","#882E72"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(colour="black",size=6),
        axis.text.x = element_text(colour="black",size=6),
        axis.title.y = element_text(colour="black",size=6),
        axis.title.x = element_text(colour="black",size=6),
        legend.position = "none")

xlim_min <- min(twisst_df_long[twisst_df_long$scaffold=="chr10",]$mid)
xlim_max <- max(twisst_df_long[twisst_df_long$scaffold=="chr10",]$mid)

candidates <- ggplot()+
  geom_rect(data=candidate_matrix,mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=0,ymax=1),
            fill="azure3",color=NA,inherit.aes = FALSE)+
    xlim(xlim_min/1e6, xlim_max/1e6)+
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 

  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # theme(panel.border = element_blank()) +
  # theme(axis.line = element_line(colour = "black")) +
  # theme(axis.title.x = element_blank()) +
  # theme(axis.text.y = element_text(colour="black",size=6),
  #       axis.text.x = element_text(colour="black",size=6),
  #       axis.title.y = element_text(colour="black",size=6),
  #       axis.title.x = element_text(colour="black",size=6),
  #       legend.position = "none")


plot_figure2 <- cowplot::plot_grid(candidates, twisst_plot_area, nrow=2,align="v")  

ggsave(plot_figure2, file="figures/chr10_stacked_weights_sprat_outgroup_fig2.pdf",
       units = "mm",
       height = 80,
       width = 180)
