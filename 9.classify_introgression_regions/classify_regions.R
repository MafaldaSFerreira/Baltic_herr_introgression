# Classifying Introgression Regions based on dxy and positive selection results
# These results are displayed in Supplementary Figure 13
# And in Supplementary Table 7 

# We start with dxy results from pixy which have been generated for group of Baltic spring-spawning homozygotes and White Sea. This is documented in step7.pixy/plot_dxy.R
# We also use the xpEHH results as documented in step 8.


## So I need to find how I classified the introgression regions based on dxy and xpEHH. Look at the code before.
##
setwd("/Users/mafaldaferreira/Dropbox/Mac/Documents/Postdoc/Project_Herring/Introgression/pixy/2025-10-18_results_homozygotes_maf5")

# First let's read the dxy files

dxy_files<-list.files(pattern = "out_dxy")
regions <- str_split_fixed(dxy_files, pattern="_", 3)[,c(1,2)]

# Great, now, we probably need to read in the selection files.
xpehh_results<-read.table("xpehh_group1+2.out", skip=1)
xpehh_results$CHROM<-str_split_fixed(xpehh_results$V2, "_", 2)[,1]
xpehh_results$POS<-as.numeric(str_split_fixed(xpehh_results$V2, "_", 2)[,2])
xpehh_results$V2<-NULL

colnames(xpehh_results)<-c("Index", "Freq", "iHH_A1", "iHH_B1", "iHH_P1", "XPEHH", "std_XPEHH", "CHROM", "POS")

# Score results
xpehh_results_sig<-xpehh_results  %>% 
  arrange(desc(abs(std_XPEHH))) %>%   # Sort SNPs by XPEHH score (descending)
  mutate(
    rank=row_number(),           # number of SNPs with higher scores
    rank.perc=rank/n(),          # fraction of SNPs with higher score
    rank.log=-log10(rank.perc)   # P-value
  )




## CLASSIFY INTROGROSSION REGIONS ####

# Classify the windows as introgressed if:
# dxy BS vs WS < 0.05
# dxy BA vs WS > 0.05
# region is under selection std_XPEHH > 2
options(scipen = 999)
classified_regions<-data.frame(chr=NA, start=NA, end=NA, dxy_BS_vs_WS=NA, dxy_pcent_BS_vs_WS=NA, dxy_BA_vs_WS=NA, dxy_pcent_BA_vs_WS=NA, XPEHH_mean=NA, std_XPEHH_mean=NA,
                               VeryGood=NA, Good=NA)

for(i in 1:nrow(regions)){
  
  chr<-regions[i,1]
  s<-as.numeric(str_split_fixed(regions[i,2], pattern="-", 2)[1,1])
  e<-as.numeric(str_split_fixed(regions[i,2], pattern="-", 2)[1,2])
  c<-str_remove(chr, "chr")
  
  print(c(chr, s, e))
  
  # These new files are called "WG". They are called maf1 but they are maf5!!
  table_dxy_homozygotes <-read.table(header=T, paste0(chr,"_",s,"-",e,"_Baltic_Spring_WhiteSea_BalticAutumn_popfile.wg.maf1.20kb.popgenpixy.out_dxy.txt"))
  
  table_dxy_homozygotes$mid <- (table_dxy_homozygotes$window_pos_1 + table_dxy_homozygotes$window_pos_2) / 2
  
  
  tmp_dxy_homo_1<-table_dxy_homozygotes %>% 
    filter(no_sites>=5000) %>%
    filter(pop1=="Baltic_Spring" & pop2=="WhiteSea") %>%
    filter(chromosome==chr)
  
  pcent_tmp_dxy_homo_1 <- quantile(tmp_dxy_homo_1$avg_dxy, 0.05)
  
  tmp_dxy_homo_2<-table_dxy_homozygotes %>% 
    filter(no_sites>=5000) %>%
    filter(pop2=="Baltic_Autumn" & pop1=="WhiteSea") %>%
    filter(chromosome==chr)
  
  pcent_tmp_dxy_homo_2 <- quantile(tmp_dxy_homo_2$avg_dxy, 0.05)
  
  tmp1<-tmp_dxy_homo_1 %>% 
    filter(window_pos_1 >= s & window_pos_2 <= e) %>% summarise(dxy_BS_vs_WS=mean(avg_dxy))
  
  tmp2<-tmp_dxy_homo_2 %>% 
    filter(window_pos_1 >= s & window_pos_2 <= e) %>% summarise(dxy_BA_vs_WS=mean(avg_dxy))
  
  tmp3<-xpehh_results %>%
    filter(CHROM== chr & POS >= s & POS <=e) %>% summarise(XPEHH_sum=mean(XPEHH), std_XPEHH_mean=mean(std_XPEHH))
  
  # Very Good
  classification1 <- (tmp1 < pcent_tmp_dxy_homo_1 && tmp2 > pcent_tmp_dxy_homo_2 && tmp3$std_XPEHH_mean < -2)
  # Good
  classification2 <- (tmp1 < pcent_tmp_dxy_homo_1 && tmp2 > pcent_tmp_dxy_homo_2 && !tmp3$std_XPEHH_mean > 1)
  
  classified_regions[i,]<-matrix(data=c(chr, s, e, tmp1, pcent_tmp_dxy_homo_1, tmp2, pcent_tmp_dxy_homo_2, tmp3, classification1, classification2), nrow=1, ncol= 11)
  
}

# Particularly convincing
introgressed_regions_very_good <- classified_regions %>% filter(VeryGood==T )

# chr    start      end dxy_BS_vs_WS dxy_pcent_BS_vs_WS dxy_BA_vs_WS dxy_pcent_BA_vs_WS XPEHH_mean std_XPEHH_mean VeryGood Good
# chr    start      end dxy_BS_vs_WS dxy_pcent_BS_vs_WS dxy_BA_vs_WS dxy_pcent_BA_vs_WS XPEHH_mean std_XPEHH_mean VeryGood Good
# 1  chr10 21380001 21400000 0.0007741627        0.001250665  0.006536459        0.001635233 -0.4069019      -2.274428     TRUE TRUE
# 2  chr10 25100001 25120000 0.0004963587        0.001391529  0.003289750        0.001635233 -0.6436167      -3.154223     TRUE TRUE
# 3  chr10 25120001 25140000 0.0007953807        0.001384884  0.003084939        0.001635233 -0.6186389      -3.093933     TRUE TRUE
# 4  chr10 25140001 25160000 0.0006846187        0.001397731  0.005882148        0.001635233 -0.8504188      -3.943150     TRUE TRUE
# 5  chr10 25160001 25180000 0.0011352964        0.001398000  0.004295774        0.001635233 -0.8957956      -4.060855     TRUE TRUE
# 6  chr10 25180001 25200000 0.0009187868        0.001420946  0.004536037        0.001635233 -0.7257754      -3.525363     TRUE TRUE
# 7  chr10 25200001 25220000 0.0005176850        0.001491551  0.003693793        0.001635233 -0.5398388      -2.812242     TRUE TRUE
# 8  chr10 26320001 26340000 0.0010355904        0.001341587  0.002258501        0.001635233 -0.5648750      -2.858833     TRUE TRUE
# 9  chr16 14580001 14600000 0.0011660070        0.001599871  0.004308363        0.001793913 -0.4105469      -2.191994     TRUE TRUE
# 10 chr16 14600001 14640000 0.0008219703        0.001599871  0.003525828        0.001793913 -0.5284055      -2.638631     TRUE TRUE
# 11 chr16 14640001 14700000 0.0011708535        0.001591224  0.002937160        0.001793913 -0.4764528      -2.507714     TRUE TRUE
# 12 chr16 15320001 15340000 0.0014398202        0.001602516  0.006670790        0.001793913 -0.4126547      -2.209229     TRUE TRUE
# 13 chr16 15340001 15360000 0.0006359993        0.001602516  0.003605828        0.001793913 -0.4521212      -2.351228     TRUE TRUE
# 14 chr19  6360001  6380000 0.0005004860        0.001582253  0.005893937        0.001744877 -0.6191529      -2.986948     TRUE TRUE

# Ok (none of these are almost under selection in the oposite direction)
introgressed_regions_good <- classified_regions %>% filter(Good==T )

# chr    start      end dxy_BS_vs_WS dxy_pcent_BS_vs_WS dxy_BA_vs_WS dxy_pcent_BA_vs_WS  XPEHH_mean std_XPEHH_mean VeryGood Good
# chr    start      end dxy_BS_vs_WS dxy_pcent_BS_vs_WS dxy_BA_vs_WS dxy_pcent_BA_vs_WS  XPEHH_mean std_XPEHH_mean VeryGood Good
# 1   chr1 28960001 28980000 0.0010416903        0.001543446  0.003800076        0.001657025 -0.06438769     -0.7942366    FALSE TRUE
# 2   chr1 28980001 29000000 0.0007264900        0.001588648  0.004029011        0.001657025  0.22238112      0.2911055    FALSE TRUE
# 3   chr1 29000001 29020000 0.0010284191        0.001586924  0.004647001        0.001657025  0.37971208      0.8671535    FALSE TRUE
# 4  chr10 21120001 21140000 0.0003970791        0.001305608  0.003345616        0.001635233 -0.14906983     -1.2345129    FALSE TRUE
# 5  chr10 21380001 21400000 0.0007741627        0.001250665  0.006536459        0.001635233 -0.40690189     -2.2744281     TRUE TRUE
# 6  chr10 21480001 21500000 0.0008827555        0.001287677  0.003517969        0.001635233  0.07223072     -0.3387552    FALSE TRUE
# 7  chr10 21500001 21580000 0.0011545090        0.001287677  0.006285313        0.001635233 -0.11761052     -1.1250222    FALSE TRUE
# 8  chr10 22040001 22060000 0.0009691367        0.001124809  0.004758567        0.001635233 -0.01867937     -0.6713740    FALSE TRUE
# 9  chr10 25100001 25120000 0.0004963587        0.001391529  0.003289750        0.001635233 -0.64361668     -3.1542229     TRUE TRUE
# 10 chr10 25120001 25140000 0.0007953807        0.001384884  0.003084939        0.001635233 -0.61863887     -3.0939334     TRUE TRUE
# 11 chr10 25140001 25160000 0.0006846187        0.001397731  0.005882148        0.001635233 -0.85041877     -3.9431505     TRUE TRUE
# 12 chr10 25160001 25180000 0.0011352964        0.001398000  0.004295774        0.001635233 -0.89579557     -4.0608548     TRUE TRUE
# 13 chr10 25180001 25200000 0.0009187868        0.001420946  0.004536037        0.001635233 -0.72577542     -3.5253628     TRUE TRUE
# 14 chr10 25200001 25220000 0.0005176850        0.001491551  0.003693793        0.001635233 -0.53983884     -2.8122418     TRUE TRUE
# 15 chr10 26320001 26340000 0.0010355904        0.001341587  0.002258501        0.001635233 -0.56487502     -2.8588327     TRUE TRUE
# 16 chr12 16220001 16240000 0.0007587875        0.001576810  0.002901449        0.001884019  0.33500935      0.7389797    FALSE TRUE
# 17 chr16 13980001 14000000 0.0004348548        0.001594672  0.003177866        0.001793913 -0.18769592     -1.2793317    FALSE TRUE
# 18 chr16 14080001 14100000 0.0005350016        0.001408944  0.002254131        0.001793913 -0.14850001     -1.1385289    FALSE TRUE
# 19 chr16 14380001 14420000 0.0009370177        0.001519033  0.003325914        0.001793913 -0.03835588     -0.7567673    FALSE TRUE
# 20 chr16 14420001 14440000 0.0011793215        0.001519033  0.003677915        0.001793913 -0.24579671     -1.6124016    FALSE TRUE
# 21 chr16 14440001 14460000 0.0014411347        0.001519033  0.002788909        0.001793913 -0.19151235     -1.3026898    FALSE TRUE
# 22 chr16 14460001 14500000 0.0010302910        0.001519033  0.002674526        0.001793913 -0.26004011     -1.6032690    FALSE TRUE
# 23 chr16 14520001 14580000 0.0015047778        0.001590416  0.005487745        0.001793913 -0.27366745     -1.7273769    FALSE TRUE
# 24 chr16 14580001 14600000 0.0011660070        0.001599871  0.004308363        0.001793913 -0.41054692     -2.1919936     TRUE TRUE
# 25 chr16 14600001 14640000 0.0008219703        0.001599871  0.003525828        0.001793913 -0.52840549     -2.6386314     TRUE TRUE
# 26 chr16 14640001 14700000 0.0011708535        0.001591224  0.002937160        0.001793913 -0.47645282     -2.5077136     TRUE TRUE
# 27 chr16 14700001 14720000 0.0014353052        0.001596850  0.004029987        0.001793913 -0.18546093     -1.3820543    FALSE TRUE
# 28 chr16 15320001 15340000 0.0014398202        0.001602516  0.006670790        0.001793913 -0.41265475     -2.2092294     TRUE TRUE
# 29 chr16 15340001 15360000 0.0006359993        0.001602516  0.003605828        0.001793913 -0.45212122     -2.3512284     TRUE TRUE
# 30 chr19  6360001  6380000 0.0005004860        0.001582253  0.005893937        0.001744877 -0.61915288     -2.9869484     TRUE TRUE

# These are the results presented in Supplementary Table 7
write.table(introgressed_regions_good, file = "introgressed_regions_good_2025-10-20.txt", col.names = T, row.names = F, quote = F)
write.table(introgressed_regions_very_good, file = "introgressed_regions_very_good_2025-10-20.txt", col.names = T, row.names = F, quote = F)
write.table(classified_regions, file = "classified_introgressed_regions_2025-10-20.txt", sep="\t", col.names = T, row.names = F, quote = F)

# Now, let's make the plots again to visualise these results ####

# Supplementary Figure 13 ####
## COMPARE DXY IN INTROGRESSED REGIONS ####
# How can we show that dxy in all regions of introgression is bellow expected for all comparisons with homozygotes?
# Let's read all the dxy files as a list:

# Based on scan1_v01_baltic_alt_ref_summary_filter2
setwd("~/Documents/Postdoc/Project_Herring/Introgression/pixy/2025-10-18_results_homozygotes_maf5/")

intro_reg<-read.table(header=T,"~/Documents/Postdoc/Project_Herring/Introgression/introgression_scan/2023-11-16_analysis/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")
head(intro_reg)
intro_reg_gr <- GRanges(seqnames=intro_reg$seqnames, IRanges(start=intro_reg$start, end=intro_reg$end))

dxy_files<-list.files(pattern = "out_dxy")
regions <- str_split_fixed(dxy_files, pattern="_", 3)[,c(1,2)]
regions_2 <- paste0(regions[,1], "_",regions[,2])

# read the tables into a list
dxy_tables_list <- lapply(dxy_files, read.table, header=T)

# Let's add a classifier to each table within list:
for (x in 1:length(dxy_tables_list)){
  dxy_tables_list[[x]]$region <-  regions_2[x]
}

library(data.table)
# merge all the tables:
dxy_tables_all <- rbindlist(dxy_tables_list)

# convert to genomic ranges to find overlap:
dxy_tables_all_gr <- GRanges(seqnames=dxy_tables_all$chromosome, IRanges(start=dxy_tables_all$window_pos_1, end=dxy_tables_all$window_pos_2),
                             avg_dxy=dxy_tables_all$avg_dxy)

overlaps<-findOverlaps(intro_reg_gr, dxy_tables_all_gr)

dxy_tables_all$type<-"NotIntrogressed"

dxy_tables_all[data.frame(overlaps)$subjectHits,]$type<-"Introgressed"

# Using the classified regions:
introgressed_regions_good
introgressed_regions_very_good

introgressed_regions_good_gr <- GRanges(seqnames=introgressed_regions_good$chr, IRanges(start=introgressed_regions_good$start, end=introgressed_regions_good$end))
introgressed_regions_very_good_gr <- GRanges(seqnames=introgressed_regions_very_good$chr, IRanges(start=introgressed_regions_very_good$start, end=introgressed_regions_very_good$end))

overlaps_irg<-findOverlaps(introgressed_regions_good_gr, dxy_tables_all_gr)

dxy_tables_all$type_irg<-"No"
dxy_tables_all[data.frame(overlaps_irg)$subjectHits,]$type_irg<-"Yes"

overlaps_irvg<-findOverlaps(introgressed_regions_very_good_gr, dxy_tables_all_gr)
dxy_tables_all$type_irvg<-"No"
dxy_tables_all[data.frame(overlaps_irvg)$subjectHits,]$type_irvg<-"Yes"

# Remove windows with low number of sites:
dxy_tables_all_filter5000 <- dxy_tables_all %>% filter(no_sites >=5000)

# Add comparison classification:
dxy_tables_all_filter5000$xy <- paste0(dxy_tables_all_filter5000$pop1, "_", dxy_tables_all_filter5000$pop2)

pcent_tmp_dxy_homo_1 <- quantile(tmp_dxy_homo_1$avg_dxy, 0.05)


dxy_plot_introgressed_regions <- ggplot()+
  geom_boxplot(data=dxy_tables_all_filter5000, aes(y=avg_dxy, x=xy, color=type)) +
  theme_classic()+
  ylab(label = expression(paste("divergence ", (d[xy]))))+ 
  scale_x_discrete(label=c("Baltic spring \nvs Baltic autumn", "Baltic spring \nvs White Sea", "Baltic autumn \nvs White Sea"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text( size=6, color="black"),
        legend.position="none")

ggsave(dxy_plot_introgressed_regions, filename="figures/dxy_plot_introgressed_regions_vs_background_three_populations_2025-10-20.pdf", 
       units="mm", height=90, width= 90)

ggplot()+
  geom_boxplot(data=dxy_tables_all_filter5000, aes(y=avg_dxy, x=xy, color=type_irg)) +
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title.y = element_text(size=6, color="black"))


# Selection score within introgressed regions: 
xpehh_results_gr <- GRanges(seqnames=xpehh_results$CHROM, IRanges(start=xpehh_results$POS, end=xpehh_results$POS)) 
overlaps_selection<-findOverlaps(intro_reg_gr, xpehh_results_gr)

xpehh_results$type<-"NoIntrogressed"
xpehh_results[data.frame(overlaps_selection)$subjectHits,]$type<-"Introgressed"

mean_selection <- ggplot()+
  geom_density(data=xpehh_results, aes(x=XPEHH, color=type)) +
  theme_classic()+
  xlab("Cross population Extended Haplotype Homozigosity (xpEHH)")+
  theme(axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title = element_text(size=6, color="black"),
        legend.position = "none")

std_selection <- ggplot()+
  geom_density(data=xpehh_results, aes(x=std_XPEHH, color=type)) +
  theme_classic()+
  xlab("Standard Deviation xpEHH Baltic spring vs Baltic autumn")+
  geom_vline(xintercept = -2, linetype="dashed", color="gray")+
  geom_vline(xintercept = 2, linetype="dashed", color="gray")+ 
  theme(axis.text.x = element_text(size=6, color="black"),
        axis.text.y = element_text(size=6, color="black"),
        axis.title = element_text(size=6, color="black"),
        legend.position = "none")

all_plots <- grid.arrange(dxy_plot_introgressed_regions , std_selection)
## suplementary figure 4 ####
ggsave(all_plots, filename="figures/dxy_and_selection_comparison_three_pop_comparison_2025-10-20.pdf", units="mm", height=90, width=90)


