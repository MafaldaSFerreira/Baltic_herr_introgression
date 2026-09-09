# Classifying Introgression Regions based on dxy and positive selection results
# These results are displayed in Supplementary Figure 13
# And in Supplementary Table 7 

# We start with dxy results from pixy which have been generated for group of Baltic spring-spawning homozygotes and White Sea. This is documented in step7.pixy/plot_dxy.R
# We also use the xpEHH results as documented in step 8.

library(tidyverse)
library(GenomicRanges)

# Input files can be found in Figshare:
setwd("~/Documents/Postdoc/Project_Herring/Introgression/Manuscript/Figshare/7.pixy/results_homozygotes_dxy/")

# First let's read the dxy files
dxy_files<-list.files(pattern = "out_dxy")
regions <- str_split_fixed(dxy_files, pattern="_", 3)[,c(1,2)]

# Great, now, we probably need to read in the selection files.
xpehh_results<-read.table("../../8.xpehh/results_xpehh/xpehh_group1+2.out", skip=1)
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
  table_dxy_homozygotes <-read.table(header=T, paste0(chr,"_",s,"-",e,"_Baltic_Spring_WhiteSea_BalticAutumn_popfile.wg.maf5.20kb.popgenpixy.out_dxy.txt"))
  
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
# Results are deposited on GitHub
write.table(introgressed_regions_good, file = "introgressed_regions_good_2025-10-20.txt", col.names = T, row.names = F, quote = F)
write.table(introgressed_regions_very_good, file = "introgressed_regions_very_good_2025-10-20.txt", col.names = T, row.names = F, quote = F)
write.table(classified_regions, file = "classified_introgressed_regions_2025-10-20.txt", sep="\t", col.names = T, row.names = F, quote = F)

# Now, let's make the plots again to visualise these results ####

# Supplementary Figure 13 ####
## COMPARE DXY IN INTROGRESSED REGIONS ####
# How can we show that dxy in all regions of introgression is bellow expected for all comparisons with homozygotes?
# Let's read all the dxy files as a list:

# Based on scan1_v01_baltic_alt_ref_summary_filter2
intro_reg<-read.table(header=T,"../../4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")
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

# 2026-08-05 #### 
# Reviews Comment 6

# Confirm the wg file has the correct dxy values, so we can extract Canada Spring dxy from it.
wg <- read.table("../results_WGS_dxy/2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_dxy.txt", header = TRUE)

wg %>% distinct("Canada_Spring", "WhiteSea")   # confirm Canada_Spring x WhiteSea exists, and the ordering

wg_ba <- wg %>%
  filter((pop1 == "Baltic_Autumn" & pop2 == "WhiteSea") |
           (pop2 == "Baltic_Autumn" & pop1 == "WhiteSea"),
         chromosome == "chr19", window_pos_1 >= 6000000, window_pos_2 <= 6800000)

region_file <- read.table(header = TRUE,
                          "chr19_6360001-6380000_Baltic_Spring_WhiteSea_BalticAutumn_popfile.wg.maf5.20kb.popgenpixy.out_dxy.txt")

region_ba <- region_file %>%
  filter((pop1 == "Baltic_Autumn" & pop2 == "WhiteSea") |
           (pop2 == "Baltic_Autumn" & pop1 == "WhiteSea"),
         chromosome == "chr19")

comparison <- inner_join(
  wg_ba,      # Baltic_Autumn vs WhiteSea, chr19, from the clustered WG run
  region_ba,  # same pair, same chromosome, from the per-region run
  by = c("chromosome", "window_pos_1", "window_pos_2"),
  suffix = c("_wg", "_reg")
)

comparison %>%
  summarise(n         = n(),
            sites_cor = cor(no_sites_wg, no_sites_reg),
            dxy_cor   = cor(avg_dxy_wg,  avg_dxy_reg),
            max_diff  = max(abs(avg_dxy_wg - avg_dxy_reg)))

#all good. so let+s extract Canada Spring
wg_cs <- wg %>%
  filter(pop1 == "Canada_Spring" & pop2 == "WhiteSea") %>% filter()

wg_cs_gr <- GRanges(wg_cs$chromosome,
                    IRanges(wg_cs$window_pos_1, wg_cs$window_pos_2))

vg_gr <- GRanges(classified_regions$chr,
                 IRanges(as.numeric(classified_regions$start),
                         as.numeric(classified_regions$end)))

wg_cs$type <- "Background"
wg_cs$type[subjectHits(findOverlaps(vg_gr, wg_cs_gr))] <- "Introgressed"

wg_cs %>% filter(type=="Introgressed")

wg_cs$introgressed_chr <- classified_regions[queryHits(findOverlaps(vg_gr, wg_cs_gr)),]$chr

wg_cs %>%
  group_by(type) %>%
  summarise(n = n(), mean_dxy = mean(avg_dxy))

left_join(classified_regions, wg_cs, by=join_by("chr"=="chromosome", "start"=="window_pos_1", "end"=="window_pos_2"))

wg %>% filter(chromosome=="chr1" & window_pos_1==28920001 & pop1 == "Canada_Spring" & pop2 == "WhiteSea" & no_sites>=5000)

wg %>%
  filter(no_sites >= 5000 & pop1 == "Canada_Spring" & pop2 == "WhiteSea") %>% summarise(mean(avg_dxy))

classified_regions_copy <- classified_regions

classified_regions_copy$dxy_CS_vs_WS <- NA
classified_regions_copy$pcent_dxy_CS_vs_WS <- NA
 

for(i in 1:nrow(regions)){

  chr<-regions[i,1]
  s<-as.numeric(str_split_fixed(regions[i,2], pattern="-", 2)[1,1])
  e<-as.numeric(str_split_fixed(regions[i,2], pattern="-", 2)[1,2])
  c<-str_remove(chr, "chr")
  
  print(c(chr, s, e))
  
 # Add dxy from Canada
  tmp_dxy_3 <- wg %>%
    filter(no_sites >= 5000 & pop1 == "Canada_Spring" & pop2 == "WhiteSea")
  
  pcent_tmp_dxy_3 <- quantile(tmp_dxy_3$avg_dxy, 0.05)
  
  tmp4 <- tmp_dxy_3 %>%
    filter(chromosome== chr & window_pos_1 >= s & window_pos_2 <= e) %>%
    summarise(dxy_CS_vs_WS = mean(avg_dxy))
  
  classified_regions_copy[i, c(12,13)] <- matrix(data=c(tmp4, pcent_tmp_dxy_3), nrow=1, ncol=2)
  
}

names(classified_regions_copy)

classified_regions_copy <- classified_regions_copy[,c(1:7,12,13,8:11)]

write.table(classified_regions_copy, file = "classified_introgressed_regions_2026-08-05.txt", sep="\t", col.names = T, row.names = F, quote = F)

# Let's try to make a figure to display this: 

# Helper: symmetric population-pair filter (pixy ordering is inconsistent)
pair_filter <- function(tab, popA, popB) {
  tab %>% filter((pop1 == popA & pop2 == popB) | (pop1 == popB & pop2 == popA))
}

# Helper: mean dxy in focal window(s) and in background windows of same chromosome,
# where background = all windows on that chromosome NOT overlapping ANY introgression region
focal_vs_background <- function(tab, chr, s, e) {
  tab <- tab %>% filter(no_sites >= 5000, chromosome == chr)
  if (nrow(tab) == 0) return(c(focal = NA, background = NA))
  gr <- GRanges(tab$chromosome, IRanges(tab$window_pos_1, tab$window_pos_2))
  hits <- GenomicRanges::findOverlaps(intro_reg_gr, gr)
  tab$is_intro <- FALSE
  tab$is_intro[S4Vectors::subjectHits(hits)] <- TRUE
  focal <- tab %>% filter(window_pos_1 >= s, window_pos_2 <= e)
  c(focal      = mean(focal$avg_dxy),
    background = mean(tab$avg_dxy[!tab$is_intro]))
}

dxy_ratios <- list()

for (i in 1:nrow(regions)) {
  
  chr <- regions[i, 1]
  s   <- as.numeric(str_split_fixed(regions[i, 2], pattern = "-", 2)[1, 1])
  e   <- as.numeric(str_split_fixed(regions[i, 2], pattern = "-", 2)[1, 2])
  
  print(c(chr, s, e))
  
  table_dxy_homozygotes <- read.table(header = TRUE, paste0(
    chr, "_", s, "-", e,
    "_Baltic_Spring_WhiteSea_BalticAutumn_popfile.wg.maf5.20kb.popgenpixy.out_dxy.txt"))
  
  bs <- focal_vs_background(pair_filter(table_dxy_homozygotes, "Baltic_Spring", "WhiteSea"), chr, s, e)
  ba <- focal_vs_background(pair_filter(table_dxy_homozygotes, "Baltic_Autumn", "WhiteSea"), chr, s, e)
  cs <- focal_vs_background(pair_filter(wg, "Canada_Spring", "WhiteSea"), chr, s, e)
  
  dxy_ratios[[i]] <- data.frame(
    chr        = chr,
    start      = s,
    end        = e,
    comparison = c("Baltic spring", "Baltic autumn", "Canada spring"),
    focal      = c(bs["focal"],      ba["focal"],      cs["focal"]),
    background = c(bs["background"], ba["background"], cs["background"]),
    row.names  = NULL
  )
}

dxy_ratios <- bind_rows(dxy_ratios) %>%
  mutate(ratio  = focal / background,
         region = paste0(chr, ":", start),
         comparison = factor(comparison,
                             levels = c("Baltic spring", "Baltic autumn", "Canada spring")))

# Attach the classification so you can plot the 14 "very good" regions
dxy_ratios <- dxy_ratios %>%
  left_join(classified_regions %>%
              mutate(start = as.numeric(start), end = as.numeric(end)) %>%
              dplyr::select(chr, start, end, VeryGood, Good),
            by = c("chr", "start", "end"))

# Summary
dxy_ratios %>%
  filter(VeryGood == TRUE) %>%
  group_by(comparison) %>%
  summarise(n = n(), mean_ratio = mean(ratio, na.rm = TRUE),
            min = min(ratio, na.rm = TRUE), max = max(ratio, na.rm = TRUE))

# Plot: one point per region, dashed line = no change from background
dxy_ratio_plot <- ggplot(dxy_ratios %>% filter(VeryGood == TRUE),
                         aes(x = comparison, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.7) +
  scale_y_log10() +
  ylab(expression(paste(d[xy], " in introgressed region / chromosome background"))) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 7, colour = "black"),
        axis.text.y  = element_text(size = 7, colour = "black"),
        axis.title.y = element_text(size = 7, colour = "black"))

ggsave(dxy_ratio_plot,
       filename = "~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Introgression/Manuscript/Figures/Supplements/dxy_ratio_introgressed_vs_background_three_populations_2026-08-05.pdf",
       units = "mm", height = 60, width = 90)

# Add this to Sup Fig. 17:

wg_cs_plot <- pair_filter(wg, "Canada_Spring", "WhiteSea") %>%
  filter(no_sites >= 5000)

gr <- GRanges(wg_cs_plot$chromosome,
              IRanges(wg_cs_plot$window_pos_1, wg_cs_plot$window_pos_2))
wg_cs_plot$type <- "NotIntrogressed"
wg_cs_plot$type[subjectHits(findOverlaps(intro_reg_gr, gr))] <- "Introgressed"
wg_cs_plot$xy <- "Canada_Spring_WhiteSea"

panel_a_data <- bind_rows(
  dxy_tables_all_filter5000 %>% dplyr::select(avg_dxy, type, xy),
  wg_cs_plot %>% dplyr::select(avg_dxy, type, xy)
) %>%
  mutate(xy = factor(xy, levels = c("Baltic_Spring_Baltic_Autumn",
                                    "Baltic_Spring_WhiteSea",
                                    "WhiteSea_Baltic_Autumn",
                                    "Canada_Spring_WhiteSea")))

dxy_plot_introgressed_regions <- ggplot(panel_a_data) +
  geom_boxplot(aes(y = avg_dxy, x = xy, colour = type), outlier.size = 0.3) +
  theme_classic() +
  ylab(expression(paste("divergence ", (d[xy])))) +
  scale_x_discrete(labels = c("Baltic spring \nvs Baltic autumn",
                              "Baltic spring \nvs White Sea",
                              "Baltic autumn \nvs White Sea",
                              "Canada spring \nvs White Sea")) +
  theme(axis.title.x = element_blank(),
        axis.text  = element_text(size = 6, colour = "black"),
        axis.title.y = element_text(size = 6, colour = "black"),
        legend.position = "none")

# We will add an extra panel with the trees illustrating the relationships here:

library(ape)
library(ggplotify)
library(cowplot)

# Depths from present, root = sprat divergence scaled to 1
#   species split (Atl/Pac, ~3 My)  = 0.30
#   introgression (~19 kya)         = 0.05
#   Baltic founding (~8 kya)        = 0.01
trees <- list(
  "Genomic background" = "(Sprat:1.0,(Pacific:0.30,(Atlantic:0.01,Baltic:0.01):0.29):0.70);",
  "Pacific to Baltic"  = "(Sprat:1.0,(Atlantic:0.30,(Baltic:0.05,Pacific:0.05):0.25):0.70);",
  "Baltic to Pacific"  = "(Sprat:1.0,((Atlantic:0.01,Baltic:0.01):0.01,Pacific:0.02):0.98);"
)

tip_cols <- c(Sprat = "grey40", Atlantic = "#1F78B4",
              Baltic = "#6A3D9A", Pacific = "#D95F02")

#Draw the titles inside the plot region with text() instead:
  

plot_tree_panel <- function() {
  op <- par(ps=6, mfrow = c(1, 3), mar = c(7, 1, 5, 1), xpd = NA)
  on.exit(par(op))
  for (nm in names(trees)) {
    tr <- read.tree(text = trees[[nm]])
    plot.phylo(tr, direction = "downwards", edge.width = 1.5,
               show.tip.label = FALSE,
               tip.color = tip_cols[tr$tip.label])
    
    mtext(nm, side = 3, line = 1, cex = 1)
    
    pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    n  <- ape::Ntip(tr)
    text(x = pp$xx[1:n], y = pp$yy[1:n] -0.015,
         labels = tr$tip.label, srt = 90, adj = c(1, 0.5),
         cex = 1, col = tip_cols[tr$tip.label])
  }
}

panel_a <- as.ggplot(plot_tree_panel)

all_plots <- grid.arrange(panel_a, dxy_plot_introgressed_regions , std_selection)

## suplementary figure 4 ####
ggsave(all_plots, filename="~/Dropbox/Mac/Documents/Postdoc/Project_Herring/Introgression/Manuscript/Figures/Supplements/dxy_and_selection_comparison_three_pop_comparison_2026-08-05.pdf", 
       units="mm", height=120, width=90, pointsize=6)

