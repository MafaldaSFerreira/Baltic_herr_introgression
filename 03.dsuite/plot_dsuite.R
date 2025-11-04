# Code for plotting D-statistics and f4-ratio results from Dsuite Dtrios output
# Results in Figure 1C of the Baltic herring introgression paper

# Libraries ####
library(tidyverse) 

setwd("Manuscript/Figshare/")

# Sprat as outgroup maf 0.05 JK1000 ####
# FOR PUBLICATION ####
D_BBAA <- read.table("3.dsuite/results_dtrios/all_vs_all_clusters_v03_maf5_JK1000_run_2209_2023-11-24_BBAA.txt", header=T, as.is=T)

D_BBAA$p.adjust <- p.adjust(D_BBAA$p.value, method="BH")
D_BBAA_adjust <- D_BBAA %>% filter(p.adjust < 0.05)

plot_order<-c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK","BarentSea","WhiteSea","Balsfjord","Japan","Vancouver")
D_BBAA_adjust$P2<- factor(D_BBAA_adjust$P2, levels=plot_order)
D_BBAA_adjust$P3<- factor(D_BBAA_adjust$P3, levels=plot_order)

D_stat <- D_BBAA_adjust %>% filter(P3 %in% c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")) %>%
  filter(P2 %in% c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK")) %>%
  group_by(P2, P3) %>% 
  select(P1, P2, P3, Dstatistic, p.adjust) %>%
  #summarise(D=max(Dstatistic)) %>%
  ggplot()+
  geom_point(aes(x=P2, y=Dstatistic, color=p.adjust < 0.0001672241), size=0.5) +
  scale_colour_manual(values=c("black", "orange")) +
  #ylim(0,0.02)+
  theme_bw()+
  theme(axis.text.x=element_blank(), 
        legend.position = "none",
        axis.title.y=element_text(size = 5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size = 5),
        strip.text.x = element_text(size = 5),
        panel.spacing = unit(0,'lines'))+
  facet_grid(cols = vars(P3))


Fstat <- D_BBAA_adjust %>% filter(P3 %in% c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")) %>%
  filter(P2 %in% c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK")) %>%
  group_by(P2, P3) %>% 
  select(P1, P2, P3, f4.ratio, p.adjust) %>%
  #summarise(D=max(Dstatistic)) %>%
  ggplot()+
  geom_point(aes(x=P2, y=f4.ratio, color=p.adjust < 0.0001672241), size=0.5) +
  scale_colour_manual(values=c("black", "orange")) +
  #ylim(0,0.02)+
  theme_bw()+
  theme(axis.title.y =element_text(size = 5),
        axis.text.y =element_text(size = 5),
        axis.title.x=element_text(size = 5),
        axis.text.x=element_text(angle=90, size=5), 
        legend.position = "none",
        strip.text.x = element_blank(),
        panel.spacing = unit(0,'lines'))+
  facet_grid(cols = vars(P3)) 

#all_plots <- gridExtra::grid.arrange(D_stat, Fstat, nrow=2)

all_plots<-plot_grid(D_stat, Fstat, nrow=2, align = "h" )

ggsave(all_plots, filename = "figures_dtrios/Dstat_f4.ratio_maf5_JK1000_run_2209_2023-11-24.pdf", units = "mm", width = 78, height = 80)

