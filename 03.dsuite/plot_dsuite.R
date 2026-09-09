# Code for plotting D-statistics and f4-ratio results from Dsuite Dtrios output
# Results in Figure 1C of the Baltic herring introgression paper

# Libraries ####
library(tidyverse) 

# Assumes your R working directory is the root of this repo, and the Figshare
# data folder is downloaded as a sibling of this repo (see
# 00.revisions_relabeling/00.population_labels.R for the convention).
setwd("../Figshare/")

# # Sprat as outgroup maf 0.05 JK1000 ####
# # FOR PUBLICATION ####
# D_BBAA <- read.table("3.dsuite/results_dtrios/all_vs_all_clusters_v03_maf5_JK1000_run_2209_2023-11-24_BBAA.txt", header=T, as.is=T)
# 
# D_BBAA$p.adjust <- p.adjust(D_BBAA$p.value, method="BH")
# D_BBAA_adjust <- D_BBAA %>% filter(p.adjust < 0.05)
# 
# plot_order<-c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK","BarentSea","WhiteSea","Balsfjord","Japan","Vancouver")
# D_BBAA_adjust$P2<- factor(D_BBAA_adjust$P2, levels=plot_order)
# D_BBAA_adjust$P3<- factor(D_BBAA_adjust$P3, levels=plot_order)
# 
# D_stat <- D_BBAA_adjust %>% filter(P3 %in% c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")) %>%
#   filter(P2 %in% c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK")) %>%
#   group_by(P2, P3) %>% 
#   select(P1, P2, P3, Dstatistic, p.adjust) %>%
#   #summarise(D=max(Dstatistic)) %>%
#   ggplot()+
#   geom_point(aes(x=P2, y=Dstatistic, color=p.adjust < 0.0001672241), size=0.5) +
#   scale_colour_manual(values=c("black", "orange")) +
#   #ylim(0,0.02)+
#   theme_bw()+
#   theme(axis.text.x=element_blank(), 
#         legend.position = "none",
#         axis.title.y=element_text(size = 5),
#         axis.title.x=element_blank(),
#         axis.text.y=element_text(size = 5),
#         strip.text.x = element_text(size = 5),
#         panel.spacing = unit(0,'lines'))+
#   facet_grid(cols = vars(P3))
# 
# 
# Fstat <- D_BBAA_adjust %>% filter(P3 %in% c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")) %>%
#   filter(P2 %in% c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK")) %>%
#   group_by(P2, P3) %>% 
#   select(P1, P2, P3, f4.ratio, p.adjust) %>%
#   #summarise(D=max(Dstatistic)) %>%
#   ggplot()+
#   geom_point(aes(x=P2, y=f4.ratio, color=p.adjust < 0.0001672241), size=0.5) +
#   scale_colour_manual(values=c("black", "orange")) +
#   #ylim(0,0.02)+
#   theme_bw()+
#   theme(axis.title.y =element_text(size = 5),
#         axis.text.y =element_text(size = 5),
#         axis.title.x=element_text(size = 5),
#         axis.text.x=element_text(angle=90, size=5), 
#         legend.position = "none",
#         strip.text.x = element_blank(),
#         panel.spacing = unit(0,'lines'))+
#   facet_grid(cols = vars(P3)) 
# 
# #all_plots <- gridExtra::grid.arrange(D_stat, Fstat, nrow=2)
# 
# all_plots<-plot_grid(D_stat, Fstat, nrow=2, align = "h" )
# 
# ggsave(all_plots, filename = "figures_dtrios/Dstat_f4.ratio_maf5_JK1000_run_2209_2023-11-24.pdf", units = "mm", width = 78, height = 80)
# 
# # Significance:
# # Use your exact plot logic to split into Orange and Black points
# pop_counts_plot <- pops_F %>%
#   mutate(significance = ifelse(p.adjust < 0.0001672241, "Significant", "Non_significant")) %>%
#   pivot_longer(cols = c(P1, P2), names_to = "Position", values_to = "Population") %>%
#   group_by(Population, significance) %>%
#   tally() %>%
#   ungroup() %>%
#   pivot_wider(names_from = significance, values_from = n, values_fill = 0)
# 
# # Print the table to confirm you have a good mix of numbers
# print(pop_counts_plot)
# 
# # Run the final statistical test
# contingency_matrix <- as.matrix(pop_counts_plot[, -1])
# rownames(contingency_matrix) <- pop_counts_plot$Population
# 
# fisher_test <- fisher.test(contingency_matrix, simulate.p.value = TRUE)
# print(fisher_test)
# 
# # Look at the residuals to prove Baltic_Spring stands out
# chisq.test(contingency_matrix)$residuals
# 
# 
# # 1. Create a clean data frame with a grouping variable
# magnitude_df <- pops_F %>%
#   mutate(Target_Group = ifelse(P2 == "Baltic_Spring", "Baltic_Spring", "Other_Atlantic"))
# 
# # 2. Test if Baltic_Spring has significantly higher f4 ratios than the rest
# wilcox_test <- wilcox.test(f4.ratio ~ Target_Group, data = magnitude_df, alternative = "greater")
# print(wilcox_test)
# 
# # 3. Get the median values to see the difference in size
# magnitude_df %>%
#   group_by(Target_Group) %>%
#   summarise(median_f4 = median(f4.ratio), mean_f4 = mean(f4.ratio))
# 
# # Norwegian Costal:
# # 1. Create a clean data frame with a grouping variable
# magnitude_df <- pops_F %>%
#   mutate(Target_Group = ifelse(P2 == "Norway_Costal", "Norway_Costal", "Other_Atlantic"))
# 
# # 2. Test if Baltic_Spring has significantly higher f4 ratios than the rest
# wilcox_test <- wilcox.test(f4.ratio ~ Target_Group, data = magnitude_df, alternative = "greater")
# print(wilcox_test)
# 
# # 3. Get the median values to see the difference in size
# magnitude_df %>%
#   group_by(Target_Group) %>%
#   summarise(median_f4 = median(f4.ratio), mean_f4 = mean(f4.ratio))


### Revisions 2026-08-16 ####
# Something seems inconsistent about this code:
# 1. I can't remember why I filter by BH level twice. I think it might be because
# dsuite does way more tests that the ones that are relevant to us; i.e., we are only interested
# in tests where P3 is one of the Pacific populations. but still, not sure why that resulted in 
# me the pvalue filter twice. So I will try to filter first for the comparisons I want, and then 
# apply the BH correction.
# 2. the wilcox test is applied to a table that I can't recreate. 


# Let's read the table and apply the bonferroni correction
D_BBAA <- read.table("3.dsuite/results_dtrios/all_vs_all_clusters_v03_maf5_JK1000_run_2209_2023-11-24_BBAA.txt", header=T, as.is=T)

# Filter to biologically relevant tests first, where P3 is a Pacific herring and P1 or P2 are Atlantic herring:
# which is the correct tree
atlantic_pops <- c("Baltic_Autumn", "Baltic_Spring", "NorthSea", "Norway_Costal",
                   "Norway_Pelagic", "Canada_Autumn", "Canada_Spring",
                   "Canada_Summer", "UK")
pacific_pops <- c("Vancouver", "Japan", "WhiteSea", "BarentSea", "Balsfjord")

pops_F <- D_BBAA %>%
  filter(P3 %in% pacific_pops) %>%
  filter(P1 %in% atlantic_pops & P2 %in% atlantic_pops)

# nrow(pops_F) is 180 tests:

# Then apply BH correction only to these tests
pops_F$p.adjust <- p.adjust(pops_F$p.value, method="BH")

# let's not apply this filter, but instead apply it in the table:

plot_order<-c("Baltic_Autumn","Baltic_Spring","NorthSea","Norway_Costal","Norway_Pelagic","Canada_Autumn","Canada_Spring","Canada_Summer","UK")
p3_order<-c("WhiteSea", "BarentSea", "Balsfjord", "Japan", "Vancouver")

pops_F$P2<- factor(pops_F$P2, levels=plot_order)
pops_F$P1<- factor(pops_F$P1, levels=plot_order)
pops_F$P3<- factor(pops_F$P3, levels=p3_order)

# Fix the populations labels consistently:
dsuite_label_lookup <- c(
  "Baltic_Autumn"  = "Baltic Autumn",
  "Baltic_Spring"  = "Baltic Spring",
  "NorthSea"       = "North Sea",
  "Norway_Costal"  = "Norway stationary",
  "Norway_Pelagic" = "Norway migratory",
  "Canada_Autumn"  = "Canada Autumn",
  "Canada_Spring"  = "Canada Spring",
  "Canada_Summer"  = "Canada Summer",
  "UK"             = "Britain and Ireland",
  "BarentSea"      = "Barents Sea",
  "WhiteSea"       = "White Sea",
  "Balsfjord"      = "Balsfjord",
  "Japan"          = "Japan",
  "Vancouver"      = "Vancouver"
)

Dstat_heatmap <- ggplot(pops_F, aes(x = P2, y = P1, fill = Dstatistic)) +
  geom_tile(color = "white") +
  geom_point(data = subset(pops_F, p.adjust < 0.05),
             aes(x = P2, y = P1), color = "white", shape = 8, size = 1,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "D-statistic") +
  scale_x_discrete(labels = dsuite_label_lookup) +
  scale_y_discrete(limits = rev(plot_order), labels = dsuite_label_lookup) +
  facet_wrap(~ P3, nrow = 1, labeller = as_labeller(dsuite_label_lookup)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 7)) +
  labs(x = NULL, y = NULL)

F4_heatmap <- ggplot(pops_F, aes(x = P2, y = P1, fill = f4.ratio)) +
  geom_tile(color = "white") +
  geom_point(data = subset(pops_F, p.adjust < 0.05),
             aes(x = P2, y = P1), color = "white", shape = 8, size = 1,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(name = "f4-ratio") +
  scale_x_discrete(labels = dsuite_label_lookup) +
  scale_y_discrete(limits = rev(plot_order), labels = dsuite_label_lookup) +
  facet_wrap(~ P3, nrow = 1, labeller = as_labeller(dsuite_label_lookup)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(size = 7)) +
  labs(x = NULL, y = NULL)

all_plots<-plot_grid(Dstat_heatmap, F4_heatmap, nrow=2, align = "h" )
all_plots

write.figure
write.table(pops_F, file = "3.dsuite/results_dtrios/pops_F_table_20260817.txt", col.names = T, row.names = F, quote = F)

# ---- Table: f4-ratio against a common reference ---------------------------
# The f4-ratio measures excess Pacific ancestry in P2 RELATIVE TO P1, so the
# values are only comparable across populations if P1 is held fixed.
# Canada_Summer is the natural reference: it has the lowest excess sharing of
# any Atlantic population and so appears as P1 in all 40 of its comparisons,
# giving one value per population per Pacific donor (n = 5 each).

ref_dat <- pops_F %>% filter(P1 == "Canada_Summer")

per_pop <- ref_dat %>%
  group_by(population = as.character(P2)) %>%
  summarise(median_f4 = median(f4.ratio),
            mean_f4   = mean(f4.ratio),
            n         = n(),
            .groups   = "drop") %>%
  mutate(population = dsuite_label_lookup[population]) %>%
  arrange(desc(median_f4))

pooled <- bind_rows(
  ref_dat %>%
    filter(P2 != "Baltic_Spring") %>%
    summarise(population = "Other Atlantic (7 populations)",
              median_f4 = median(f4.ratio), mean_f4 = mean(f4.ratio), n = n()),
  ref_dat %>%
    filter(!P2 %in% c("Baltic_Spring", "Norway_Costal")) %>%
    summarise(population = "Other Atlantic, excl. Norway stationary (6)",
              median_f4 = median(f4.ratio), mean_f4 = mean(f4.ratio), n = n())
)

f4_reference_table <- bind_rows(per_pop, pooled)
print(f4_reference_table)

# # A tibble: 10 × 4
# population                                  median_f4 mean_f4     n
# <chr>                                           <dbl>   <dbl> <int>
#   1 Baltic Spring                                 0.0141  0.0157      5
# 2 Norway stationary                             0.0104  0.0118      5
# 3 Baltic Autumn                                 0.00718 0.00884     5
# 4 Norway migratory                              0.00495 0.00612     5
# 5 Canada Spring                                 0.00490 0.00675     5
# 6 Canada Autumn                                 0.00417 0.00597     5
# 7 North Sea                                     0.00375 0.00512     5
# 8 Britain and Ireland                           0.00321 0.00457     5
# 9 Other Atlantic (7 populations)                0.00495 0.00703    35
# 10 Other Atlantic, excl. Norway stationary (6)   0.00467 0.00623    30

write.table(f4_reference_table,
            "3.dsuite/results_dtrios/f4_reference_table_20260817.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Reviewer comments #6:
# Are there any statistically significant differences in the fraction of introgression among Atlantic herring populations?
# Whenever P1 or P2 is Baltic_spring, how many times is the test significant?
pops_F %>% filter(P1=="Baltic_Spring" | P2=="Baltic_Spring") %>% summarise(n=n())
# This should be 40 tests
pops_F %>% 
  filter(P2=="Baltic_Spring") %>% 
  mutate(win = ifelse(p.adjust <= 0.05, "significant", "non_significant")) %>%
  group_by(P1, P2, win) %>%
  summarise(n=n())
# Baltic_spring is signficiant in 39 out of 40 tests
# A tibble: 9 × 4
# Groups:   P1, P2 [8]
# P1             P2            win                 n
# <fct>          <fct>         <chr>           <int>
#   1 Baltic_Autumn  Baltic_Spring significant         5
# 2 NorthSea       Baltic_Spring significant         5
# 3 Norway_Costal  Baltic_Spring non_significant     1
# 4 Norway_Costal  Baltic_Spring significant         4
# 5 Norway_Pelagic Baltic_Spring significant         5
# 6 Canada_Autumn  Baltic_Spring significant         5
# 7 Canada_Spring  Baltic_Spring significant         5
# 8 Canada_Summer  Baltic_Spring significant         5
# 9 UK             Baltic_Spring significant         5

pops_F %>% 
  filter(P1=="Norway_Costal" | P2=="Norway_Costal") %>% 
  mutate(win = ifelse(p.adjust <= 0.05, "significant", "non_significant")) %>%
  group_by(P1, P2, win) %>%
  summarise(n=n())

# # A tibble: 9 × 4
# # Groups:   P1, P2 [8]
# P1             P2            win                 n
# <fct>          <fct>         <chr>           <int>
#   1 Baltic_Autumn  Norway_Costal significant         5
# 2 NorthSea       Norway_Costal significant         5
# 3 Norway_Costal  Baltic_Spring non_significant     1
# 4 Norway_Costal  Baltic_Spring significant         4
# 5 Norway_Pelagic Norway_Costal significant         5
# 6 Canada_Autumn  Norway_Costal significant         5
# 7 Canada_Spring  Norway_Costal significant         5
# 8 Canada_Summer  Norway_Costal significant         5
# 9 UK             Norway_Costal significant         5
