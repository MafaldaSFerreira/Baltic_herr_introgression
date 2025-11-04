
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggsci)

setwd("Manuscript/Figshare/2.population_structure/admixture/results")


# run_admixture_plots <- function(indTable_file, base_pattern, out_pdf="admixture_plots.pdf") {
#  
#   # indTable_file = "sampleinfo125.txt.ordered.mafalda"
#   # base_pattern  = "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv"
#   # out_pdf       = "figure/Top50angela_woInv.allPlots.pdf"
#   # 
#   # === 1. Read individual table ===
#   indTable <- read.table(indTable_file, header = TRUE)
#   colnames(indTable) <- c("id1","id2","sampling_location","name_region",
#                           "longitude","latitude","collection_season",
#                           "species","order_in_vcf","order_in_admixture")
#   
#   # === 2. Build regex pattern dynamically ===
#   base_pattern_escaped <- gsub("\\.", "\\\\.", base_pattern)
#   pattern <- paste0("^", base_pattern_escaped, "\\.geno\\.[0-9]+\\.10\\.Q$")
#   
#   # === 3. List files ===
#   files <- list.files(path=".", pattern=pattern, full.names=TRUE)
#   
#   # Extract the K values
#   geno_nums <- as.integer(sub(".*\\.geno\\.([0-9]+)\\.10\\.Q$", "\\1", basename(files)))
#   
#   # Sort files by K
#   ord <- order(geno_nums)
#   files <- files[ord]
#   geno_nums <- geno_nums[ord]
#   
#   # === 4. Read Q files ===
#   data_list <- lapply(files, function(f) read.table(f, header=FALSE))
#   
#   # === 5. Prepare color palette ===
#   population_colors <- c("#10345C","#73888C","#0D79F2","#6DB2FF","#404FBF",
#                          "#D95F02","#DF9A65","#72FEFF","#1B9E77","#6DC7AC",
#                          "#C5ACFF","#7570B3")
#   
#   # === 6. Make plots for each K ===
#   list_plots <- list()
#   
#   for (i in seq_along(data_list)) {
#   #  i=1
#     tbl <- data_list[[i]]
#     no <- geno_nums[i]
#     
#     # Merge with indTable, but preserve tbl order
#     mergedAdmixtureTable <- merge(tbl, indTable, by.x="V1", by.y="id2")
#     ordered <- mergedAdmixtureTable[order(mergedAdmixtureTable$order_in_admixture),]
#     ordered$V1 <- factor(ordered$V1, levels=ordered$V1)
#     
#     # 🔑 Ensure consistent color/legend order based on tbl order
#     ordered$name_region <- factor(ordered$name_region, levels=unique(ordered$name_region))
#     
#     # Build legend from this ordering
#     legend <- ordered %>%
#       ggplot(aes(x=V1, y=1, fill=name_region)) +
#       geom_tile() +
#       scale_fill_manual(values=population_colors) +
#       theme_void() +
#       theme(legend.position="none")
#     
#     plot_1 <- ordered %>%
#       pivot_longer(cols=paste0("V",2:(no+1)), 
#                    names_to="K", values_to="admix_proportion") %>%
#       ggplot(aes(x=V1, y=admix_proportion, fill=K)) +
#       geom_bar(position="stack", stat="identity") +
#       ylab(paste0("Ancestry (K", no, ")")) +
#       scale_fill_observable() +
#       theme_classic() +
#       theme(axis.title.y=element_text(size=6),
#             axis.text.y=element_text(size=6),
#             axis.text.x=element_blank(),
#             axis.title.x=element_blank(),
#             legend.position="none",
#             axis.line.x=element_blank(),
#             axis.ticks.x=element_blank())
#     
#     plot_2 <- cowplot::plot_grid(legend, plot_1, align="v", nrow=2, rel_heights=c(1,6))
#     list_plots[[i]] <- plot_2
#   }
#   
#   # === 7. Arrange all plots together ===
#   all_plots <- gridExtra::grid.arrange(grobs=list_plots, nrow=length(list_plots))
#   
#   # === 8. Save PDF ===
#   ggsave(all_plots, filename=out_pdf, units="mm", width=180, height=180)
#   
#   return(invisible(list(files=files, geno_nums=geno_nums, plots=list_plots)))
# }
# 
# 
# 
# # Top50Angela without inversion
# run_admixture_plots(
#   indTable_file = "sampleinfo125.txt.ordered.mafalda",
#   base_pattern  = "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv",
#   out_pdf       = "figure/Top50angela_woInv.allPlots.pdf"
# )
# 
# # Top50Angela with inversion
# run_admixture_plots(
#   indTable_file = "sampleinfo125.txt.ordered.mafalda",
#   base_pattern  = "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela",
#   out_pdf       = "figure/Top50angela_withInv.allPlots.pdf"
# )
# 
# # One SNP every 10kb
# run_admixture_plots(
#   indTable_file = "sampleinfo125.txt.ordered.mafalda",
#   base_pattern  = "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.1SNPevery10000",
#   out_pdf       = "figure/1SNPevery10kb.allPlots.pdf"
# )
# 




#input file gives the sample info in order as they occur in tbl

indTable=read.table("../../plotting_files/sampleinfo125.txt.ordered.mafalda", header=T, sep="\t", as.is = T)
colnames(indTable)=c("id1","id2","sampling_location","name_region","longitude","latitude","collection_season","species","order_in_vcf","order_in_admixture", "color")
#indTable$order_in_admixture <- 1:nrow(indTable)

#input file is the .Q file
# Define the base pattern of your files (everything before ".geno.")

# Top50Angela without inversion
base_pattern <- "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv"

# Top50Angela with Inversion:
base_pattern <- "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela"

# 1 SNP every 10kb
base_pattern <- "herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.1SNPevery10000"


# Escape dots in the base pattern (regex needs \\.)
base_pattern_escaped <- gsub("\\.", "\\\\.", base_pattern)

# Build the regex pattern dynamically
pattern <- paste0("^", base_pattern_escaped, "\\.geno\\.[0-9]+\\.10\\.Q$")


# 1. Get the list of matching files in your working directory
files <- list.files(
  path = ".",  
  pattern = pattern,
  full.names = TRUE
)

# Extract the K value
geno_nums <- sub(
  pattern = ".*\\.geno\\.([0-9]+)\\.10\\.Q$", 
  replacement = "\\1", 
  x = basename(files)
)

geno_nums <- as.integer(geno_nums)  # convert to numeric


# Check you got the right files
print(files)

# 2. Read them all into a list
# Assuming these are whitespace-separated text files without headers
data_list <- lapply(files, function(f) {
  read.table(f, header = FALSE)
})

# plot

list_plots <- list()

for(i in 1:length(data_list)){

  tbl=data_list[[i]]
  no=geno_nums[i]
  
  mergedAdmixtureTable <- merge(tbl, indTable, by.x="V1", by.y="id2")
  ordered <- mergedAdmixtureTable[order(mergedAdmixtureTable$order_in_admixture),]
  ordered$V1 <- factor(ordered$V1, levels=ordered$V1)
  
  if(no==10){
    
    plot_1 <- ordered %>% 
      pivot_longer(cols=c(V2:paste0("V",no+1)),
                   names_to = "K", values_to = "admix_proportion") %>%
      ggplot(aes(x=V1, y=admix_proportion, fill=K))+
      geom_bar(position="stack", stat="identity")+
      #scale_fill_brewer(palette = "Paired")+
      ylab(paste0("Ancestry (K", no, ")"))+
      scale_fill_observable()+
      theme_classic()+
      theme(axis.title.y=element_text(size=6),
            axis.text.y=element_text(size=6),
            axis.text.x = element_blank(),
            # axis.text.x = 
            #   element_text(size=6, angle=90, hjust=1, vjust=0.5),
            axis.title.x=element_blank(),
            legend.position="none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
     
      #      legend.position = "top",
      #       legend.title = element_blank(), 
      #       legend.text = element_text(size=6),
      #       legend.key.size = unit(2, 'mm'), #change legend key size
      #       legend.key.height = unit(2, 'mm'), #change legend key height
      #       legend.key.width = unit(2, 'mm'))+ 
      # guides(fill = guide_legend(nrow = 1))
    
    # add the populations:
    plot_2 <- cowplot::plot_grid(legend, plot_1, align = "v", nrow=2, rel_heights = c(1,4))
    
    list_plots[[i]] <- plot_2
    
  }else{
    
    plot_1 <- ordered %>% 
      pivot_longer(cols=c(V2:paste0("V",no+1)), 
                   names_to = "K", values_to = "admix_proportion") %>%
      ggplot(aes(x=V1, y=admix_proportion, fill=K))+
      geom_bar(position="stack", stat="identity")+
      #scale_fill_brewer(palette = "Paired")+
      ylab(paste0("Ancestry (K", no, ")"))+
      scale_fill_observable()+
      theme_classic()+
      theme(axis.title.y=element_text(size=6),
            axis.text.y=element_text(size=6),
            axis.text.x = element_blank(),
            axis.title.x=element_blank(),
            legend.position="none",
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank())
      # 
      #       legend.position = "top",
      #       legend.title = element_blank(), 
      #       legend.text = element_text(size=6),
      #       legend.key.size = unit(2, 'mm'), #change legend key size
      #       legend.key.height = unit(2, 'mm'), #change legend key height
      #       legend.key.width = unit(2, 'mm')) + 
      # guides(fill = guide_legend(nrow = 1))
     
    # add the populations: 
    plot_2 <- cowplot::plot_grid(legend, plot_1, align = "v", nrow=2, rel_heights = c(1,4))
    
    list_plots[[i]] <- plot_2
  }
}


all_plots <- gridExtra::grid.arrange(grobs = list_plots[c(2,3,4,5,6,7,8,9,1)], nrow=9)

ggsave(all_plots, filename=paste0("figure/",base_pattern,".allPlots.pdf"),
       units="mm",
       width=210,
       height=270)

# Plot the population legend on top:
population_colors=c("#10345C",
                    "#73888C",
                    "#0D79F2",
                    "#6DB2FF",
                    "#404FBF",
                    "#D95F02",
                    "#DF9A65",
                    "#72FEFF",
                    "#1B9E77",
                    "#6DC7AC",
                    "#C5ACFF",
                    "#7570B3")

ordered$name_region <- factor(ordered$name_region, levels=unique(ordered$name_region))

legend <- ordered %>% 
  pivot_longer(cols=c(V2:paste0("V",no+1)), 
               names_to = "K", values_to = "admix_proportion") %>%
  ggplot(aes(x=V1, y=1, fill=name_region))+
  geom_tile()+
  #scale_fill_brewer(palette = "Paired")+
  ylab(paste0("Ancestry (K", no, ")"))+
  scale_fill_manual(values=population_colors)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x =element_blank(),
        axis.title.x=element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

