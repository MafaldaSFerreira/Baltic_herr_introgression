### Find baltic individuals introgressed for pacific herring haplotypes ###

library(IRanges)
library(GenomicRanges)
library(stringr)
library(tidyverse)

# load the data ####
setwd("Manuscript/Figshare/")

load(file = "4.introgression_scan/intermediate_files/scan1_v01_baltic_alt_ref_summary_filter2.Rdata")

# load introgressed regions ####

intro_reg<-read.table(header=T, "4.introgression_scan/introgression_regions/scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.txt")

## Conver to list ####
scan1_v01_baltic_alt_ref_summary_filter2_cov7<-scan1_v01_baltic_alt_ref_summary_filter2$cov[scan1_v01_baltic_alt_ref_summary_filter2$cov[, "cov"] > 7,]
names(scan1_v01_baltic_alt_ref_summary_filter2_cov7)[2]<-"seqname"
regions_GRange<-makeGRangesFromDataFrame(scan1_v01_baltic_alt_ref_summary_filter2_cov7)

regions.list<-rep(list(list()), nrow(scan1_v01_baltic_alt_ref_summary_filter2_cov7))
names(regions.list)<-paste0(scan1_v01_baltic_alt_ref_summary_filter2_cov7$seqname,":",scan1_v01_baltic_alt_ref_summary_filter2_cov7$start,"-", scan1_v01_baltic_alt_ref_summary_filter2_cov7$end)

# I want to know which regions the particular individual overlaps in the introgression regions

# First, retrieve the individuals in this pacific group
individuals<-names(scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2$pac_red)

# This loop will populate the list of list aboves with the individuals that support each region:
for(i in 1:length(individuals)){
  
  individual<-individuals[i]
  
  overlap<-data.frame(subsetByOverlaps(regions_GRange, scan1_v01_baltic_v_subarctic_alt_ref_intro_20k_lists_filter2$pac_red[individual]))
  
  for(r in 1:nrow(overlap)){
    window<-paste0(overlap[r,1],":", overlap[r,2], "-", overlap[r,3])
    print(window)
    regions.list[[window]]<-append(regions.list[[window]], individual)
    
  }
}

# Now we need to try and consolidate this:

# For each region, find the names that are duplicated. We can split by Spring, and then find which name is duplicated
# These are the individuals we want to use for each region:

target.individuals<-rep(list(list()), nrow(scan1_v01_baltic_alt_ref_summary_filter2_cov7))
names(target.individuals)<-paste0(scan1_v01_baltic_alt_ref_summary_filter2_cov7$seqname,":",scan1_v01_baltic_alt_ref_summary_filter2_cov7$start,"-", scan1_v01_baltic_alt_ref_summary_filter2_cov7$end)


# Create outputs:
for(region in 1:length(target.individuals)){
 
  # This line finds the names that are duplicated in each region and removes the haplotype code, assigning just the duplicated individual names to the list
  target.individuals[[region]]<-str_remove(unlist(regions.list[[region]])[duplicated(str_split_fixed(unlist(regions.list[[region]]), "Spring", 2)[,1])], "_2")
  
  if(length(target.individuals[[region]])!=0) {
    # Create an output data.frame for pixy:
    df<-data.frame(V1=str_remove(unlist(regions.list[[region]])[duplicated(str_split_fixed(unlist(regions.list[[region]]), "Spring", 2)[,1])], "_2"),
                   V2="Baltic_Spring")
    
    #df_output<-rbind(populations, df)
    
    outputfilename<-str_replace(paste0("scan1_v01_baltic_alt_ref_summary_filter2_cov7/",names(target.individuals)[region], "_Baltic_Spring_popfile.txt"), ":", "_")
    
    print(outputfilename)
    
    write.table(file = outputfilename, x= df, quote=F, col.names = F, row.names = F, sep="\t")
  }
}