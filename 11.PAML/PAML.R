setwd("~/Dropbox/Mac (2)/Documents/Postdoc/Project_Herring/SEC16B/PAML/")

library(ape)


#Rainbow trout
onc_mykiss_cds<-read.FASTA("Oncorhynchus_mykiss_sec16b_202_CDS.fa",type="DNA")
onc_mykiss_aa<-read.FASTA("Oncorhynchus_mykiss_sec16b_202_AA.fa",type="AA")

# Zebrafish
dan_rer_cds<-read.FASTA("Danio_rerio_sec16b_201_CDS.fa",type="DNA")
dan_rer_aa<-read.FASTA("Danio_rerio_sec16b_201_AA.fa",type="AA")

# Clupeidae alignment
Align_CDS<-read.FASTA("SEC16B_CDS.fasta",type="DNA")
Align_AA<-read.FASTA("SEC16B_CDS.translated.translated.fasta",type="AA")

# Atlantic herring Baltic spring
baltic_spring_cds<-Align_CDS[7]
baltic_spring_aa<-Align_AA[7]

# Atlantic herring spring 
atlantic_spring_cds <- Align_CDS[2]
atlantic_spring_aa <- Align_AA[2]

  # European sprat
sprat_sprat_cds<-Align_CDS[4]
sprat_sprat_aa<-Align_AA[4]

###### 4 species
# Only in this alignment the aa are all included.

sec16b_align_aa<-clustalomega(c(onc_mykiss_aa,dan_rer_aa,baltic_spring_aa,atlantic_spring_aa,sprat_sprat_aa),
                              exec = "~/Documents/Postdoc/Software/clustalo")


write.FASTA(sec16b_align_aa,file="4_species/sec16b_4species_atlantic.aa.fasta")
write.FASTA(c(onc_mykiss_cds,dan_rer_cds,baltic_spring_cds,atlantic_spring_cds,sprat_sprat_cds),file="4_species/sec16b_4species_atlantic.cds.fasta")

# I converted the AA alignment into codons using PAL2NAL
SEC16B_PAL2NAL <- read.dna("4_species/sec16b_4species_atlantic.cds.algn.fasta", as.matrix = T, format = "fasta")
image(SEC16B_PAL2NAL)
image(trans(SEC16B_PAL2NAL, codonstart = 1))
dnds(SEC16B_PAL2NAL, codonstart = 1)

SEC16B_PAL2NAL_dist<-dist.dna(SEC16B_PAL2NAL)
SEC16B_PAL2NAL_tree<-bionj(SEC16B_PAL2NAL_dist)
SEC16B_PAL2NAL_tree$tip.label <- c("Rainbow_trout", "Zebrafish", "Baltic_spring", "Baltic_autumn", "European_sprat")

tip_color_vec <- rep("black", times = length(SEC16B_PAL2NAL_tree$tip.label))
tip_color_vec[3] <- "dodgerblue1"
tip_color_vec[4] <- "orange"
tree_angle <- 30

plot.phylo(SEC16B_PAL2NAL_tree, tip.color = tip_color_vec, type = "unrooted", rotate.tree = tree_angle , edge.width = 2, main = "Nucleotide distance", lab4ut = "axial")

#Generating files need for PAML analysis
dim(SEC16B_PAL2NAL)
write(x = paste0(dim(SEC16B_PAL2NAL)[1], " ", dim(SEC16B_PAL2NAL)[2]), file = "4_species/SEC16B_BS_ATL_PAML.nuc")
for(i in 1:dim(SEC16B_PAL2NAL)[1]){
  write(x = rownames(SEC16B_PAL2NAL)[i], file = "4_species/SEC16B_BS_ATL_PAML.nuc", append = T)
  write(x = toupper(paste(as.character(SEC16B_PAL2NAL[i,]), collapse = "")), file = "4_species/SEC16B_BS_ATL_PAML.nuc", append = T)
}

#Generating the guide tree
F_ToL <- read.tree("actinopt_12k_raxml.tre")
Guide_tree_SEC16B <- keep.tip(F_ToL, grep("Clupea_harengus|Oncorhynchus_mykiss|Danio_rerio|Sprattus_sprattus", F_ToL$tip.label))
Guide_tree_SEC16B$tip.label <- c( "Rainbow_trout", "Zebrafish", "Herring","European_sprat")
plot.phylo(Guide_tree_SEC16B, type = "unrooted", tip.color = tip_color_vec, rotate.tree = 90, edge.width = 2, main = "RAxML")
write.tree(phy = Guide_tree_SEC16B, file = "4_species/SEC16B_BS_ATL_PAML.tre")
#Added herring sequences and removed branch lengths
Guide_tree_SEC16B_topo <- read.tree("4_species/SEC16B_BS_ATL_PAML_edited.tre")
plot.phylo(Guide_tree_SEC16B_topo, ,  tip.color = tip_color_vec, rotate.tree = 90, edge.width = 2, main = "Guide")

unroot(Guide_tree_SEC16B_topo)->Guide_tree_SEC16B_topo_unroot
plot.phylo(Guide_tree_SEC16B_topo_unroot,  tip.color = tip_color_vec, rotate.tree = 90, edge.width = 2, main = "Guide")
write.tree(Guide_tree_SEC16B_topo_unroot, file="4_species/SEC16B_BS_ATL_PAML_edited_topo_unrooted.tre")


