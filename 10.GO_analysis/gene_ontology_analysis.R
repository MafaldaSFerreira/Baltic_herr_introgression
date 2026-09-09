# R script to conduct gene ontology enrichment analysis with the set of genes that are
# introgressed in Baltic herring


#if (!requireNamespace("biomaRt", quietly = TRUE)) {
#  install.packages("BiocManager")
#  BiocManager::install("biomaRt")
#}

library(biomaRt)
library(GenomicRanges)
library(dplyr)

# Input files for this can be found on Github
setwd("~/Baltic_herr_introgression")

# Ensembl 109 archive (Feb 2023 release, Atlantic herring available)
mart109 <- useMart("ensembl",
                   dataset = "charengus_gene_ensembl",
                   host = "https://feb2023.archive.ensembl.org")


# paste ids from my clipboard:
old_ids <- read.delim(pipe("pbpaste"), header=F)

old_ids <- read.table("10.GO_analysis/list_of_genes_introgressed_regions.txt")

# get attributes:
anno109 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                "chromosome_name", "start_position", "end_position"),
                 filters = "ensembl_gene_id",
                 values = old_ids,
                 mart = mart109)

# Gene ontology enrichment analysis
install.packages("topGO")
BiocManager::install("topGO")

library(topGO)
library(dplyr)

# Suppose go_terms is your BioMart table:
# ensembl_gene_id | go_id

# Query GO terms
go_terms <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "go_id", "name_1006", "namespace_1003"),
  mart = mart109
)

head(go_terms)


# Build gene2GO mapping
gene2GO <- by(go_terms$go_id, go_terms$ensembl_gene_id, function(x) as.character(unique(x)))

# Define gene universe (all genes with GO terms)
allGenes <- names(gene2GO)

# Define your interesting genes (subset)
myGenes <- old_ids$V1
geneList <- factor(as.integer(allGenes %in% myGenes))
names(geneList) <- allGenes


# Check levels
table(geneList)


# Build topGO object
GOdataBP <- new("topGOdata",
              ontology = "BP",                 # BP, MF, or CC
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = gene2GO)
+
# Let's run for all types of GO terms:
# -----------------------------
# Function to run topGO for one ontology
# -----------------------------
run_topGO <- function(geneList, gene2GO, ontology, what_alg, what_stat){
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = geneList,
                annot = annFUN.gene2GO,
                gene2GO = gene2GO)
  
  # Fisher test (elim algorithm)
  result_test <- runTest(GOdata, algorithm = what_alg, statistic = what_stat)
  
  # Extract top 100 terms
  topGO <- GenTable(GOdata, PValue = result_test, orderBy = "PValue", topNodes = 100)
  topGO$PValue <- as.numeric(topGO$PValue)
  
  # Add FDR-adjusted p-value and ontology label
  topGO <- topGO %>%
    mutate(PValue_FDR = p.adjust(PValue, method = "BH"),
           Ontology = ontology)
  
  return(topGO)
}

# -----------------------------
# Run for BP, MF, CC
# -----------------------------

topBP_weight01_Fisher <- run_topGO(geneList, gene2GO, "BP", "weight01", "fisher")
topMF_weight01_Fisher <- run_topGO(geneList, gene2GO, "MF", "weight01", "fisher")
topCC_weight01_Fisher <- run_topGO(geneList, gene2GO, "CC", "weight01", "fisher")

go_all <- rbind(topBP_weight01_Fisher, topMF_weight01_Fisher, topCC_weight01_Fisher)

# add back full gene names:
full_terms <- getBM(attributes=c("go_id", "name_1006"),
                    filters="go",
                    values=go_all$GO.ID,
                    mart=mart109)


go_all_terms_full <- go_all %>%
  left_join(full_terms, by=c("GO.ID"="go_id")) %>%
  mutate(Term_full = name_1006)


# Select top N enriched GO terms
topBP_to_plot_classic <- topBP_classic_Fisher %>%
  filter(!is.na(PValue_FDR)) %>%
  slice_min(PValue_FDR, n = 10)  # top 20 terms

# Select top N enriched GO terms
topBP_to_plot_elim <- topBP_elim_Fisher %>%
  filter(!is.na(PValue_FDR)) %>%
  slice_min(PValue_FDR, n = 10)  # top 20 terms

topBP_to_plot_weight <- topBP_weight01_Fisher %>%
  filter(!is.na(PValue_FDR)) %>%
  slice_min(PValue_FDR, n = 10)  # top 20 terms

# Dotplot
p1 <- ggplot(topBP_to_plot_weight, aes(x = reorder(Term, -PValue_FDR), y = -log10(PValue_FDR), size = Significant)) +
  geom_point(color = "steelblue") +
  coord_flip() +  # horizontal labels
  labs(x = "GO Term",
       y = "-log10(FDR-adjusted p-value)",
       size = "Genes in your list",
       title = "Top Enriched GO Terms") +
  theme_minimal(base_size = 12)

p2 <- ggplot(topBP_to_plot_classic, aes(x = reorder(Term, -PValue_FDR), y = -log10(PValue_FDR), size = Significant)) +
  geom_point(color = "steelblue") +
  coord_flip() +  # horizontal labels
  labs(x = "GO Term",
       y = "-log10(FDR-adjusted p-value)",
       size = "Genes in your list",
       title = "Top Enriched GO Terms") +
  theme_minimal(base_size = 12)

gridExtra::grid.arrange(p1, p2)

# I will go with weight01 and Fisher test:


# Assume go_all is already created (from topBP, topMF, topCC)
# Select top N per ontology
topGo_terms_weight01_Fisher <- go_all_terms_full %>%
  filter(!is.na(PValue_FDR)) %>%
  group_by(Ontology) %>%
  slice_min(PValue_FDR, n = 10) %>%   # top 15 per ontology
  ungroup() %>%
  mutate(logFDR = -log10(PValue_FDR))

topGo_terms_weight01_Fisher <- topGo_terms_weight01_Fisher %>%
  mutate(Term_wrapped = str_wrap(Term_full, width = 50))

# Faceted dotplot
p1 <- ggplot(topGo_terms_weight01_Fisher, aes(x = reorder(Term_wrapped, logFDR),
                        y = logFDR,
                        size = Significant)) +
  geom_point(aes(color = Ontology)) +
  coord_flip() +
  facet_wrap(~Ontology, scales = "free_y", ncol = 1) +  # separate panels
  labs(x = "GO Term",
       y = "-log10(FDR-adjusted p-value)",
       size = "Genes in list",
       color = "Ontology",
       title = "Top Enriched GO Terms per Ontology") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 7, hjust=1))

ggsave(p1, filename="figures/GO_term_enrichment_weight01_fisher.pdf", units = "mm", height = 180, width = 180)
ggsave(p1, filename="figures/GO_term_enrichment_weight01_fisher.png", units = "mm", height = 180, width = 180)


setwd("~/Documents/Postdoc/Project_Herring/Introgression/gene_ontology_analysis/")

# Expor~/Documents/Postdoc/Project_Herring/Introgression/admixtools/# Export all results
write.csv(go_all, "results_topGO/GO_enrichment_all_results_weight01_fisher.csv", row.names = FALSE)


# Suppose topGO_terms is your GenTable output
# gene2GO is a list: names are gene IDs, values are GO IDs assigned to that gene
# geneList is a factor: 1 = interesting gene, 0 = background

# Get only genes in your list
sig_genes <- names(geneList)[geneList == 1]

# Function to get genes per GO term
getGenesPerGO <- function(go_id, gene2GO, sig_genes){
  # For each gene in the significant list, check if it has this GO ID
  genes <- sig_genes[sapply(sig_genes, function(g) go_id %in% gene2GO[[g]])]
  return(paste(genes, collapse = ", "))
}

# Add gene names to topGO_terms
go_all_terms_full_names <- go_all_terms_full %>%
  rowwise() %>%
  mutate(Genes = getGenesPerGO(GO.ID, gene2GO, sig_genes)) %>%
  ungroup()


convertGenesToExternal <- function(genes_string, id_map){
  # Split string into individual Ensembl IDs
  gene_list <- str_split(genes_string, ",\\s*")[[1]]
  
  # Map each Ensembl ID to external name
  ext_names <- id_map$external_gene_name[match(gene_list, id_map$ensembl_gene_id)]
  
  # Replace missing names with original Ensembl IDs
  ext_names[is.na(ext_names)] <- gene_list[is.na(ext_names)]
  
  # Collapse back into a single string
  paste(ext_names, collapse = ", ")
}

# let's annotate with gene names
ensembl_ids <- unique(unlist(gene2GO))
id_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart109
)

# Apply to your topGO_terms
go_all_terms_full_names_gene_names <- go_all_terms_full_names %>%
  rowwise() %>%
  mutate(Genes_external = convertGenesToExternal(Genes, anno109)) %>%
  ungroup() 

write.csv(go_all_terms_full_names_gene_names, "results_topGO/GO_enrichment_all_results_weight01_fisher_with_names.csv", row.names = FALSE)
