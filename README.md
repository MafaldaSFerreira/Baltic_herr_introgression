# Adaptive introgression from Pacific herring to Atlantic herring in the brackish Baltic Sea

Here I document the code used for analysis of genomic data that was used to study the role of adaptive introgression from Pacific herring to Atlantic herring in adaptation to the Baltic sea environment. 

The work is now available as a pre-print in [BioRxiv](https://www.biorxiv.org/content/10.1101/2025.11.02.685776v1).

The code is divided according to the steps of analysis, as follows:

1. Mapping and Genotype calling
    - i. Atlantic and Pacific herring high coverage whole genome data
    - ii. European sprat long-read sequencing data
2. Population structure analysis
    - i. PCA
    - ii. Admixture (sNMF)
3. Hybridization (Dsuite)
4. Introgression scan
5. Processing Pool-sequencing data
6. Topology weighting (TWISST)
7. Genetic divergence scans (pixy)
8. Positive selection (xpEHH)
9. Classifying introgression regions using dxy and xpEHH
10. Gene Ontology Enrichment analysis of genes in introgressed regions
11. Signatures of acelerated evolution (codeml in PAML)
12. Plotting Figure 3 SEC16B
13. Plotting Figure 4 THRB


All R code links to a FigShare repository containing all the input data necessary to re-plot the figures. At the moment, this FigShare data is only available to reviewers. Upon publication, the link to the repository will be made available here. 

Figures and Tables mentioned in the code are numbered according to the BioRxiv publication.


#### Contact

Mafalda Sousa Ferreira: sferreira.mafalda [at] gmail [.] com