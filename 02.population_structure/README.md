# Population structure 


This analysis with the following vcf files:

> Pacific and Atlantic herring: biallelic SNPs, maf 5 and miss 0.20

`herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz`

> Only Atlantic herring: biallelic SNPs, maf 5 and miss 0.20

`herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.gz`


These vcf files were filtered to generate three datasets:
- 1 SNP every 10 kb using `script_2023-08-07_filter_vcf_1SNPeveryKkb.py`
- The top 50 SNPs in selection regions defined in `selected_Top50angela.bed`
- The top 50 SNPs in selection regions and excluding inversion regions defined in `selected_Top50angela_woInv.bed`

Check PCA.md to see how the vcf were generated. They were reused for the sNMF analysis. 

