# Introgression Scan 

We start the introgression scan with biallelic genotypes that have been phased.

The code to generate this file is in `1.mapping_calling_variants`

Phased file: `/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz`

~~~bash
python /Library/Frameworks/R.framework/Resources/library/HaploDistScan/vcfToGenotype.py -v herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5
~~~

This will generate the following files:
`herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5_GT.txt`

