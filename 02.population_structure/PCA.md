# Principle Component Analysis



Code written by Sabine Felkel

~~~bash
#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J plink_PCA
#SBATCH -e plink_PCA_%J_%A_%a.err
#SBATCH -o plink_PCA_%J_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabine.felkel@imbim.uu.se


#print script to log
cat script_plink_PCA.sh;


#load modules
module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages/3.6.0 eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;

pwd;
#/proj/snic2021-2-4/private/herring/users/sabine/2023-07-13_introgression/PCA

#01.09.2023

######
# PCA

#get full set of samplenames in order as they occur in vcf
gunzip -cd  /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz| grep -m 1 "CHROM" | cut -f10- | tr '\t' '\n' > samplenames_125vcf.args;
#samplenames of only Atlantic herring (without Pacific herring)
grep -v "White" samplenames_125vcf.args  | grep -v "Pecho" | grep -v "Japan" | grep -v "Bals" | grep -v "Pacif" > samplenames_95woPacvcf.args;

#downsample the vcf to 1 var per 10kbp
python script_2023-08-07_filter_vcf_1SNPeveryXkb.py /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz 10000;
#67566 variants left for PCA
#downsample the vcf to selected variants once with and once without inversion variants
for sel in Top50angela Top50angela_woInv;
do
bcftools view /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz -R ../selected_"$sel".bed -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5."$sel".vcf.gz;
done;
#15231 and 11492 variants left for PCA

#the same two files for without Pacific herring:
gatk IndexFeatureFile --input /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz;
gatk SelectVariants --variant /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz --sample-name samplenames_95woPacvcf.args --exclude-non-variants --max-nocall-fraction 0.2 --output herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.gz;
#3225850 variants
#downsample the vcf to 1 var per 10kbp
python script_2023-08-07_filter_vcf_1SNPeveryXkb.py herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.gz 10000;
#65957 variants left for PCA
#downsample the vcf to selected variants once with and once without inversion variants
for sel in Top50angela Top50angela_woInv;
do
bcftools view herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.gz -R ../selected_"$sel".bed -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac."$sel".vcf.gz;
done;
#15019 and 11302 variants left for PCA

#run PCA with plink
for vcf in *vcf.gz;
do
plink --vcf "$vcf" --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca \
--out "$(basename "$vcf" .vcf.gz)";
done;

#plot PCA in R
for file in *.eigenvec; do sed "s/INFILE/"$(basename "$file" .eigenvec)"/g" script_2023-08-25_plink_PCA_visualise.R | Rscript -; done;

# PCA
######
~~~