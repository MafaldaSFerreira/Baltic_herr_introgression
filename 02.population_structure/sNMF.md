# Ancestry analysis with sNMF



Code written by Sabine Felkel
~~~bash
#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J sNMF_admix
#SBATCH -e sNMF_admix_%J_%A_%a.err
#SBATCH -o sNMF_admix_%J_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabine.felkel@imbim.uu.se


#print script to log
cat script_sNMF_admixture.sh;


#load modules
module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages/3.6.0 eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;

pwd;
#/proj/snic2021-2-4/private/herring/users/sabine/2023-07-13_introgression/sNMF_new

#01.09.2023

############
# admixture

for vcf in ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela.vcf ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.Top50angela_woInv.vcf ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.95AtlwoPac.vcf.1SNPevery10000.vcf ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela.vcf ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv.vcf ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.1SNPevery10000.vcf;
do 
#bgzip "$vcf";
#tabix "$vcf".gz;
bcftools query -f '[ %GT]\n' "$vcf".gz | sed "s/^ //g; s/|/\//g" | sed "s/\.\/\./9/g; s/0\/0/0/g; s/1\/1/2/g; s/0\/1/1/g" | sed "s/ //g" > "$(basename "$vcf" .vcf).geno";
for k in $(seq 2 10);
do for seed in 10 20 30 40 50;
do
#log output not tested...manually created
../../software/sNMF_CL_v1.2/bin/sNMF -x "$(basename "$vcf" .vcf).geno" -K "$k" -s "$seed" -p 1 -m 2 -q "$(basename "$vcf" .vcf).geno.${k}.${seed}.Q" -g "$(basename "$vcf" .vcf).geno.${k}.${seed}.G" > "$(basename "$vcf" .vcf).geno.${k}.${seed}.LeastSquare";
#calculate cross-entropy criteria to get best k - all in one solution:
../../software/sNMF_CL_v1.2/bin/sNMF -x "$(basename "$vcf" .vcf).geno" -K "$k" -s "$seed" -p 1 -m 2 -q "$(basename "$vcf" .vcf).geno.${k}.${seed}.Q" -g "$(basename "$vcf" .vcf).geno.${k}.${seed}.G" -c > "$(basename "$vcf" .vcf).geno.${k}.${seed}.log";
done;
done;
done;

#define best k
for vcf in 95AtlwoPac.Top50angela 95AtlwoPac.Top50angela_woInv 95AtlwoPac.vcf.1SNPevery10000 Top50angela Top50angela_woInv vcf.1SNPevery10000;
do for k in $(seq 2 10);
do
#or use clumpak and put in the least squares: http://clumpak.tau.ac.il/bestK.html (Log Probability table file)
grep "Least-square error" herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5."$vcf".geno.${k}.*.LeastSquare | cut -d " " -f3 | sed "s/^/$k /g" >> herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5."$vcf"_LeastSquare.log;
for seed in 10 20 30 40 50;
do
#the lowest entropy is the best k (or the k where entropy does not any further decrease)
grep "Cross-Entropy (masked data)" herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5."$vcf".geno.${k}.${seed}.log | sed "s/^/$k /g" | sed "s/ /\t/g" | cut -f1,6 | awk -v seed="$seed" '($3="run"seed)' >>  herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5."$vcf"_CrossEntropy.log;
done;
done;
done;

#saved clumpak results to sNMF
mkdir plots;
ls plots/Best* *CLUMP*;

#plot cross entropy
for INFILE in *_CrossEntropy.log; 
do sed "s/INFILE/$INFILE/g" script_2023-08-30_sNMF_plot_CrossEntropy.R | Rscript -; 
done;

#get metadata info for samples in admixture for plotting in the following format for R script
#name_in_vcf	id	sampling_location	region	longitude	latitude	collection_season	species	order_in_vcf	order_in_admixture
#AAL1_CelticSea_Atlantic_Winter	AAL1_CelticSea_Atlantic_Winter	Celtic_Sea	Northeast_Atlantic_Ocean	51.330612	-8.909912	Winter	C.harengus	86
cp ../sNMF/herring_sentieon_125ind_230722_filter_setGT_miss0.2.vcf.minDP4.0maxDP3.0avg_miss0.2_biallelic_MAF5.1SNPevery10000.info sampleinfo125.txt;
cp ../sNMF/herring_sentieon_125ind_230722_filter_setGT_miss0.2.vcf.minDP4.0maxDP3.0avg_miss0.2_biallelic_MAF5.1SNPevery10000_woCpallasii.info sampleinfo95.txt;
#reorder such that it fits to new vcf
for file in sampleinfo*; 
do head -1 "$file" > "$file".ordered; 
gunzip -cd ../PCA_new/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.Top50angela_woInv.vcf.gz | grep "CHROM" | cut -f10- | tr '\t' '\n' | while read in; 
do grep -P "\t$in\t" "$file" >> "$file".ordered; 
done; 
done;

#use customised order for data .... reorder the Q files according to wanted order (mafalda file.. that also goes into the R script afterwards)
tail -n +2 sampleinfo125.txt.ordered | cut -f2 > sampleinfo125.txt.ordered.names;
for file in *maf5.[v,T]*.Q;
do
paste sampleinfo125.txt.ordered.names "$file" | sed "s/ /\t/g" > "$file".temp;
tail -n +2 sampleinfo125.txt.ordered.mafalda | cut -f2 | while read in;
	do
		grep -P "^$in\t" "$file".temp >> "$file".temp2;
	done;
mv "$file".temp2 "$file";
rm "$file".temp;
done;

tail -n +2 sampleinfo95.txt.ordered | cut -f2 > sampleinfo95.txt.ordered.names;
for file in *maf5.9*.Q;
do
paste sampleinfo95.txt.ordered.names "$file" | sed "s/ /\t/g" > "$file".temp;
tail -n +2 sampleinfo95.txt.ordered.mafalda | cut -f2 | while read in;
        do
          	grep -P "^$in\t" "$file".temp >> "$file".temp2;
        done;
mv "$file".temp2 "$file";
rm "$file".temp;
done;


#plot admixtures - if redo, submit both.. it is a quick and dirty solution..
for file in *maf5.[v,T]*.Q; do k_temp="$(rev <<< "$file" | cut -d "." -f3 | rev)"; k=$(( $k_temp+1 )); sed "s/\.2\.10\.Q/\.${k}\.10\.Q/g; s/V5/V${k}/g; s/(5)/(${k})/g; s/INFILE/$file/g" script_2023-08-30_sNMF_plot_Admixtures.R | Rscript -; done;
for file in *woPac*.Q; do k_temp="$(rev <<< "$file" | cut -d "." -f3 | rev)"; k=$(( $k_temp+1 )); sed "s/\.2\.10\.Q/\.${k}\.10\.Q/g; s/V5/V${k}/g; s/(5)/(${k})/g; s/INFILE/$file/g" script_2023-08-30_sNMF_plot_Admixtures_woPac.R | Rscript -; done;



# admixture
############
~~~

