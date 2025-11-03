# TWISST

We ran the Twisst analysis using two different outgroups: Pacific herring and European sprat. While we present both results in supplements, the final accepted results overall are with the European sprat (Featured in Figure 2).


### Pacific herring as outgroup
Sabine Felkel's original code for Pacific herring is as follows:

Sabine generated different datasets, but we used `samplenames_100_soft.arg`

~~~bash
#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J twisst
#SBATCH -e twisst_%J_%A_%a.err
#SBATCH -o twisst_%J_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabine.felkel@imbim.uu.se


#print script to log
cat script_twisst.sh;


#load modules
module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages/3.6.0 eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;
cp ../twisst/*py .;
cp ../twisst/plot*R .;

#########
# twisst

#twisst for the classic contrast (A-B...P-O) using 100 snp windows
#the group affiliations discussed with mats (0..pacific_vancouver=OG, P..pacific, B..baltic, A..atlantic) but I modified them to take similar samples (PCA and mainly admixture) and one relaxed set with more admixed individuals

pwd;
#/proj/snic2021-2-4/private/herring/users/sabine/2023-07-13_introgression/twisst

#data needed:

#provide list of desired samplenames
ls groups*.tsv samplenames*.args;

#head -1 groups_*
#==> groups_100_soft.tsv <==
#Pacific16_Vancouver_Pacific_A	O

#==> groups_100.tsv <==
#Pacific16_Vancouver_Pacific_A	O

#==> groups_50.tsv <==
#Pacific16_Vancouver_Pacific_A	O

#head -1 samplenames_*
#==> samplenames_100.args <==
#Pacific16_Vancouver_Pacific

#==> samplenames_100_soft.args <==
#Pacific16_Vancouver_Pacific

#==> samplenames_50.args <==
#Pacific16_Vancouver_Pacific

#provide vcf file, phased
ln -s /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz .;
tabix herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz;

#analysis - go through each dataset
for run in 50 100 100_soft;
do 
wind="$(cut -d "_" -f1 <<< "$run")";
#filter vcf for the desired samples
bcftools view herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz -S samplenames_"$run".args \
	-o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_twisst_"$run".vcf.gz -O z;
#generate input geno and tree file - filter * variants
python parseVCF.py -i herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_twisst_"$run".vcf.gz | grep -v "*" | gzip \
	> herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_"$run".geno.gz;
python phyml_sliding_windows.py -T 1 -g herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_"$run".geno.gz --prefix output.phyml_bionj.w"$run" \
	 -w "$wind" --windType sites --model GTR --optimise n;
#run twisst 
python twisst.py -t output.phyml_bionj.w"$run".trees.gz -w output.weights_"$run".csv.gz --outputTopos topologies_"$run".trees -g A -g B -g O -g P --method complete --groupsFile groups_"$run".tsv;
echo "get first idea of windows with B-P introgression (topology 2):";
gunzip -cd output.weights_"$run".csv.gz | nl | awk  '($3>$4)' | awk '($3>$2)' | tail -n +4;
#the window line number can be extracted from output.phyml_bionj.w"$run".data.tsv to get the actual location on chr
echo -e 'AB-OP_raw\tAO-BP_raw\tAP-BO_raw\ttotal\tAB-OP\tAO-BP\tAP-BO' > output.weights_"$run".csv.ext;
#replace zeros with some tiny epsilon
gunzip -cd output.weights_"$run".csv.gz | sed "s/\t0\t/\t0\.000001\t/g; s/^0\t/0\.000001\t/g; s/\t0$/\t0\.000001/g"|tail -n +4 |awk '($4=$1+$2+$3)' | awk '($5=$1/$4)'| awk '($6=$2/$4)'| awk '($7=$3/$4)' | \
	sed "s/ /\t/g" >> output.weights_"$run".csv.ext;
paste output.phyml_bionj.w"$run".data.tsv output.weights_"$run".csv.ext > output.weights_"$run".csv.ext.Rinput;
#NO NEED#sed "s/weights/weights_${run}/g" script_2023-08-28_twisst_plot_mats.R | Rscript -;
done;

#plot twisst results for the candidate regions provided by mafalda:

#cat ../treemix/scan1_v01_baltic_alt_ref_intro_regions_cov7_min10kb.bed | while read CHR START END;
#new regions:
tail -n +2 ../scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed | awk '($2=$2+1000000)' | awk '($3=$3-1000000)' | sed "s/ /\t/g; s/^/chr/g" > scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed;
tail -n +2 ../scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed | sed "s/^/chr/g" | while read CHR START END;
do
for run in 50 100 100_soft;
do
#add here another column that gives info about candidate or non-candidate region used to get some summary statistics
tail -n +2 output.weights_"$run".csv.ext.Rinput > temp;
head -1 output.weights_"$run".csv.ext.Rinput | tr '\n' ' '  | sed "s/$/\tcategory\n/g" > output.weights_"$run"_mod.csv.ext.Rinput;
bedtools intersect -a temp -b scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed | sed "s/$/\tintrogressed/g" >> output.weights_"$run"_mod.csv.ext.Rinput; 
bedtools intersect -v -a temp -b scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed | sed "s/$/\tnot_introgressed/g" >> output.weights_"$run"_mod.csv.ext.Rinput;
sed "s/weights/weights_${run}_mod/g; s/CHR/$CHR/g; s/START/$START/g; s/END/$END/g" script_2023-08-28_twisst_plot.R | Rscript -;
done;
done;
~~~


### European sprat as outgroup

Generate the input `.geno` file with the European sprat.

Check how we generated this file in `1.mapping_calling_variants`

~~~bash
#!/bin/bash
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J twisst
#SBATCH -e twisst_%A_%a.err            # File to which STDERR will be written
#SBATCH -o twisst_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

# adapting the code that Sabine wrote:
#load modules
module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages/3.6.0 eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;

#########
# twisst

#twisst for the classic contrast (A-B...P-O) using 100 snp windows
#the group affiliations discussed with mats (0..pacific_vancouver=OG, P..pacific, B..baltic, A..atlantic) but I modified them to take similar samples (PCA and mainly admixture) and one relaxed set with more admixed individuals

pwd;
#/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/twisst

#provide vcf file, phased
tabix herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.sprat.vcf.gz;

#analysis - go through each dataset
run="100_soft"
wind="$(cut -d "_" -f1 <<< "$run")";

#filter vcf for the desired samples
bcftools view herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.sprat.vcf.gz -S samplenames_"$run".args \
	-o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_twisst_"$run".vcf.gz -O z;

#generate input geno and tree file - filter * variants

python parseVCF.py -i herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_twisst_"$run".vcf.gz | grep -v "*" | gzip \
	> herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_"$run".geno.gz;

~~~


And now we run phyml and Twisst:

run_all_GW_twisst_Mi10.sh
~~~bash
#!/bin/bash
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J twisst
#SBATCH -e twisst_%A_%a.err            # File to which STDERR will be written
#SBATCH -o twisst_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages/3.6.0 eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;

# Base
BD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/twisst"
# Input
ID="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/twisst/data"
# Output
OD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/twisst/results"
# Scripts
#SW="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/twisst/scripts"


echo "Generate trees with phyml"

python ${BD}/scripts/phyml_sliding_windows.py -T 1 -g ${ID}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased_100_soft.geno.gz -Mi 10 -w 100 --windType sites --model GTR --optimise n --prefix ${OD}/output.wg.phyml.w100.Mi10

echo "Phyml: Done"
echo "Run Twisst"

python ${BD}/scripts/twisst.py -t ${OD}/output.wg.phyml.w100.Mi10.trees.gz -w ${OD}/output.wg.phyml.w100.Mi10.weights.csv.gz --outputTopos ${OD}/output.wg.phyml.w100.Mi10.topos -g A -g B -g O -g P --method complete --groupsFile ${ID}/groups_100_soft.tsv

echo "Twisst: Done"
~~~



