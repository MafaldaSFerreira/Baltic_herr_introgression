# Selection scan using xpEHH

Code below writen by Sabine Felkel.

Check Figshare for input and output files of this script. 

The output file used for plotting is called `xpehh_group1+2.out`. 

~~~bash
#!/bin/bash -l
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 0-10:00:00
#SBATCH -J xpehh
#SBATCH -e xpehh_%J_%A_%a.err
#SBATCH -o xpehh_%J_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sabine.felkel@imbim.uu.se


#print script to log
cat script_xpehh_selectionscan.sh;


#load modules
module load bioinfo-tools samtools vcftools bcftools/1.14 picard/2.20.4 GATK/4.1.4.1 BEDTools/2.29.2 R_packages eigensoft/7.2.1 Beagle ANGSD/0.940-stable PCAngsd/1.11 biopython/1.78 pixy PhyML/3.3.20190321 plink2 TreeMix/1.12 SMC++/1.15.4;

#pwd;
#/proj/snic2021-2-4/private/herring/users/sabine/2023-07-13_introgression/xpehh_analysis_new

##############
# xpehh 


#input is two haplotype files - each is for one pop; each line is a variant and the phased HT per individual one after the other coded as 1 0
#second file is snp info with chr - chr_pos - 0 - position

#run it for group 1 (atlantic autumn) vs group 2 (baltic spring) and eventually only for the introgressed region (or extract afterwards if runs quickly on whole genome)
ls group1+2.args;

#phase vcf
ln -s /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz .;
#done by madalfa - old command here as example
#java -jar $BEAGLE_ROOT/beagle.jar gt=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/dsuite/data/herring_sentieon_125ind_230722_filter_setGT_miss0.2.vcf.minDP4.0maxDP3.0avg_miss0.2_biallelic_MAF1.outgroup.nomissSprat.miss0.2.bia.vcf.gz out=herring_sentieon_125ind_230722_filter_setGT_miss0.2.vcf.minDP4.0maxDP3.0avg_miss0.2_biallelic_MAF1.outgroup.nomissSprat.miss0.2.bia.vcf_phased nthreads=8;

#first generate vcf for both populations to be compared
gatk --java-options "-Xmx29g" IndexFeatureFile -I herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz;

gatk --java-options "-Xmx29g" SelectVariants \
    -R /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta \
    -V herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz \
    --exclude-non-variants --sample-name group1+2.args \
    -O input.vcf.gz;

#filter to biallelic sites only
bcftools view input.vcf.gz -o input_group1+2.vcf.gz -O z -M2;

#map file for variant info
bcftools view -H input_group1+2.vcf.gz | cut -f1,2 | awk '{printf "%s\t%s_%s\t%f\t%i\n", $1, $1, $2, $2/2000000, $2}'  > input_group1+2.map;

#then subset both and generate the input files for xpehh
bcftools query -f '[ %GT]\n' input_group1+2.vcf.gz -S group1.args | sed "s/|/ /g; s/ //1" > input_group1.hap;
bcftools query -f '[ %GT]\n' input_group1+2.vcf.gz -S group2.args | sed "s/|/ /g; s/ //1" > input_group2.hap;

#run xpehh - columns are 1) the SNP name, 2) the physical position, 3) the integrated haplotype homozygosity in population 1 4) the integrated haplotype homozygosity in population 2 and 5) the unstandardized XP-EHH score
#works only if I execute it directly in terminal for my data, lol?
/proj/snic2020-2-19/private/herring/users/mafalda/software/hapbin/build/xpehhbin --hapA input_group1.hap --hapB input_group2.hap --map input_group1+2.map -o xpehh_group1+2.out > xpehh_group1+2.log;

#plot for mafaldas introgressed candidate regions
sed "s/_/\t/1" xpehh_group1+2.out | grep -e "chr" -e "Index" |  sed "s/ID/chr\tpos/g; s/iHH\tA1/iHH_A1/g" > xpehh_group1+2_mod.out;
#new regions are already collapsed and added 1MB to flanks
tail -n +2 ../scan1_v01_baltic_alt_ref_intro_regions_cov7_min50kb.collapsed.1Mb.bed | sed "s/^chr//g" | sed "s/^/chr/g" | while read CHR START ENDE; 
do 
head -1 xpehh_group1+2_mod.out | sed "s/std XPEHH/std_XPEHH\tflank/g"  > xpehh_group1+2_mod_"$CHR"_"$START"_"$ENDE".out;
grep -P "\t$CHR\t" xpehh_group1+2_mod.out | awk -v START="$START" -v ENDE="$ENDE" '($3 <= START+1000000)' | sed "s/$/\tflank/g" >> xpehh_group1+2_mod_"$CHR"_"$START"_"$ENDE".out;
grep -P "\t$CHR\t" xpehh_group1+2_mod.out | awk -v START="$START" -v ENDE="$ENDE" '($3 >= START+1000000 && $3 <= ENDE-1000000)' | sed "s/$/\ttarget/g" >> xpehh_group1+2_mod_"$CHR"_"$START"_"$ENDE".out; 
grep -P "\t$CHR\t" xpehh_group1+2_mod.out | awk -v START="$START" -v ENDE="$ENDE" '($3 >= ENDE-1000000)' | sed "s/$/\tflank/g" >> xpehh_group1+2_mod_"$CHR"_"$START"_"$ENDE".out;
sed "s/INPUT/xpehh_group1+2_mod_${CHR}_${START}_${ENDE}.out/g" script_xpehh_selectionscan.R | Rscript -;
done;


# xpehh
#############
~~~