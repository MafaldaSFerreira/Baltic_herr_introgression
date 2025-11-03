# Genotype calling

We are starting with previously generated .bam files from Han et al (2020), Pettersson et al (2023), Jamsandekar et al (2024) and Fuentes-Pardo et al (2024). 

Joint genotype calling had not been done for all these data before, so I re-called genotypes using sentieon.

### Haplotyper
~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p devcore
#SBATCH -n 10
#SBATCH -M rackham
#SBATCH --array=1-125:1
#SBATCH -t 24:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J haplotyper
#SBATCH -e haplotyper_%A_%a.err
#SBATCH -o haplotyper_%A_%a.out

conda activate /proj/snic2020-2-19/private/herring/users/mafalda/hifi_mapping

bam_dir=/crex/proj/snic2020-2-19/private/herring/alignment/125_individuals

# Starting with alignments:
BAM=$(ls ${bam_dir}/*.MD.RG.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $BAM

sh /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/scripts/herring_sention_125i.sh ${BAM}
~~~

### Genotyper 

First, we make a list of vcf files:

~~~bash
for GVCF in */*.g.vcf.gz
    do LST=gvcf_count.txt
    ls $GVCF >> gvcf_count.txt
    echo -ne " -v $GVCF" >> gvcf_list.txt
done
~~~


Run `GVCFtyper`. This will also emit variant and invariantn sites with the flag `--emit_mode CONFIDENT`

~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p node
#SBATCH -t 20
#SBATCH -t 24:00:00
#SBATCH -J genotyper
#SBATCH -e genotyper_%A_%a.err            # File to which STDERR will be written
#SBATCH -o genotyper_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com


TOPDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call
OUTPUTDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/joint_call

cd $TOPDIR

REFERENCE=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#sentieon options
# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files_v02/Uppsala_University_node-51.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files_v02/sentieon-genomics-202112.07

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REFERENCE -t 20 --algo GVCFtyper --emit_mode CONFIDENT ${OUTPUTDIR}/herring_sentieon_125ind_231031.vcf.gz $(<gvcf_list.txt)
~~~

Rename individuals in the vcf file to avoid future confusion.

 ~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 24:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J rename_split
#SBATCH -e rename_split_%A_%a.err
#SBATCH -o rename_split_%A_%a.out

ml load bioinfo-tools bcftools/1.17

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/joint_call"
output=${WD}"/chromosomes"

cd ${WD}

bcftools reheader -s old_names_new_names.txt -o herring_sentieon_125ind_231031.newID.vcf.gz herring_sentieon_125ind_231031.vcf.gz

bcftools index herring_sentieon_125ind_231031.newID.vcf.gz

# # Split files by chromosome:
# echo "Run by chromosome"

mkdir ${output}

for c in $(seq 1 26); do
    bcftools view -O z -o ${output}/herring_sentieon_125ind_231031.newID.chr${c}.vcf.gz -r chr${c} ${WD}/herring_sentieon_125ind_231031.newID.vcf.gz

    bcftools index ${output}/herring_sentieon_125ind_231031.newID.chr${c}.vcf.gz
done
~~~

### Genotype filtering

Applyt hard filters

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH --array=1-26:1
#SBATCH -t 10:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J filter1
#SBATCH -e filter2_%A_%a.err
#SBATCH -o filter3_%A_%a.out

# Load GATK
ml load bioinfo-tools GATK/4.3.0.0

# Select Chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/joint_call/chromosomes"
OUTPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes"

gatk IndexFeatureFile --input ${WD}/herring_sentieon_125ind_231031.newID.${ChrName}.vcf.gz

gatk --java-options "-Xmx29g" VariantFiltration  \
    -R /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta \
    -V ${WD}/herring_sentieon_125ind_231031.newID.${ChrName}.vcf.gz \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0))" \
    --filter-name "RPRS8" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('QD') && QD < 2.0))" \
    --filter-name "QD2" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('FS') && FS > 60.0))" \
    --filter-name "FS60" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('SOR') && SOR > 3.0))" \
    --filter-name "SOR3" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQ') && MQ < 40.0))" \
    --filter-name "MQ40" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
    --filter-name "MQ12.5" \
    -G-filter "vc.isSNP() && DP < 3" \
    -G-filter-name "gtDP3" \
    -G-filter "vc.isSNP() && GQ < 20" \
    -G-filter-name "gtGQ10" \
    -O ${OUTPUT_DIR}/herring_sentieon_125ind_231031.newID.filter.${ChrName}.vcf.gz

gatk IndexFeatureFile --input ${OUTPUT_DIR}/herring_sentieon_125ind_231031.newID.filter.${ChrName}.vcf.gz

#set filtered GT to no call

gatk --java-options "-Xmx29g" SelectVariants  \
    -R /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta \
    -V ${OUTPUT_DIR}/herring_sentieon_125ind_231031.newID.filter.${ChrName}.vcf.gz \
    --set-filtered-gt-to-nocall \
    -O ${OUTPUT_DIR}/herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.vcf.gz

~~~

Run depth filter

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 10
#SBATCH -M rackham
#SBATCH --array=1-26:1
#SBATCH -t 10:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J depth
#SBATCH -e depth_%A_%a.err
#SBATCH -o depth_%A_%a.out

## Load 
ml load bioinfo-tools bcftools/1.17

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

## Working directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes"

cd ${WD}

## Do a final filtering, keeping only SNPs that PASS, remove indels, missing ref positions and wildcard ALT
bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.vcf.gz herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.vcf.gz

bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.vcf.gz

bcftools stats -s - herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.vcf.gz > herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.stats

bgzip -cd herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.vcf.gz > herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.vcf

# Depth filter:
python /proj/snic2020-2-19/private/herring/variants/herring_125individuals/script_2023-10-30_beta_vcf_filter_DP.py herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.vcf /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/herring_sentieon_125ind_231031.IndvAvgDepth 3 3

bgzip herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.minDP3.0maxDP3.0avg.vcf
bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.${ChrName}.minDP3.0maxDP3.0avg.vcf.gz
~~~


### Minor allele frequency and missing data filters

We then filtered these vcf files for to only include positions with genotype information at for at least 80% of the individuals.

We also filtered biallelic variants with minor allele frequency < 0.05 or 0.01. In the manuscript, we present results with minor allele frequency > 0.05. (FILTER 3 and FILTER 5)

The code also generates allSites vcf files.

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH --array=1
#SBATCH -t 5:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J depth
#SBATCH -e depth_%A_%a.err
#SBATCH -o depth_%A_%a.out

## Load 
ml load bioinfo-tools bcftools/1.17

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

## Working directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes"

cd ${WD}

## Various filters:

# FILTER 1: miss 0.2
mkdir filter_minmiss20

bcftools view -e "F_MISSING > 0.2" -O z -o ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.vcf.gz

bcftools index ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz > ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.stats

# FILTER 2: miss 0.2 biallelic snps maf 0.01
mkdir filter_minmiss20_biallelic_snps_maf1

bcftools view -m 2 -M 2 -v snps --min-af 0.01:minor -O z -o ${WD}/filter_minmiss20_biallelic_snps_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf1.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools index ${WD}/filter_minmiss20_biallelic_snps_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf1.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20_biallelic_snps_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf1.vcf.gz > ${WD}/filter_minmiss20_biallelic_snps_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf1.stats 

# FILTER 3: miss 0.2 biallelic snps maf 0.05
mkdir filter_minmiss20_biallelic_snps_maf5

bcftools view -m 2 -M 2 -v snps --min-af 0.05:minor -O z -o ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools index ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz ${WD}/filter_minmiss20_biallelic_snps_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.stats

# FILTER 4: miss 0.2 all sites snps maf 0.01
mkdir filter_minmiss20_allsites_maf1

## here we also want multiallelic sites, so we don't set the -M option. The -m will make sure the site has at least two alleles
bcftools view -m 2 --min-af 0.01:minor -O z -o ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.var.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

## to get invariant sites (i.e., ref calls, set --max-af of non-reference alleles to 0)
bcftools view --max-af 0:nref -O z -o ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.inv.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools index ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.var.vcf.gz

bcftools index ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.inv.vcf.gz

bcftools concat --allow-overlaps ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.var.vcf.gz ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.inv.vcf.gz -O z -o ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.allSites.vcf.gz

bcftools index ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.allSites.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.allSites.vcf.gz ${WD}/filter_minmiss20_allsites_maf1/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf1.allSites.stats

# FILTER 5: miss 0.2 all sites snps maf 0.05
mkdir filter_minmiss20_allsites_maf5

bcftools view -m 2 --min-af 0.05:minor -O z -o ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.var.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools view --max-af 0:nref -O z -o ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.inv.vcf.gz ${WD}/filter_minmiss20/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.vcf.gz 

bcftools index ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.var.vcf.gz

bcftools index ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.inv.vcf.gz

bcftools concat --allow-overlaps ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.var.vcf.gz ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.inv.vcf.gz -O z -o ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.allSites.vcf.gz

bcftools index ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.allSites.vcf.gz

bcftools stats -s - ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.allSites.vcf.gz > ${WD}/filter_minmiss20_allsites_maf5/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels."${ChrName}".minDP3.0maxDP3.0avg.miss0.2.maf5.allSites.stats
~~~



### Phasing

First, I concatenated all the chromosome files together.

~~~bash
#!/bin/bash
## Load 
ml load bioinfo-tools bcftools/1.17

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5"

cd ${WD}

bcftools concat --allow-overlaps -f vcf_files.txt -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz
bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz
bcftools stats -s - herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz > herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.stats
~~~

Then I phased the variants with Beagle.

~~~bash
#!/bin/bash

## Load 
ml load Beagle/5.4

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5"


java -jar $BEAGLE_ROOT/beagle.jar gt=${WD}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz out=${WD}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased nthreads=8
~~~


### Annotating Variants 

Annotate the file with SNPEff:

~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J snpEff
#SBATCH -e snpEff_%A_%a.err            # File to which STDERR will be written
#SBATCH -o snpEff_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

start=$(date +%H:%M_%F)

ml load java/OracleJDK_11.0.9

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes"
cd ${WD}

## Annotate SNPs:
java -jar /proj/snic2020-2-19/private/herring/users/mafalda/software/snpEff/snpEff.jar Ch_v2.0.2.105 herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.vcf.gz > herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.vcf

bgzip herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.vcf

bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.snpEff.vcf.vcf.gz

end=$(date +%H:%M_%F)

echo "Start and End:"
echo $start
echo $end
~~~

### Concatenate with sprat outgroup

~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J addoutmaf5
#SBATCH -e addoutmaf5_%A_%a.err            # File to which STDERR will be written
#SBATCH -o addoutmaf5_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5"

cd ${WD}

bcftools merge -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.vcf.gz herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.vcf.gz /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.vcf.gz 

bcftools view -e "FMT/GT[125]=\"./.\"" herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.vcf.gz | bcftools view -e 'F_MISSING > 0.99' -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.vcf.gz

bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.vcf.gz

bcftools stats -s - herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.vcf.gz > herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.stats
~~~

There are less SNPs in these files where I filtered the vcf such that the sprat is not missing, but this is what we need for the Dsuite runs anyway since the sprat is only one indiviudal.


maf1: 10314436 SNPs
maf1 with sprat:4206888 SNPs

maf5 : 4427424 SNPs
maf5 with sprat: 1807880 SNPs

# 2025-02-14

Merge with sprat

run_sprat_merge.sh
~~~bash
#!/bin/bash
#SBATCH -A uppmax2025-2-114
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J addSprat
#SBATCH -e addSprat_%A_%a.err            # File to which STDERR will be written
#SBATCH -o addSprat_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=sferreira.mafalda@gmail.com

ml load bioinfo-tools bcftools/1.17

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_biallelic_snps_maf5"

cd ${WD}

bcftools index herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz

bcftools merge -O z -o herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.sprat.vcf.gz herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.phased.vcf.gz /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.vcf.gz
~~~

