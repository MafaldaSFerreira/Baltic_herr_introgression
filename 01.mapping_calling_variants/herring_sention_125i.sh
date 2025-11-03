#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# starting with BAM files
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
data_dir=/crex/proj/snic2020-2-19/private/herring/alignment/125_individuals
# MF: I dont know what this dir is, but assuming is the output, I'll give 
# a potential output.
TOPDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call

# Get the bam file from the command line
BAM=$1
# Set the sample name
SAMPLE=$(basename $BAM)
SAMPLE_NAME=${SAMPLE/.MD.RG.bam/}
OUTPUT_SAMPLE=${SAMPLE/.MD.RG.bam/.g.vcf.gz}

# Update with the location of the reference data files
FASTA=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files_v02/Uppsala_University_node-52.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files_v02/sentieon-genomics-202112.07

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

platform="ILLUMINA"

# Other settings
nt=10 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$TOPDIR/$SAMPLE_NAME #ede edited this to be sample
mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 6. HC Variant caller
# Note: Sentieon default setting matches versions before GATK 3.7.
# Starting GATK v3.7, the default settings have been updated multiple times.
# Below shows commands to match GATK v3.7 - 4.1
# Please change according to your desired behavior.
# ******************************************

# Matching GATK 4.1
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r ${FASTA} -t $nt -i ${BAM} --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 ${OUTPUT_SAMPLE}.g.vcf.gz
