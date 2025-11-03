# Long-read mapping and genotype calling for European sprat

Here we map and call genotypes with long-read data for the European sprat data published in Pettersson et al (2023).

I opted for creating a mamba environment in this case. 

~~~yml
channels:
- bioconda
dependencies:
- samtools=1.17
- bcftools=1.17
- bedtools=2.31.0
- minimap2=2.26
- bedtools=2.29.2
~~~

nstalled the environment:
~~~bash
mamba env create -f 01_mapping.yaml --prefix /proj/snic2020-2-19/private/herring/users/mafalda/hifi_mapping
mamba env update -f 01_mapping.yaml
~~~

### Mapping with minimap2

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 3
#SBATCH -M rackham
#SBATCH -t 24:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J minimap2
#SBATCH -e minimap2_%A_%a.err
#SBATCH -o minimap2_%A_%a.out

conda activate /proj/snic2020-2-19/private/herring/users/mafalda/hifi_mapping

i="m64077_sprat"
genome="/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta"
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup"

minimap2 -ayYL --MD --eqx -x asm20 -R @RG\\tID:${i}\\tSM:${i}\\tPL:SequelII ${genome} ${WD}/raw/m64077_201204_132418.Q20.fastq.gz | samtools sort -o ${WD}/output_minimap2/m64077_201204_132418.bam 2> ${WD}/output_minimap2/m64077_201204_132418.log
~~~

### Genotype calling

I considerd using DNAScope to call genotypes for long-read data but this software does not output invariant sites. So, in the end, I used GATK.

#### Haplotyper
~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 5
#SBATCH -M rackham
#SBATCH -t 72:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J senteion
#SBATCH -e senteion_%A_%a.err
#SBATCH -o senteion_%A_%a.out

ml python/2.7.15

# *******************************************
# Script to perform DNA seq variant calling
# adapted from senteion original script
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
data_dir=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_minimap2
TOPDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup

# Update with the location of the reference data files
fasta="/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta"

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files_v02/Uppsala_University_node-51.lic 
echo $SENTIEON_LICENSE

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files_v02/sentieon-genomics-202112.07

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
sample="m64077_201204_132418"

# Other settings
nt=5 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$TOPDIR/output_gatk_variant_call
#mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# Matching GATK 4.1
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i $data_dir/${sample}.bam --algo Haplotyper --ploidy 2 --call_conf 30 --pcr_indel_model aggressive --emit_mode gvcf --emit_conf 30 ${sample}.HC30.conf.vcf.gz
~~~

#### Genotyper

~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -t 1
#SBATCH -t 72:00:00
#SBATCH -J genotyper
#SBATCH -e genotyper_%A_%a.err            # File to which STDERR will be written
#SBATCH -o genotyper_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com


TOPDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup
OUTPUTDIR=/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call

cd $TOPDIR

REFERENCE=/crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

#sentieon options
# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE=/domus/h1/mafaldaf/sentieon_files_v02/Uppsala_University_node-52.lic

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR=/domus/h1/mafaldaf/sentieon_files_v02/sentieon-genomics-202112.07

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

# wc -l gvcf_count_V2.txt

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REFERENCE -t 1 --algo GVCFtyper -v ${OUTPUTDIR}/m64077_201204_132418.HC30.conf.vcf.gz --call_conf 60 --emit_conf 30 --emit_mode CONFIDENT ${OUTPUTDIR}/m64077_201204_132418.HC30.ALLSITES.vcf.gz 
~~~


### Genotype filter

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH -t 5:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J set_filters
#SBATCH -e set_filters_%A_%a.err
#SBATCH -o set_filters_%A_%a.out

## Load required modules
ml load bioinfo-tools GATK/4.3.0.0 bcftools/1.17

gatk --java-options "-Xmx29g" VariantFiltration  \
    -R /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta \
    -V /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.vcf.gz \
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
    -G-filter "DP < 6" \
    -G-filter-name "gtDPMin6" \
    -G-filter "DP > 69" \
    -G-filter-name "gtDPMax63" \
    -G-filter "vc.isSNP() && GQ < 20" \
    -G-filter-name "gtGQ20" \
    -O /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.vcf.gz

#set filtered GT to no call
gatk --java-options "-Xmx29g" SelectVariants  \
   -R /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta \
   -V /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.vcf.gz \
   --set-filtered-gt-to-nocall \
   -O /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.vcf.gz


bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.vcf.gz /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.vcf.gz 

bcftools index /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.vcf.gz
~~~


### Create a consensus file

~~~bash
#!/bin/bash -l
 
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 1
#SBATCH -M rackham
#SBATCH --array=1-26:1
#SBATCH -t 3:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J consensus
#SBATCH -e consensus_%A_%a.err
#SBATCH -o consensus_%A_%a.out

## Load required modules
conda activate /proj/snic2020-2-19/private/herring/users/mafalda/hifi_mapping

ml load bioinfo-tools BEDTools/2.29.2

## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

INPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_gatk_variant_call"
OUTPUT_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/outgroup/output_consensus"

bcftools view -O z -o ${INPUT_DIR}/chromosome_vcfs/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.${ChrName}.vcf.gz ${INPUT_DIR}/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.vcf.gz ${ChrName} 

bcftools index ${INPUT_DIR}/chromosome_vcfs/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.${ChrName}.vcf.gz

bcftools view ${INPUT_DIR}/chromosome_vcfs/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.${ChrName}.vcf.gz | bcftools view -H -g ^miss | awk '{OFS="\t"}{print $1,$2-1,$2}' > ${OUTPUT_DIR}/${ChrName}.sprat.called.bed

grep -w ${ChrName} /proj/snic2020-2-19/private/herring/users/mafalda/Inversion_project/Consensus/Ch_v2.0.2.chromsizes > ${OUTPUT_DIR}/${ChrName}.chromsize

bedtools complement -i ${OUTPUT_DIR}/${ChrName}.sprat.called.bed -g ${OUTPUT_DIR}/${ChrName}.chromsize > ${OUTPUT_DIR}/${ChrName}.sprat.nocall.bed

samtools faidx /crex/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta ${ChrName} | bcftools consensus ${INPUT_DIR}/chromosome_vcfs/m64077_201204_132418.HC30.ALLSITES.VF.setGT.noIndels.${ChrName}.vcf.gz > ${OUTPUT_DIR}/EuSprat.${ChrName}.fasta

bedtools maskfasta -fi ${OUTPUT_DIR}/EuSprat.${ChrName}.fasta -bed ${OUTPUT_DIR}/${ChrName}.sprat.nocall.bed -fo ${OUTPUT_DIR}/EuSprat.${ChrName}.consensus.fa
~~~

