### Reference Documentation

* [JunctionSeq Example Walkthrough PDF](http://hartleys.github.io/JunctionSeq/doc/example-walkthrough.pdf)

### Array Mapping Test

```bash
# Initial test to check array task ID mapping against fastq files
FW_READS=(./*e.fastp.t1.bbduk.t2.fastq.gz)
FW_READ=${FW_READS[$SLURM_ARRAY_TASK_ID]}

echo ${FW_READ}
echo ${${SLURM_ARRAY_TASK_ID}}

```

### Final Trimming Pipeline (SLURM Batch Script)

```bash
#!/bin/bash -l
 
#SBATCH -A snic2021-5-400
#SBATCH -p core
#SBATCH -n 2
#SBATCH --array=1-12:1
#SBATCH -t 10:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J fastqc
#SBATCH -e fastqc_%A_%a.err
#SBATCH -o fastqc_%A_%a.out

# Target specific trimmed read files using the SLURM Task ID
FW_READ=$(ls *e.fastp.t1.bbduk.t2.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo ${FW_READ}

# Run initial FastQC check
fastqc ${FW_READ}

# ----------------------------------------------------
# Core Preprocessing Pipeline (fastp -> BBDuk -> Trimmomatic)
# ----------------------------------------------------

# Load required cluster modules
ml load bioinfo-tools fastp
ml load bioinfo-tools FastQC/0.11.9
ml load bioinfo-tools bbmap/38.61b

# Define Input Reads
# Read 1
FW_READ=$(ls *_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Read 2
RV_READ=${FW_READ/_R1.fastq.gz/_R2.fastq.gz}

#### FASTP DEFINITIONS ####
# Output read 1
FW_READ_OUT_1=${FW_READ/_R1.fastq.gz/_R1_FASTP.fastq.gz}
# Output read 2
RV_READ_OUT_1=${FW_READ/_R1.fastq.gz/_R2_FASTP.fastq.gz}
# HTML
HTML=${FW_READ/_R1.fastq.gz/_fastp.html}
# JSON
JSON=${FW_READ/_R1.fastq.gz/_fastp.json}

#### BBDUK DEFINITIONS ####
# Output read 1
FW_READ_OUT_2=${FW_READ/_R1.fastq.gz/_R1_FASTP_BBDUK.fastq.gz}
# Output read 2
RV_READ_OUT_2=${FW_READ/_R1.fastq.gz/_R2_FASTP_BBDUK.fastq.gz}
# Output read 1 failed
FW_READ_OUT_2_FAILED=${FW_READ/_R1.fastq.gz/_R1_FASTP_BBDUK_FAILED.fastq.gz}
# Output read 2 failed
RV_READ_OUT_2_FAILED=${FW_READ/_R1.fastq.gz/_R2_FASTP_BBDUK_FAILED.fastq.gz}
# Stats
STATS=${FW_READ/_R1.fastq.gz/_FASTP_BBDUK.stats}

#### TRIMMOMATIC DEFINITIONS ####
# Output read 1 PE
FW_READ_OUT_PE_3=${FW_READ/_R1.fastq.gz/_R1_FASTP_BBDUK_TRIMM_PE.fastq.gz}
# Output read 2 PE
RV_READ_OUT_PE_3=${FW_READ/_R1.fastq.gz/_R2_FASTP_BBDUK_TRIMM_PE.fastq.gz}
# Output read 1 SE
FW_READ_OUT_SE_3=${FW_READ/_R1.fastq.gz/_R1_FASTP_BBDUK_TRIMM_SE.fastq.gz}
# Output read 2 SE
RV_READ_OUT_SE_3=${FW_READ/_R1.fastq.gz/_R2_FASTP_BBDUK_TRIMM_SE.fastq.gz}


#### VERIFY FILES ####
echo ${FW_READ}
echo ${RV_READ}
echo ${FW_READ_OUT_1}
echo ${RV_READ_OUT_1}
echo ${FW_READ_OUT_2}
echo ${RV_READ_OUT_2}
echo ${FW_READ_OUT_2_FAILED}
echo ${RV_READ_OUT_2_FAILED}
echo ${FW_READ_OUT_PE_3}
echo ${RV_READ_OUT_PE_3}
echo ${FW_READ_OUT_SE_3}
echo ${RV_READ_OUT_SE_3}

#### RUN FASTP ####
#fastp --in1 ${FW_READ} --in2 ${RV_READ} --out1 ${FW_READ_OUT_1} --out2 ${RV_READ_OUT_1} -h ${HTML} -j ${JSON} --trim_poly_g

#fastqc ${FW_READ_OUT_1}
#fastqc ${RV_READ_OUT_1}

#### RUN BBDUK ####
#bbduk.sh in=${FW_READ_OUT_1} in2=${RV_READ_OUT_1} out=bbduk/${FW_READ_OUT_2} out2=bbduk/${RV_READ_OUT_2} outm=bbduk/${FW_READ_OUT_2_FAILED} outm2=bbduk/${RV_READ_OUT_2_FAILED} stats=bbduk/${STATS} tbo=t tpe=t ktrim=r ref=/sw/bioinfo/bbmap/38.61b/rackham/resources/adapters.fa

# Post-BBDuk Quality Check
fastqc bbduk/${FW_READ_OUT_2}
fastqc bbduk/${RV_READ_OUT_2}

#### RUN TRIMMOMATIC ####
java -jar /proj/snic2020-2-19/private/herring/users/mafalda/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 bbduk/${FW_READ_OUT_2} bbduk/${RV_READ_OUT_2} trimmomatic/${FW_READ_OUT_PE_3} trimmomatic/${FW_READ_OUT_SE_3} trimmomatic/${RV_READ_OUT_PE_3} trimmomatic/${RV_READ_OUT_SE_3} SLIDINGWINDOW:5:15 LEADING:3 TRAILING:3 MINLEN:40

# Post-Trimmomatic Quality Check
fastqc trimmomatic/${FW_READ_OUT_PE_3}
fastqc trimmomatic/${RV_READ_OUT_PE_3}
fastqc trimmomatic/${FW_READ_OUT_SE_3}
fastqc trimmomatic/${RV_READ_OUT_SE_3}

```

---
### RSEM Reference Genome Preparation

The reference index was built using both STAR and Bowtie2 backends.

```bash
# Load required bioinformatics software suites
ml load bioinfo-tools rsem/1.3.3 star/2.7.9a bowtie2/2.3.5.1

# 1. Prepare reference index utilizing STAR alignment options
rsem-prepare-reference --gtf /proj/snic2020-2-19/private/herring/users/mafalda/Herring_ISOSeq/Ensembl_gff3/Clupea_harengus.Ch_v2.0.2.104.gtf --star /proj/snic2020-2-19/private/herring/users/mafalda/Herring_ISOSeq/Ensembl_gff3/Ch_v2.0.2.mod.fasta star_alignments/Ch_v2_0_2

# 2. Prepare reference index utilizing Bowtie2 alignment options
rsem-prepare-reference --gtf /proj/snic2020-2-19/private/herring/users/mafalda/Herring_ISOSeq/Ensembl_gff3/Clupea_harengus.Ch_v2.0.2.104.gtf --bowtie2 /proj/snic2020-2-19/private/herring/users/mafalda/Herring_ISOSeq/Ensembl_gff3/Ch_v2.0.2.mod.fasta bowtie_alignments/Ch_v2_0_2

```

---

> 💡 **Library Specification Note**
> Per communication from Leif, libraries were prepared using the **TruSeq stranded mRNA library protocol** with polyA selection. This dictates a reverse strandedness configuration down the line.

### File Management & Permissions

To ensure data integrity and avoid accidental deletion/modifications, intermediate files were secured:

```bash
# Set permissions of all final trimmomatic fastq outputs to read-only
chmod a=r *fastq.gz

```

*Note: Linked files were structured inside the `rsem/` project folder.*

### Alignment & Expression Quantification

Quantification step performed via RSEM using the STAR alignment option.

```bash
# Load environment modules
ml load bioinfo-tools rsem/1.3.3 star/2.7.9a bowtie2/2.3.5.1

# Parse array task inputs from fully trimmed paired-end reads
# Read 1
FW_READ=$(ls *_R1_FASTP_BBDUK_TRIMM_PE.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Read 2
RV_READ=${FW_READ/_R1_FASTP_BBDUK_TRIMM_PE.fastq.gz/_R2_FASTP_BBDUK_TRIMM_PE.fastq.gz}
# Output Prefixes
OUT_STAR=${FW_READ/_R1_FASTP_BBDUK_TRIMM_PE.fastq.gz/_star_rsem_results}

# Calculate expression values using reverse strandedness logic
rsem-calculate-expression --star-gzipped-read-file --paired-end --calc-ci --strandedness reverse -p 16 --star ${FW_READ} ${RV_READ} star_alignments/Ch_v2_0_2 star_results/${OUT_STAR}

```

