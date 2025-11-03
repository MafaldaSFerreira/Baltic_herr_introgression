# ABBA-BABA analysis

run_dtrios_maf5_23-11.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -t 1
#SBATCH -t 10:00:00
#SBATCH -J dsuite    
#SBATCH -e dsuite_%A_%a.err            # File to which STDERR will be written
#SBATCH -o dsuite_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/dsuite/"
DATA_WD=${WD}"/data"
OUTPUT_WD=${WD}"/results_2023-11_maf0.05_dtrios_sprat"
POP_WD=${WD}"/populations"

#JKNUM=$ARG1

mkdir ${OUTPUT_WD}
cd ${OUTPUT_WD}

run=$(date +%H%M_%F)

/proj/snic2020-2-19/private/herring/users/mafalda/software/Dsuite/Build/Dsuite Dtrios -o all_vs_all_clusters_v03_maf5_JK1000 -n run_$run --JKwindow 1000 ${DATA_WD}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.vcf.gz ${POP_WD}/clusters_v03.txt
~~~

# Dinvestigate

run_dinvestigate_v02_maf5_2023-11-30.sh
~~~bash
#!/bin/bash
#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -t 1
#SBATCH -t 10:00:00
#SBATCH -J dinv   
#SBATCH -e dinv_%A_%a.err            # File to which STDERR will be written
#SBATCH -o dinv_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=amafaldasferreira@gmail.com

WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/dsuite/"
DATA_WD=${WD}"/data"
OUTPUT_WD=${WD}"/results_2023-11-30_maf5_dinvestigate_sprat_outgroup"
POP_WD=${WD}"/populations"

JKNUM=$ARG1
STEP=$ARG2

mkdir ${OUTPUT_WD}
cd ${OUTPUT_WD}

run=$(date +%H%M_%F)

/proj/snic2020-2-19/private/herring/users/mafalda/software/Dsuite/Build/Dsuite Dinvestigate -n run_${run}_w${JKNUM}_s${STEP} --window=${JKNUM},${STEP} ${DATA_WD}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.minDP3.0maxDP3.0avg.miss0.2.biallelic.maf5.sprat.filtered.vcf.gz ${POP_WD}/clusters_v03.txt ${POP_WD}/trios_v02.txt
~~~


sbatch --export=ALL,ARG1="50",ARG2="25" run_dinvestigate_v02_maf5_2023-11-30.sh
