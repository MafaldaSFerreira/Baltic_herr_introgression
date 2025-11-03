# Genetic divergence scans (dxy) with pixy

We used pixy to run genetic divergence scans in two steps:

- 1. Using all individuals assigned to a specific population group
- 2. Using homozygotes for each introgression region, to determine dxy between Baltic spring and White Sea individuals

## 1. dxy scan

For this analysis, we start with all sites vcf files, including invariant sites and variants sites (biallelic SNPs, maf > 0.05 and %miss 20).

Populations are defined as in `cluster_v03.txt`

~~~bash
ln -s /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/variant_call/filtered_vcfs/chromosomes/filter_minmiss20_allsites_maf5/* ./
~~~

~~~bash
#!/bin/bash -l

#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -M rackham
#SBATCH --array=1-26:1
#SBATCH -t 5:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J pixy
#SBATCH -e pixy_%A_%a.err
#SBATCH -o pixy_%A_%a.out

## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1 bcftools/1.17
## Determine chromosome
ChrName=chr${SLURM_ARRAY_TASK_ID}

# Read arguments from the command line
DIR=$ARG1
FILTER=$ARG2

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy"
PopDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy/population_files"
input_vcf_dir="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy/"${DIR}

cd ${WD}

run=$(date +%H%M_%F)

tabix ${input_vcf_dir}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.minDP3.0maxDP3.0avg.miss0.2.${FILTER}.allSites.vcf.gz

pixy --stats pi fst dxy --populations ${PopDIR}/clusters_v03.txt --vcf ${input_vcf_dir}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.minDP3.0maxDP3.0avg.miss0.2.${FILTER}.allSites.vcf.gz --window_size 20000 --n_cores 2 --output_folder results/clusters_v03_${FILTER} --output_prefix ${run}.clusters_v03.${ChrName}.${ARG2}.20kb.popgenpixy.out
~~~

We run this slurm script as:

~~~bash
sbatch --export=ALL,ARG1="chromosomes_maf0.05",ARG2="maf5" run_pixy_2023-11.sh
~~~

Use `wg_tables.R` to create wg files:

~~~bash
wg_tables.R 20kb.popgenpixy.out_fst.txt 2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_fst.txt
wg_tables.R 20kb.popgenpixy.out_dxy.txt 2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_dxy.txt
wg_tables.R 20kb.popgenpixy.out_pi.txt 2313_2023-11-24.clusters_v03.wg.maf5.20kb.popgenpixy.out_pi.txt
~~~


## 2. dxy scan for homozygotes

I will add Baltic Autumn and White Sea individuals as contrast for dxy and fst.

I extracted the homozygotes in R using `find_homozygote_indv_dxy.R`

Then I added the following lines to each file:

HWS41_KandalakshaBay_WhiteSea_Spring	WhiteSea
HWS42_KandalakshaBay_WhiteSea_Spring	WhiteSea
HWS43_KandalakshaBay_WhiteSea_Spring	WhiteSea
HWS44_KandalakshaBay_WhiteSea_Spring	WhiteSea
HWS31_WhiteSea_WhiteSea	WhiteSea
HWS32_WhiteSea_WhiteSea	WhiteSea
HWS33_WhiteSea_WhiteSea	WhiteSea
HWS34_WhiteSea_WhiteSea	WhiteSea
HWS51_KandalakshaBay_WhiteSea_Summer	WhiteSea
HWS52_KandalakshaBay_WhiteSea_Summer	WhiteSea
HWS53_KandalakshaBay_WhiteSea_Summer	WhiteSea
HWS54_KandalakshaBay_WhiteSea_Summer	WhiteSea
Fehmarn3_Fehmarn_Baltic_Autumn	Baltic_Autumn
Fehmarn44_Fehmarn_Baltic_Autumn	Baltic_Autumn
Fehmarn6_Fehmarn_Baltic_Autumn	Baltic_Autumn
Gavle100_Gavle_Baltic_Autumn	Baltic_Autumn
Gavle54_Gavle_Baltic_Autumn	Baltic_Autumn
Gavle98_Gavle_Baltic_Autumn	Baltic_Autumn

~~~bash
for i in $(ls chr*.txt); do cat $i extra_individuals_for_homoz_files.txt > ${i/_popfile.txt/_WhiteSea_BalticAutumn_popfile.txt}; done

# Removed extra files:
rm chr*_Baltic_Spring_popfile.txt
~~~

run_pixy_homozygotes_2024-01-19.sh
~~~bash
#!/bin/bash -l

#SBATCH -A naiss2023-5-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -M rackham
#SBATCH --array=1-67:1
#SBATCH -t 5:00:00
#SBATCH --mail-user=amafaldasferreira@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -J pixy
#SBATCH -e pixy_%A_%a.err
#SBATCH -o pixy_%A_%a.out

## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1 bcftools/1.17

# Read arguments from the command line
DIR=$ARG1
FILTER=$ARG2

# Define directories
WD="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy"
POPDIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy/population_files"
INPUT_VCF_DIR="/proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy/"${DIR}

## Determine File and Chromosome to run:
FILECHR=$(cat ${POPDIR}/clusters_v05_homozygotes.txt | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}')

FILE=$(echo $FILECHR | cut -f1 -d" ")
CHR=$(echo $FILECHR | cut -f2 -d" ")

echo $FILECHR
echo $FILE 
echo $CHR

# Run Pixy:
cd ${WD}

pixy --stats pi fst dxy --populations ${POPDIR}/scan1_v01_baltic_alt_ref_summary_filter2_cov7_homozygotes/${FILE} --vcf ${INPUT_VCF_DIR}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${CHR}.minDP3.0maxDP3.0avg.miss0.2.${FILTER}.allSites.vcf.gz --window_size 20000 --n_cores 2 --output_folder results/scan1_v01_baltic_alt_ref_summary_filter2_cov7_homozygotes/${FILE/.txt/}_${FILTER} --output_prefix ${FILE/.txt/}.${CHR}.${FILTER}.20kb.popgenpixy.out
~~~



~~~bash
cd /proj/snic2020-2-19/private/herring/users/mafalda/Introgression/pixy/scripts

sbatch --export=ALL,ARG1="chromosomes_maf0.05",ARG2="maf5" run_pixy_homozygotes_2024-01-19.sh
~~~
