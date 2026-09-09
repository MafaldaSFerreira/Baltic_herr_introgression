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

## 3. revisions Table 1 with Fst 

fst_calculations_Table_1.R

I need to recalculate Fst because I have note it by population, instead of collapsing all populations sharing spawning time in the Atlantic or Baltic sea. I will follow the same logic as Han et al did to decide these populations. 

Atlantic Spring contains populations from the Atlantic that spawn in the spring time, but it excludes spring spawners from the Norwegian Fjords. Includes transition zone Baltic <-> Atlantic. However, we don't have these transition individuals in the high coverage data. 
Atlantic Autumn contains populations from the Atlantic that spawn in the autumn time, including British or Irish autumn spawners.
Baltic spring contains populations from the Baltic spawning in spring.
Baltic autumn contains populations from the Baltic spawning in autumn.

These are the populations I am defining:
~~~
Fehmarn3_Fehmarn_Baltic_Autumn	Baltic_Autumn
Fehmarn44_Fehmarn_Baltic_Autumn	Baltic_Autumn
Fehmarn6_Fehmarn_Baltic_Autumn	Baltic_Autumn
Gavle100_Gavle_Baltic_Autumn	Baltic_Autumn
Gavle54_Gavle_Baltic_Autumn	Baltic_Autumn
Gavle98_Gavle_Baltic_Autumn	Baltic_Autumn
BF16_HastKar_Baltic_Spring	Baltic_Spring
BF18_HastKar_Baltic_Spring	Baltic_Spring
BF19_HastKar_Baltic_Spring	Baltic_Spring
BF21_HastKar_Baltic_Spring	Baltic_Spring
BM14_HastKar_Baltic_Spring	Baltic_Spring
BM15_HastKar_Baltic_Spring	Baltic_Spring
BM16_HastKar_Baltic_Spring	Baltic_Spring
BM19_HastKar_Baltic_Spring	Baltic_Spring
F1_HastKar_Baltic_Spring	Baltic_Spring
F2_HastKar_Baltic_Spring	Baltic_Spring
F3_HastKar_Baltic_Spring	Baltic_Spring
F4_HastKar_Baltic_Spring	Baltic_Spring
F5_HastKar_Baltic_Spring	Baltic_Spring
F6_HastKar_Baltic_Spring	Baltic_Spring
NorthSea13_NorthSea_Atlantic_Autumn	Atlantic_Autumn
NorthSea19_NorthSea_Atlantic_Autumn	Atlantic_Autumn
NorthSea34_NorthSea_Atlantic_Autumn	Atlantic_Autumn
NSSH33_Norway_Atlantic_Spring	Atlantic_Spring
NSSH34_Norway_Atlantic_Spring	Atlantic_Spring
NSSH36_Norway_Atlantic_Spring	Atlantic_Spring
144Sbs344_Canada_Atlantic_Autumn	Atlantic_Autumn
144Sbs349_Canada_Atlantic_Autumn	Atlantic_Autumn
14F4TL404_Canada_Atlantic_Autumn	Atlantic_Autumn
14F4TL415_Canada_Atlantic_Autumn	Atlantic_Autumn
14F4WK306_Canada_Atlantic_Autumn	Atlantic_Autumn
14F4WK316_Canada_Atlantic_Autumn	Atlantic_Autumn
15F2J602_Canada_Atlantic_Autumn	Atlantic_Autumn
15F2J603_Canada_Atlantic_Autumn	Atlantic_Autumn
15F2J606_Canada_Atlantic_Autumn	Atlantic_Autumn
15F2J616_Canada_Atlantic_Autumn	Atlantic_Autumn
15F3K601_Canada_Atlantic_Autumn	Atlantic_Autumn
15F3K602_Canada_Atlantic_Autumn	Atlantic_Autumn
15F3K609_Canada_Atlantic_Autumn	Atlantic_Autumn
15F3K621_Canada_Atlantic_Autumn	Atlantic_Autumn
15F4XRsb625_Canada_Atlantic_Autumn	Atlantic_Autumn
15F4XRsb626_Canada_Atlantic_Autumn	Atlantic_Autumn
15F5Y514-617_Canada_Atlantic_Autumn	Atlantic_Autumn
15F5Y514-620_Canada_Atlantic_Autumn	Atlantic_Autumn
F3L312_Canada_Atlantic_Autumn	Atlantic_Autumn
F3L337_Canada_Atlantic_Autumn	Atlantic_Autumn
F4TH327_Canada_Atlantic_Autumn	Atlantic_Autumn
F4TH337_Canada_Atlantic_Autumn	Atlantic_Autumn
F4XQgb402_Canada_Atlantic_Autumn	Atlantic_Autumn
F4XQgb408_Canada_Atlantic_Autumn	Atlantic_Autumn
12S4Rsv38_Canada_Atlantic_Spring	Atlantic_Spring
12S4Rsv41_Canada_Atlantic_Spring	Atlantic_Spring
14S3Ps233_Canada_Atlantic_Spring	Atlantic_Spring
14S3Ps265_Canada_Atlantic_Spring	Atlantic_Spring
15S3K404_Canada_Atlantic_Spring	Atlantic_Spring
15S3K410_Canada_Atlantic_Spring	Atlantic_Spring
15S3K430_Canada_Atlantic_Spring	Atlantic_Spring
15S3K449_Canada_Atlantic_Spring	Atlantic_Spring
16S4Pla78_Canada_Atlantic_Spring	Atlantic_Spring
16S4Pla79_Canada_Atlantic_Spring	Atlantic_Spring
16S6Pla71_Canada_Atlantic_Spring	Atlantic_Spring
16S6Pla72_Canada_Atlantic_Spring	Atlantic_Spring
16SBDO106_Canada_Atlantic_Spring	Atlantic_Spring
16SBDO107_Canada_Atlantic_Spring	Atlantic_Spring
S3Ps246_Canada_Atlantic_Spring	Atlantic_Spring
S3Ps259_Canada_Atlantic_Spring	Atlantic_Spring
S4TH231_Canada_Atlantic_Spring	Atlantic_Spring
S4TH244_Canada_Atlantic_Spring	Atlantic_Spring
S4TM205_Canada_Atlantic_Spring	Atlantic_Spring
S4TM211_Canada_Atlantic_Spring	Atlantic_Spring
Z12_IsleofMan_Atlantic_Autumn	Atlantic_Autumn
Z14_IsleofMan_Atlantic_Autumn	Atlantic_Autumn
Z4_IsleofMan_Atlantic_Autumn	Atlantic_Autumn
~~~

I might use pixy but also vcftools because with the later I can calculate a per SNP Fst.

run_pixy_rev_260208.sh
~~~
#!/bin/bash -l

#SBATCH -A uppmax2025-2-453
#SBATCH -n 2
#SBATCH -M pelle
#SBATCH --array=1-26:1
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL
#SBATCH -J pixy
#SBATCH -e pixy_%A_%a.err
#SBATCH -o pixy_%A_%a.out

## Load required modules
ml load bioinfo-tools pixy/1.2.5.beta1
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

pixy --stats pi fst dxy --populations ${PopDIR}/clusters_rev_260208.txt --vcf ${input_vcf_dir}/herring_sentieon_125ind_231031.newID.filter.setGT.noIndels.${ChrName}.minDP3.0maxDP3.0avg.miss0.2.${FILTER}.allSites.vcf.gz --window_size 20000 --n_cores 2 --output_folder results/clusters_rev_260208_${FILTER} --output_prefix ${run}.clusters_rev_260208.${ChrName}.${ARG2}.20kb.popgenpixy.out
~~~

sbatch --export=ALL,ARG1="chromosomes_maf0.05",ARG2="maf5" run_pixy_rev_260208.sh


