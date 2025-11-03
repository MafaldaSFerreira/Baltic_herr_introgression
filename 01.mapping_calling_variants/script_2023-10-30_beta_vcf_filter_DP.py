#!/usr/bin/python



"""
usage: provide vcf to be filtered for DP per sample (set GT to ./. if min < DP < max not in the desired range) - provide also list of avg genome cov per sample in same order as in vcf (alphabetically) and value for min cov (at least 2?) and value to multiply avg cov with for max filter (*3?)

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 cpDNA mtDNA unplaced; 
do grep "$chr" average_coverage.txt |sort -u > 629_avgcov_"$chr".txt;
python3 script_2022-12-15_beta_vcf_filter_DP.py calling_all/629beta_"$chr"_filter.gt.vcf 629beta_avgcov_"$chr".txt 2 2.5;
#<INT:INT=minDP> <INT:avgcov*INT=maxDP>;
done;

Mafalda modified 2023-10-30: the script will through an error if the DP field has a . - which it shouldn't have but these exist! So, let's the script to ignore these positions
"""

import sys

#read in vcf file to be filtered
infile=open(sys.argv[1],'r')

#minDP
minDP=float(sys.argv[3])

#maxDP
maxDP=float(sys.argv[4])

#generate output filename
outfile=".".join(sys.argv[1].split(".")[0:-1])+".minDP"+str(minDP)+"maxDP"+str(maxDP)+"avg.vcf"

#read in avg genome cov per sample - same order as in input vcf!
cov_in=open(sys.argv[2],'r')
cov=[]
for line in cov_in:
	cov.append(float(line.strip()))

#print(cov,len(cov))

#process data - filter vcf with min and max cov per sample
with open(outfile,'w') as out:
	for line in infile.readlines():
		if line.startswith("#CHROM"):
			out.write(line)
		elif line.startswith("#"):
			out.write(line)
		else:
			newline=["\t".join(line.split("\t")[0:9])]
			genotypes=line.strip().split("\t")[9:]
			print(genotypes)
			idx=-1
			for gt in genotypes:
				#print(genotypes)
				#print(len(genotypes))
				idx=idx+1
				if gt=="./.":
					newline.append(gt)	
				else:
					dp=gt.split(":")[2]
					if dp==".":
						newline.append(gt)
					elif int(dp)>=minDP:
						#print(idx)
						if int(dp) > maxDP*cov[idx]:
							newline.append("./."+gt[3:])
						else:
							newline.append(gt)
					else:
						newline.append("./."+gt[3:])
			out.write("\t".join(newline)+"\n")


