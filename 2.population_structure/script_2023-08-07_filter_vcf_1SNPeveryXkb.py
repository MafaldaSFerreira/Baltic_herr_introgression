#!/usr/bin/python

#07.08.23
#sf
#script to filter vcf for 1 SNP every x bases - x is provided as 2nd argument, vcf is 1st argument

import sys, gzip

#read vcf
#inf=open(sys.argv[1],'r')
inf=gzip.open(sys.argv[1], mode='rt')#, compresslevel=9, encoding=None, errors=None, newline=None)

#assign snp distance to variable
win=int(sys.argv[2])

#name output file
outfile=".".join(sys.argv[1].split("/")[-1].split(".")[0:-1])+".1SNPevery"+str(win)+".vcf.gz"

#filter vcf:
old="NONE"
pos=0
#print(str(win))
with gzip.open(outfile, "wt") as out:
#	for line in inf.readlines():
	for line in inf.readlines():
		#print(line)
		if line.startswith("#"):
			#print("A",line)
			out.write(line)
		elif old=="NONE":
			old=line.split("\t")[0]
			pos=int(line.split("\t")[1])
		elif line.split("\t")[0]!=old:
			old=line.split("\t")[0]
			pos=int(line.split("\t")[1])
			out.write(line)
			#print("B",line)
		elif int(line.split("\t")[1])>pos+win:
			pos=int(line.split("\t")[1])
			out.write(line)
			#print("C",line)
