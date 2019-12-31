#! /usr/bin/env python

# fasta_subsetter.py
# February 2015 by C Tepolt, last updated April 2018 by C Tepolt [small edit by Z Tobias December 2019]
# 
# Reads in a list of loci and a master fasta file,
# and creates a new fasta file containing
# only loci in the subset file.

Usage = "USAGE: fasta_subsetter.py MASTER_FASTA_FILE LIST_TO_SELECT [SELECT/REMOVE]"

import sys

if (len(sys.argv) < 4) or (len(sys.argv) > 4):
	print("Please provide two input files and SELECT or REMOVE after the script name.")
	print(Usage)
	exit()
else:
	MasterFileName = sys.argv[1]
	SubsetListName = sys.argv[2]
	select_type = sys.argv[3]

OutFileName = SubsetListName.rstrip(".txt") + "_" + select_type + ".fasta" #changed to .rstrip(), was removing leading tx in path

SubsetFile = open(SubsetListName, 'r')
MasterFile = open(MasterFileName, 'r')
OutFile = open(OutFileName, 'w')

SubsetCounter = 0
SubsetLoci = []

# Read in the desired loci, and print a count.

for Line in SubsetFile:
	SubsetLoci.append(Line.strip())
	SubsetCounter += 1

SubsetFile.close()	# Clean up after your damn self.

print("Read in %d lucky loci to %s from a subset file." % (SubsetCounter, select_type))

Switch = 0

for Line in MasterFile:
	if ">" in Line:
		Switch = 0
		Contig = Line.strip().split(" ")

		if select_type == "REMOVE":
			if Contig[0] in SubsetLoci:
				continue
			else:
				OutFile.write(Line)
				Switch = 1
		if select_type == "SELECT":
			if Contig[0] in SubsetLoci:
				OutFile.write(Line)
				Switch = 1
			else:
				continue

	else:
		if Switch == 1:
			OutFile.write(Line)

MasterFile.close()	# Clean up after your damn self.
OutFile.close()

print("Subset file is complete! Share and enjoy.")
print(OutFileName)
