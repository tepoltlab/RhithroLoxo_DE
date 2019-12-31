#! /usr/bin/env python

# list-trinity-dups-to-remove.py
# Last updated: 20 March 2019 by C Tepolt
# 
# Reads in a raw Trinity transcriptome assembly file in .fasta format
# Creates a list of "duplicate" isoforms to be removed (retaining only the longest for each gene)`

Usage = "USAGE: list-trinity-dups-to-remove.py TRINITY_FASTA_FILE"

import sys

if (len(sys.argv) < 2) or (len(sys.argv) > 2):
	print("Please provide one input file after the script name.")
	print(Usage)
	exit()
else:
	MasterFileName = sys.argv[1]

OutFileName = MasterFileName + "_dups.txt"

MasterFile = open(MasterFileName, 'r')
OutFile = open(OutFileName, 'w')

GeneID = ""
GeneHolder = ""
Name = ""
Length = ""
LengthHolder = 0
DupCounter = 0 

for Line in MasterFile:
	if '>' in Line:
		Elements = Line.strip().split(' ')
		Elems = Elements[0].split('_')
		GeneID = '_'.join(Elems[0:4])

		Elemenos = Elements[1].split('=')
		Length = int(Elemenos[1])

		if GeneID != GeneHolder:
			GeneHolder = GeneID
			Name = Elements[0]
			LengthHolder = Length

		else:
			if Length < LengthHolder:
				OutFile.write(Elements[0] + '\n')
				DupCounter += 1
			else:
				OutFile.write(Name + '\n')
				LengthHolder = Length
				Name = Elements[0]
				DupCounter += 1

MasterFile.close()	# Clean up after your damn self.
OutFile.close()

print("Listed %d trinity duplicates to remove." % (DupCounter))
