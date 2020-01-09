#! /usr/bin/env python

#Jan 2020 by Zachary Tobias. Reformats UniRef90 database header taxonomic information so it can be parsed by EnTAP for use in contaminant filtering and taxonomic weighting

Usage = "USAGE: reformat_uniref.py PATH_IN PATH_OUT"

import sys
import re

if (len(sys.argv) < 3) or (len(sys.argv) > 3):
	print("Please provide input file path and output file path")
	print(Usage)
	exit()
else:
	infile = sys.argv[1]
	outfile = sys.argv[2]

unformatted = open(infile, 'r')
formatted = open(outfile, 'w')

for line in unformatted:
    if line.startswith('>',0,1):
        line = re.sub('Tax\\=','[',line)
        line = re.sub(' TaxID','] TaxID', line)
    formatted.write(line)
    
unformatted.close()
formatted.close()

