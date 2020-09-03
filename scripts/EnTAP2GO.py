#! /usr/bin/env python

#For parsing the final annotations from EnTAP. Output table from EnTAP Beta v0.9.0-beta
#6 Jan 2020 by Zachary Tobias

import sys
import re

Usage = "EnTAP2GO.py path_to_final_annotations_lvl0.tsv path_to_output evalue_threshold[num] filter_contamination?[y/n] eggNOG_taxonomic_scope[opt]"

#recode to use argparse at some point
if (len(sys.argv) < 5):
	print("See usage:")
	print(Usage)
	exit()
elif (len(sys.argv) == 5):
    infile = sys.argv[1]
    outfile = sys.argv[2]
    eval_cutoff = float(sys.argv[3])
    exclude_contam = sys.argv[4]
elif (len(sys.argv) == 6):
    infile = sys.argv[1]
    outfile = sys.argv[2]
    eval_cutoff = float(sys.argv[3])
    exclude_contam = sys.argv[4]
    eggDB = sys.argv[5]
else:
	print("See usage:")
	print(Usage)
	exit()

print("Excluding matches above eval of " + str(eval_cutoff))

if exclude_contam == 'y':
    print("Excluding transcripts identified by EnTAP as contamination")
elif exclude_contam == 'n':
    print("Including transcripts identified by EnTAP as contamination")
else:
    print("filter_contamination argument must be 'y' or 'n', case sensitive")
    exit()
    
if "eggDB" not in locals():
    eggDB = None
    print("No eggNOG taxonomic scope selected. Retrieving hits from any matching COGs")
else:
    print("Restricting GO term retrieval to matches against " + eggDB + " COGs")


with open(infile, "r") as f_in:
    with open(outfile,"w") as f_out:
        next(f_in)
        for line in f_in:
            cols = line.split("\t") #save column values to object
            if cols[25] == "" or float(cols[25]) >= eval_cutoff: #if eval of eggNOG match exceeds threshold or is missing
                continue #skip line
            elif exclude_contam == "y" and cols[16] == "Yes": #if contaminant exclusion specified and transcript is a contaminant
                continue #skip line
            elif eggDB != None and cols[28] != eggDB: #if taxonomic scope is specified and COG doesn't match
                continue
            else:
                ID = cols[0] #get first as ID
                GO = "" #initialize empty string
                for i in range(32,35): #loop through three columns containing the GO terms
                    GO = GO + cols[i] #append column contents to string
                matches = re.findall("GO:.......", GO) #find all GO terms
                out = "" #initialize empty string
                for j, k in enumerate(matches): #loop through matches, enumerating to specify behavior on first pass 
                    if j == 0: #if first item in matches
                        out = k #append match to out
                    else: #otherwise
                        out = out + ";" + k #append match to string with ; as separator
                if len(out) == 0: #if it mapped to a KOG but no GO terms
                    continue #skip line
                else:
                    f_out.writelines(ID + "\t" + out + "\n") #write ID and GO matches to two columns
            
f_in.close()
f_out.close()
