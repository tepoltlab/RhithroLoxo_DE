#! /usr/bin/env python

import sys
import re

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as f_in:
    with open(outfile,"a") as f_out:
        next(f_in)
        for line in f_in:
            cols = line.split("\t")
            ID = cols[0]
            GO = ""
            for i in range(0,3):
                GO = GO + cols[i+9]
            matches = re.findall("GO:.......", GO)
            out = ""
            for j, k in enumerate(matches):
                if j == 0:
                    out = k
                out = out + ";" + k
            f_out.writelines(ID + "\t" + out + "\n")
            
f_in.close()
f_out.close()
