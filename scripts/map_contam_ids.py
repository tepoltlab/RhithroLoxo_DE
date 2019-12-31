import pandas as pd
import sys
from ete3 import NCBITaxa
ncbi = NCBITaxa()

INPUT = sys.argv[1]
OUTPUT = sys.argv[2]

hits = pd.read_table(INPUT,header=None)

candidate_contam = [2,2157,10239,6157,6231,4751,33630,33090,2763,554915,543769,33634,37909,1458140] #list contam taxids

contam = []
for i in hits.index: 
    contig = hits.loc[i,0]
    taxid = hits.loc[i,13]
    lineage = ncbi.get_lineage(taxid)
    if bool(set(lineage) & set(candidate_contam)):
        contam.append(contig)

with open(OUTPUT, 'w') as file:
    for i in contam:
        file.write(">"+'%s\n' % i)