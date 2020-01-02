#!/bin/bash
#SBATCH --partition=compute         # Queue selection
#SBATCH --job-name=EnTAP      # Job name
#SBATCH --mail-type=ALL             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ztobias@whoi.edu  # Where to send mail
#SBATCH --ntasks=35                  # Run on a single CPU
#SBATCH --mem=186gb                   # Job memory request
#SBATCH --time=23:59:59             # Time limit hrs:min:sec
#SBATCH --output=EnTAP_%j.log  # Standard output/error

pwd; hostname; date

export OMP_NUM_THREADS=35
 
module load anaconda 

source ~/.bash_profile
source activate EnTAP

EnTAP --config -d /vortexfs1/scratch/ztobias/RhithroLoxo_DE/db/nr.gz \
    -d /vortexfs1/scratch/ztobias/RhithroLoxo_DE/db/uniprot_trembl.fasta.gz \
    -d /vortexfs1/scratch/ztobias/RhithroLoxo_DE/db/refseq_complete.faa.gz \
    -d /vortexfs1/scratch/ztobias/RhithroLoxo_DE/db/uniprot_sprot.fasta.gz \
    -d /vortexfs1/scratch/ztobias/RhithroLoxo_DE/db/uniref90.fasta.gz \
    -t 35

date
