#!/bin/bash
#SBATCH --partition=scavenger         # Queue selection
#SBATCH --qos=scavenger
#SBATCH --requeue
#SBATCH --job-name=EnTAP      # Job name
#SBATCH --mail-type=ALL             # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=ztobias@whoi.edu  # Where to send mail
#SBATCH --ntasks=35                  # Run on a single CPU
#SBATCH --mem=186gb                   # Job memory request
#SBATCH --time=48:00:00             # Time limit hrs:min:sec
#SBATCH --output=EnTAP_%j.log  # Standard output/error

pwd; hostname; date

export OMP_NUM_THREADS=35
 
module load anaconda 

source ~/.bash_profile
source activate EnTAP

EnTAP --runN -i /vortexfs1/scratch/ztobias/RhithroLoxo_DE/txms/rhithro/rhithro_txm_long.fasta \
    -d entap_outfiles/bin/nr.dmnd \
    -d entap_outfiles/bin/refseq_complete.dmnd \
    -d entap_outfiles/bin/uniprot_sprot.dmnd \
    -d entap_outfiles/bin/uniprot_trembl.dmnd \
    -d entap_outfiles/bin/uniref90.dmnd \
    -t 35 \
    -c bacteria \
    -c archaea \
    -c viruses \
    -c platyhelminthes \
    -c nematoda \
    -c fungi \
    -c alveolata \
    -c viridiplantae \
    -c rhodophyta \
    -c amoebozoa \
    -c rhizaria \
    -c stramenopiles \
    -c rhizocephala \
    -c entoniscidae \
    --taxon crustacea

date
