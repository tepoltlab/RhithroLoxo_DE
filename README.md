## Differential expression analysis of _Rhithropanopeus harrisi_'s evolved response to _Loxothylacus panopaei_ infection

#### by Zachary Tobias

This GitHub repository contains all necessary scripts and metadata for the execution of a differential expression analysis pipeline of an RNA-seq dataset from experimental infections of various populations of the mud crab _Rhithropanopeus harrisi_ with differing degrees of historical exposure to the parasitic barnacle _Loxothylacus panopaei_. This repository serves as a supplement to my first-year research report entilted "TITLE HERE!". 

The analysis starts from previously demultiplexed, trimmed, and cleaned 50bp, single-end sequence reads from an Illumina HiSeq 2000. It also utilizes a previously generated transcriptome for the parasite. These files can be made available upon request. Aside from these two exceptions, all steps of the pipeline are included herein and should be readily repeatable by readers familiar with basic programming in bash, R, and python. Most commands up through the generation of the differential expression matrix are executed using the workflow engine 'Snakemake', relying on packages installed in conda environments and a few borrowed or custom scripts. The major exception to this are those commands with the functional annotation package EnTAP v0.9.0-beta, which lacks a conda distribution and must be installed manually. Instructions for installation are included below in the section entitled "EnTAP Setup"; additional instructions can be found at <https://entap.readthedocs.io/en/latest/introduction.html>. 

### EnTAP Setup

Start interactive session in RhithroLoxo_DE directory

```
srun -p scavenger --time=04:00:00 --ntasks-per-node=36 --mem=100gb --pty bash
conda activate EnTAP
```

Install EnTAP and compile according to instructions 

```
git clone https://gitlab.com/enTAP/EnTAP.git
cd EnTAP
cd libs/diamond-0.9.9
mkdir bin
cd bin
cmake ..
make #don’t do global install
#skip RSEM install because not expression filtering. throws an error anyway, because it's looking for usr/bin/env as if it were on local computer
cd ../../../
cmake CMakeLists.txt
make
```

Edit bash profile and source. Add $SCRATCH/RhithroLoxo_DE/EnTAP to path

```
source ~/.bash_profile
conda activate EnTAP
```

Copy the configure_EnTAP.sh and run_EnTAP.sh scripts from scripts/ to EnTAP/, and run configure_EnTAP.sh

```
cp ../scripts/*EnTAP.sh ./
sbatch configure_EnTAP.sh
```

Check out the log files and output to make sure the databases were downloaded and indexed.

Then run the run_EnTAP.sh script

```
sbatch run_EnTAP.sh
```


