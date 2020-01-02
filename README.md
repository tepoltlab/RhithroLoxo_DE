## Differential expression analysis of _Rhithropanopeus harrisi_'s evolved response to _Loxothylacus panopaei_ infection

#### by Zachary Tobias

This GitHub repository contains all necessary scripts and metadata for the execution of a differential expression analysis pipeline of an RNA-seq dataset from experimental infections of various populations of the mud crab _Rhithropanopeus harrisi_ with differing degrees of historical exposure to the parasitic barnacle _Loxothylacus panopaei_. This repository serves as a supplement to my first-year research report entilted "TITLE HERE!". 

The analysis starts from previously demultiplexed, trimmed, and cleaned 50bp, single-end sequence reads from an Illumina HiSeq 2000. It also utilizes a previously generated transcriptome for the parasite. These files can be made available upon request. Aside from these two exceptions, all steps of the pipeline are included herein and should be readily repeatable by readers familiar with basic programming in bash, R, and python. Most commands up through the generation of the differential expression matrix are executed using the workflow engine 'Snakemake', relying on packages installed in conda environments and a few borrowed or custom scripts. The major exception to this are those commands with the functional annotation package EnTAP v0.9.0-beta, which lacks a conda distribution and must be installed manually. Instructions for installation are included below in the section below entitled "EnTAP Setup"; additional instructions can be found at <https://entap.readthedocs.io/en/latest/introduction.html>. 

All steps in the pipeline were carried out on Poseidon, the high performance cluster (HPC) at Woods Hole Oceanographic Institution. 

### Reproducibility statement



### How to execute Snakemake pipeline



### How to run DESeq2 analysis in a jupyter notebook on interactive node

Running the DESeq2 differential expression analysis interactively is great for fine-tuning code and exploring the data. As such, I chose to perform this analysis from within a jupyter notebook instead of writing a .R script to be executed within the Snakemake pipeline. This notebook is called `DESeq2_RhithroLoxo.ipynb` and can be found in the `jupyter_notebooks/` directory within this repo. To harness the computational power of Poseidon (more RAM, multithreading, etc.), I launched a jupyter notebook from within an interactive session on a compute node, instead of from one of the two login nodes. Big thanks to Harriet Alexander for providing the instructions on her blog! <https://alexanderlabwhoi.github.io/post/2019-03-08_jpn_slurm/> I use a slightly modified approach, outlined below.

From within the main directory, lauch an interactive session on a compute node and activate the `deseq2` environment. (If you haven't already, use the `deseq2.yaml` file provided in `envs/` directory for creating the deseq2 conda environment within your home directory on the cluster, i.e.`conda env create -f envs/deseq2.yaml`.)

```
srun -p compute --time=04:00:00 --ntasks-per-node 8 --mem 40gb --pty bash
conda activate deseq2
export XDG_RUNTIME_DIR=""
jupyter notebook --no-browser
```

This will lauch a jupyter notebook at a port listed at `http://localhost:8888/`. The number following 'localhost:' may vary. Remember or copy this number! You also need to take note of the name of the node you are running on. It should be in your prompt, i.e. if your prompt looks like this: `(deseq2) [ztobias@pn039 RhithroLoxo_DE]$`, then the node is pn039.

On your home computer, add the following function to your `.bash_profile` to streamline opening the port.

```
function jptnode(){
    # Forwards port $1 from node $2 into port $1 on the local machine and listens to it
        ssh -t -t ztobias@poseidon.whoi.edu -L $1:localhost:$1 ssh $2 -L $1:localhost:$1
        open -a "/Applications/Google Chrome.app" "http://localhost:$1"
}
```

Don't forget to source your `.bash_profile`! (Type `source ~/.bash_profile`, or just close and reopen a terminal window.)

Now you can get the notebook running in a browser window on your local machine! Here you'll use the port and node numbers from above.

```
jptnode 8888 pn039
```

If you get an error saying it can't listen because the port is busy, try starting the jupyter notebook specifying the port, i.e. `jupyter notebook --no-browser --port XXXX`, using a different number. Sometimes the ports are busy. I think they take a while to actually close if you're doing this repeatedly. Who knows. 

Also, for some reason the browser might not launch automatically. Just type `localhost:XXXX` in the browser address bar, where XXXX is the port number, and you should be good to go. Also, when you type `jptnode ...`, it seems to actually log you into that node in the Terminal. Not sure why. Just minimize the window and ignore.

Another thing to be aware of. The `deseq2` conda environment does not include the DESeq2 conda distribution. It has a lot of package conflicts. Instead, from within the `deseq2` environment, launch R and download DESeq2 using `biocmanager`. This only has to be done once. It will take a while and is quite verbose. Also download the R packages `apeglm`,`pheatmap`, and `VennDiagram`, which are either dependencies of DESeq2 or will be useful for plotting, etc. If it asks you to update packages, JUST SAY NO! The environment is already set up as we want it; no need to go muck it up.

```
R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
install.packages("pheatmap")
install.packages("VennDiagram")
```

Okay now you're all set to actually run the DESeq2 analysis from the jupyter notebook! 
 

### EnTAP setup

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

This creates a myriad of output files in the `EnTAP/entap_outfiles` directory. The results are parsed in an accompanying jupyter notebook entitled `parse_annot.ipynb`.
