## Differential expression analysis of _Rhithropanopeus harrisi_'s evolved response to _Loxothylacus panopaei_ infection

#### by Zachary Tobias

This GitHub repository contains all necessary scripts and metadata for the execution of a differential expression analysis pipeline of an RNA-seq dataset from experimental infections of various populations of the mud crab _Rhithropanopeus harrisi_ with differing degrees of historical exposure to the parasitic barnacle _Loxothylacus panopaei_. This repository serves as a supplement to my first-year research report entilted "TITLE HERE!". 

The analysis starts from previously demultiplexed, trimmed, and cleaned 50bp, single-end sequence reads from an Illumina HiSeq 2000. It also utilizes a previously generated transcriptome for the parasite. These files can be made available upon request. Aside from these two exceptions, all steps of the pipeline are included herein and should be readily repeatable by readers familiar with basic programming in bash, R, and python. Most commands up through the generation of the differential expression matrix are executed using the workflow engine 'Snakemake', relying on packages installed in conda environments and a few borrowed or custom scripts. The major exception to this are those commands with the functional annotation package EnTAP v0.9.0-beta, which lacks a conda distribution and must be installed manually. Instructions for installation are included below in the section below entitled "EnTAP Setup"; additional instructions can be found [here](https://entap.readthedocs.io/en/latest/introduction.html). 

All steps in the pipeline were carried out on Poseidon, the high performance cluster (HPC) at Woods Hole Oceanographic Institution. 

### Reproducibility statement

This analysis was performed with ease of reproducibility in mind, both for my own sake in organizing the steps of the pipeline and for those interested in understanding the nuts and bolts or even replicating the results independently. That being said, there are some steps that may be influences by idiosyncracies of the configuration of the Poseidon HPC, specific versions of software, or release dates of reference sequence databases. 

### Overview of steps


### Executing the `Snakemake` pipeline



### Running `DESeq2` and `WCGNA`  from within a jupyter notebook on HPC compute node

Running the downstream differential expression and associated analyses interactively is great for fine-tuning code and exploring the data. As such, I chose to perform this analysis from within jupyter notebooks instead of writing .R scripts to be executed within the Snakemake pipeline. DESeq2 is performed in a notebook  called `DESeq2_RhithroLoxo.ipynb` and WGCNA in two notebooks, one for contrasting module expression with infection status and sex, `WGCNA_FP.ipynb`, and another for contrasting module expression with infection status and range, `WGCNA_noFP.ipynb` . These can be found in the `jupyter_notebooks/` directory within this repo. To harness the computational power of Poseidon (more RAM, parallelization, etc.), I launched the jupyter notebooks from within an interactive session on a compute node, instead of from one of the two login nodes. All of the notebooks operate within the `deseq2` conda environment in the `envs/` directory. Big thanks to Harriet Alexander for providing the instructions for this on her [blog](https://alexanderlabwhoi.github.io/post/2019-03-08_jpn_slurm/)!

From within the main directory, lauch an interactive session on a compute node and activate the `deseq2` environment. (If you haven't already, use the `deseq2.yaml` file provided in `envs/` directory for creating the deseq2 conda environment within your home directory on the cluster, i.e.`conda env create -f envs/deseq2.yaml`.)

```
srun -p compute --time=04:00:00 --ntasks-per-node 12 --mem 180gb --pty bash
conda activate deseq2
export XDG_RUNTIME_DIR=""
jupyter notebook --no-browser --port 8888
```

This will lauch a jupyter notebook at a port listed at `http://localhost:8888/`. It will use a different port number if that one is already in use. Remember this number! You also need to take note of the name of the node you are running on. It should be in your prompt, i.e. if your prompt looks like this: `(deseq2) [ztobias@pn039 RhithroLoxo_DE]$`, then the node is pn039.

On your home computer, add the following function to your `.bash_profile` to streamline opening the port.

```
function jptnode(){
    # Forwards port $1 from node $2 into port $1 on the local machine and listens to it
        ssh -t -t ztobias@poseidon.whoi.edu -L $1:localhost:$1 ssh $2 -L $1:localhost:$1
}
```

Don't forget to source your `.bash_profile`! (Type `source ~/.bash_profile`, or just close and reopen a terminal window.)

Now you can get the notebook running in a browser window on your local machine! Here you'll use the port and node numbers from above.

```
jptnode 8888 pn039
```

If you get an error saying it can't listen because the port is busy, try starting the jupyter notebook specifying the port, i.e. `jupyter notebook --no-browser --port XXXX`, using a different number. Sometimes the ports are busy. I think they take a while to actually close if you're doing this repeatedly.

Then go to your preferred web browser and type `localhost:8888` in the search bar (or whatever port number was assigned above). If you have configured your password, it will ask you for it. If not, you will have to set it on Poseidon by typing `jupyter notebook password`.

Another thing to be aware of. The `deseq2` conda environment does not include the `DESeq2` conda distribution. It has a lot of package conflicts. Instead, from within the `deseq2` environment, launch R and download `DESeq2` and `WGCNA` using `biocmanager` or `install.packages()` from base R. This only has to be done once. It will take a while and is quite verbose. Also download a host of other R packages, which are either dependencies of one of the three main analysis packages or will be useful for plotting, etc. If it asks you to update packages, JUST SAY NO! The environment is already set up as we want it; no need to go muck it up.

```
R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
install.packages("pheatmap")
install.packages("VennDiagram")
install.packages("ashr")
BiocManager::install("WGCNA")
install.packages("gdata")
install.packages("UpSetR")
install.packages("flashClust")
```

Okay now you're all set to actually run the DESeq2 and WGCNA analyses from the jupyter notebook! 

### Running `GO_MWU`

`GO_MWU` is not available in a CRAN or BioConductor install. Instead you just pull the scripts directly from the GitHub [repository](https://github.com/z0on/GO_MWU). 

In the main folder, clone the repository. 

````
git clone https://github.com/z0on/GO_MWU.git
````

This will have a bunch of other files, but all you need are: `GO_MWU.R`, `gomwu_a.pl`, `gomwu_b.pl`, and `gomwu.functions.R`. Feel free to remove others or use them as examples. You should replace the `go.obo` with the most recent release. (Note that the link in the GO_MWU README.md is no longer functional. Use link in code block below.)

```
rm go.obo
wget http://current.geneontology.org/ontology/go.obo
```

Then add all the necessary files, copy code blocks from the `GO_MWU.R` script into a notebook entitled `GO_MWU.ipynb` in that same folder (not the main `juptyer_notebooks/` directory), edit as [instructed](https://github.com/z0on/GO_MWU/blob/master/README.md), and run! Note that this notebook has to be launched from within the `deseq2` conda environment, or you'll lack the R kernel and won't be able to run it. 

### `EnTAP` setup

Navigate to the main repository directory and activate the `EnTAP` conda environment (previously loaded from `EnTAP.yaml` file in `envs/`). This environment contains all dependencies for installing/compiling. 

```
conda activate EnTAP
```

Install EnTAP and compile according to [instructions](https://entap.readthedocs.io/en/latest/installation.html). 

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

Edit bash profile and source. Add `$SCRATCH/RhithroLoxo_DE/EnTAP` to path

```
source ~/.bash_profile
conda activate EnTAP
```

Copy the `configure_EnTAP.sh` and `run_EnTAP.sh` scripts from `scripts/` to `EnTAP/`. Then replace the existing `entap_config.txt` file with the one provided in `metadata/`. Then run `configure_EnTAP.sh`

```
cp ../scripts/*EnTAP.sh ./
rm entap_config.txt
cp ../metadata/entap_config.txt ./
sbatch configure_EnTAP.sh
```

Check out the log files and output to make sure the databases were downloaded and indexed.

Then run the run_EnTAP.sh script

```
sbatch `run_EnTAP.sh`
```

This creates a myriad of output files in the `EnTAP/entap_outfiles` directory. The eggNOG mapping results are reformatted for input into `GO_MWU` using the script `EnTAP2GO.py` in the `scripts/` folder.

RUNNING THE ENTAP2GO.PY SCRIPT HERE!

```

```