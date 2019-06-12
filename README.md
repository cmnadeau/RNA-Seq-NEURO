# RNA-seq analysis pipeline 

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.13.3-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io) [![singularity-hub](https://img.shields.io/badge/install%20with-singularity--hub-red.svg)](https://singularity-hub.org/collections/3066) [![miniconda](https://img.shields.io/badge/install%20with-conda-green.svg)](https://docs.conda.io/en/latest/miniconda.html)

This is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) based pipeline for RNA-seq used in the [Tumor Genome Core Analysis](http://www.tgac.nl/) housed in the [Cancer Center Amsterdam](https://www.vumc.com/departments/cancer-center-amsterdam.htm), at [Amsterdam UMC location VUmc](https://www.vumc.nl/) and part of the Department of Pathology.

The pipeline processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), aligns the reads ([STAR](https://github.com/alexdobin/STAR)), generates gene counts ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/)) and performs quality-control on the results ([MultiQC](https://multiqc.info/)).


## Installation

The pipeline is preliminary used in linux environment with conda/singularity available.

### Using Conda
### Step 1: Installing Miniconda 3
First, please open a terminal or make sure you are logged into your Linux VM. Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with:

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
On MacOS X, download and install with:

```
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

### Step 2: Downloading repository & creating environment

```
mkdir snakemake_RNAseq
cd snakemake_RNAseq
git clone https://github.com/tgac-vumc/RNA-seq
conda env create --name RNAseq --file env.yaml
```
All the softwares used in the pipeline are installed by conda or executed in wrapper. We recommend to run the pipeline from a different location than the pipeline path, like the example below:

```
snakemake -s PATH_TO_PIPELINE/Snakefile --use-conda --cores=24
```
With --use-conda option, the pipeline will create environments to run rules based on *env.yaml*. **Note** the pipeline assumes that *config.yaml* is available at the location where the pipeline is executed. *config.yaml* is configurable by the user, please check it before running the pipeline.

### Using Singularity

The singularity container holds a virtual environment of CentOS 7 and it's available with:
```
singularity pull shub://tgac-vumc/RNA-seq
```
