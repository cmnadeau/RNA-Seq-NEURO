# RNA-seq analysis pipeline

[![Snakemake](https://img.shields.io/badge/snakemake==5.25.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io) [![miniconda](https://img.shields.io/badge/install%20with-conda-green.svg)](https://docs.conda.io/en/latest/miniconda.html)

This is an altered [Snakemake](https://snakemake.readthedocs.io/en/stable/) based pipeline for RNA-seq originally in the [Tumor Genome Core Analysis](http://www.tgac.nl/) housed in the [Cancer Center Amsterdam](https://www.vumc.com/departments/cancer-center-amsterdam.htm), at [Amsterdam UMC location VUmc](https://www.vumc.nl/) and part of the Department of Pathology for use with RNA-Seq data for Alzheimer's and Parkinson's disease.

The pipeline processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)), aligns the reads ([HISAT2](http://daehwankimlab.github.io/hisat2/)), generates gene counts ([HTSeq/DESeq2](https://htseq.readthedocs.io/en/master/)) and performs quality-control on the results ([MultiQC](https://multiqc.info/)). Paired-end (PE) and single read (SR) are supported.


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
git clone https://github.com/cmnadeau/RNA-Seq-NEURO
conda env create --name RNAseq --file env.yaml
```


## Path Configuration & Running the pipeline

Before attempting to run the pipeline, please open *config.yaml*. Inside, you will encounter **Path Configuration** and **Software Options**.

1. On **Path configuration**, first, you have to choose whether your data is PE or SR and after change the fastq path to the path where your fastq files are actually stored.
2. On **Software Options**, you will find several options that can be modified by the user. Please, have a look at it before running the pipeline.

All the software used in the pipeline is installed by conda or executed in a wrapper. We recommend to run the pipeline from a different location than the pipeline path, like the example below:

```
snakemake -s PATH_TO_PIPELINE/Snakefile --use-conda --cores=24
```
With --use-conda option, the pipeline will create environments to run rules based on *env.yaml*.
**Note** the pipeline assumes that *config.yaml* is available at the location where the pipeline is executed.
