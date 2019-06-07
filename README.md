# RNA-seq analysis pipeline 

This is a snakemake based pipeline for RNA-seq used in the [Tumor Genome Core Analysis](http://www.tgac.nl/) housed in the [Cancer Center Amsterdam](https://www.vumc.com/departments/cancer-center-amsterdam.htm), at [Amsterdam UMC location VUmc](https://www.vumc.nl/) and part of the Department of Pathology.

## Instalation

The pipeline is preliminary used in linux environment with conda/singularity available.

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

### Step 2: Downloading repository, creating environment & singularity image:

```
mkdir snakemake_RNAseq
cd snakemake_RNAseq
git clone https://github.com/tgac-vumc/RNA-seq
conda env create --name snakemake --file env.yaml
singularity pull shub://tgac-vumc/RNA-seq
```
All the softwares used in the pipeline are installed by conda or executed in wrapper. We recommend to run the pipeline from a different location than the pipeline path, like the example below:

```
snakemake -s PATH_TO_PIPELINE/Snakefile --use-conda --use-singularity --cores=24
```

With --use-conda option, the pipeline will create environments to run rules based on env.yaml. The singularity container holds a virtual environment of CentOS 7.
Note that the pipeline assumes that config.yaml is available at the location where the pipeline is executed.


