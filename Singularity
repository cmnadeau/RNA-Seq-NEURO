Bootstrap: yum
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-7/7/os/x86_64/
Include: yum

%environment    
    export PATH=/usr/local/bin:$PATH

%post
    ./environment

    yum -y update
    yum -qq -y install curl tar bzip2 git zip
    curl -sSL https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -bfp /usr/local
    
    git clone https://github.com/tgac-vumc/RNA-seq
    conda install -c conda-forge mamba -y 
    mamba create -c conda-forge -c bioconda -n RNA-seq snakemake==5.25.0 subread -y


    rm -rf /tmp/miniconda.sh
