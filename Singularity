Bootstrap: yum
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-7/7/os/x86_64/
Include: yum

%environment    
    export PATH=/usr/local/bin:$PATH

%post
    ./environemnt

    apt-get -y update
    apt-get -qq -y install curl
    curl -sSL https://repo.continuum.io/archive/Anaconda2-5.0.1-Linux-x86_64.sh -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -bfp /usr/local
    rm -rf /tmp/miniconda.sh
    
    git clone https://github.com/barbarabarbosa/RNA-seq
    conda env create -n RNA-seq -f RNA-seq/env.yaml
