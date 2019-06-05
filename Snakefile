import os.path as path
import pandas as pd
import os, sys, glob
import numpy as np

configfile: 'config.yaml'

PATH_FASTQ = config['path']['fastq']
PATH_LOG = config['path']['log']
PATH_QC = config['path']['qc']
PATH_BAM = config['path']['bam']
PATH_OUT = config['path']['out']

PATH_STARINDEX = config['star']['index']
PATH_FASTA = config['star']['fasta']
PATH_GTF = config['star']['gtf']

PLATFORM = config['platform']['SEorPE'] # SE or PE
PREFIX = config['platform']['prefix']

Files = []
RNAIDs = []
for p in PATH_FASTQ:
    if PLATFORM in ['SE', 'se']:
        for prefix in PREFIX:
            NewFile = glob.glob(path.join(p, '*'+prefix+'.fastq.gz'))
            RNAIDs = RNAIDs + [f.split('/')[-1].split(prefix)[0] for f in NewFile]
            Files = Files + NewFile
    if PLATFORM in ['PE', 'pe']:
        R1 = glob.glob(path.join(p, '*'+PREFIX[0]+'.fastq.gz'))
        R2 = glob.glob(path.join(p, '*'+PREFIX[1]+'.fastq.gz'))
        ID_R1 = [f.split('/')[-1].split(PREFIX[0])[0] for f in R1]
        ID_R2 = [f.split('/')[-1].split(PREFIX[1])[0] for f in R2]
        if not set(ID_R1) == set(ID_R2):
            Mismatch = list(set(ID_R1)-set(ID_R2)) + list(set(ID_R2)-set(ID_R1))
            raise Exception("Missing pairs: " + ' '.join(Mismatch))
        RNAIDs = RNAIDs + ID_R1
        Files = Files + R1 + R2


def ID2Fasta(ID, PLATFORM):
    if PLATFORM in ['SE', 'se']:
        return [s for s in Files if ID in s] 
    if PLATFORM in ['PE', 'pe']:
        return [s for s in Files if ID in s and PREFIX[1] in s] + [s for s in Files if ID in s and PREFIX[2] in s]
        

def ID2FastqPath(ID):
    return '/'.join([s for s in Files if ID in s][0].split('/')[:-1])

print(RNAIDs)
#print(Files)

def normalize_counts(counts):
    """Normalizes expression counts using DESeq's median-of-ratios approach."""

    with np.errstate(divide="ignore"):
        size_factors = estimate_size_factors(counts)
        return counts / size_factors


def estimate_size_factors(counts):
    """Calculate size factors for DESeq's median-of-ratios normalization."""

    def _estimate_size_factors_col(counts, log_geo_means):
        log_counts = np.log(counts)
        mask = np.isfinite(log_geo_means) & (counts > 0)
        return np.exp(np.median((log_counts - log_geo_means)[mask]))

    log_geo_means = np.mean(np.log(counts), axis=1)
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col, axis=0,
        arr=counts, log_geo_means=log_geo_means)

    return size_factors

rule all:
    input:
        PATH_OUT + 'featurecounts.log2.txt',
        PATH_OUT + 'multiqc.html',
        PATH_OUT + 'compress_fastq.zip',
        PATH_OUT + 'compress_bam.zip'

rule zip_fastq:
    input:
        Files
    output:
        PATH_OUT + 'compress_fastq.zip'
    shell:
        """
        zip {output} {input} 
        """
        

rule zip_bams:
    input:
        expand(PATH_BAM + '{sample}.mq20.bam', sample=RNAIDs)
    output:
        PATH_OUT + 'compress_bam.zip'
    shell:
        """
        zip {output} {input} 
        """

     
#rule index:
#        input:
#            fa = config['fa_star'], # provide your reference FASTA file
#            gtf = config['gtf_star'] # provide your GTF file
#        output: PATH_STARINDEX
#        threads: 20 
#        shell:"""
#            STAR --runThreadN {threads} 
#            --runMode genomeGenerate 
#            --genomeDir {output} 
#            --genomeFastaFiles {input.fa} 
#            --sjdbGTFfile {input.gtf} 
#            --sjdbOverhang 100
#            # cleanup
#            rm -rf _STARtmp
#            mkdir -p log/star
#            mv Log.out log/star/star_index.log
#            """        

rule normalize_counts:
    input:
        PATH_OUT + "featurecounts.txt"
    output:
        PATH_OUT + "featurecounts.log2.txt"
    run:
        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)), skiprows=1)
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_counts.to_csv(output[0], sep="\t", index=True)

rule featurecount:
    input:
        expand(PATH_BAM + '{sample}.mq20.bam', sample=RNAIDs)
    output:
        PATH_OUT + 'featurecounts.txt'
    log:
        PATH_LOG + 'featurecount_log'
    threads: config['ftcount']['threads']
    params:
        annot = PATH_GTF 
    shell:
        """
        featureCounts -p -g gene_id -a {params.annot} -o {output} -T {threads} {input} 2> {log}
        """

#if PLATFORM in ['SE', 'se']:
#    rule trimmomatic:
#        input:
#            '{path}/{fastq}.fastq.gz'
#        output:
#            temp('{path}/{fastq,[a-zA-Z0-9-_]+}_trimmed.fastq.gz')
#        params:
#            config['trim']['params']['trimmer']
#        threads:
#            config['trim']['threads']
#        log:
#            PATH_LOG + '{fastq}.log'
#        wrapper:
#            "0.34.0/bio/trimmomatic/se" # Trim single-end reads

#    rule star_alignment:
#        input:
#           fq1 = lambda wildcards: IDtoPath[wildcards.sample] + wildcards.sample + '_trimmed.fastq.gz'
#        output:
#            PATH_BAM + '{sample}/Aligned.out.bam'
#        log:
#            PATH_LOG + '{sample}_star.log'
#        params:
#            index = PATH_STARINDEX
#        threads:
#            config['star']['threads']
#        wrapper:
#            "0.34.0/bio/star/align" # Map SE reads with STAR

if PLATFORM in ['PE', 'pe']:
    rule trimmomatic:
        input:
            r1 = lambda wildcards: path.join(wildcards.path, wildcards.sample + PREFIX[0] + '.fastq.gz'),
            r2 = lambda wildcards: path.join(wildcards.path, wildcards.sample + PREFIX[1] + '.fastq.gz')
        output:
            r1= temp(path.join('{path}', '{sample}' + PREFIX[0] + '.trimmed.fastq.gz')),
            r2= temp(path.join('{path}', '{sample}' + PREFIX[1] + '.trimmed.fastq.gz')),
            r1_unpaired= temp(path.join('{path}', '{sample}' + PREFIX[0] + '.trimmed_unpaired.fastq.gz')),
            r2_unpaired= temp(path.join('{path}', '{sample}' + PREFIX[1] + '.trimmed_unpaired.fastq.gz'))
        params:
            trimmer = config['trim']['params']['trimmer']
        threads:
            config['trim']['threads']
        log:
            PATH_LOG + '{sample}.trimmomatic.log'
        wrapper:
            "0.34.0/bio/trimmomatic/pe" # Trim single-end reads

    rule star_alignment:
        input:
           fq1 = lambda wildcards: path.join(ID2FastqPath(wildcards.sample), wildcards.sample + PREFIX[0] + '.trimmed.fastq.gz'),
           fq2 = lambda wildcards: path.join(ID2FastqPath(wildcards.sample), wildcards.sample + PREFIX[1] + '.trimmed.fastq.gz')
        output:
            temp(PATH_BAM + '{sample}/Aligned.out.bam')
        log:
            PATH_LOG + '{sample}_star.log'
        params:
            index = PATH_STARINDEX 
        threads:
            config['star']['threads']
        wrapper:
            "0.34.0/bio/star/align" # Map PE reads with STAR

   
rule filter_bam:
    input:
        PATH_BAM+'{sample}/Aligned.out.bam'
    output:
        PATH_BAM+'{sample}.mq20.bam'
    params:
        config['bam']['params']
    wrapper:
        '0.34.0/bio/samtools/view' # convert or filter SAM/BAM;

rule fastqc_fastq:
    input:
        PATH_BAM + '{sample}.mq20.bam'
    output:
        html=PATH_QC+"{sample}_fastqc.html",
        zip=PATH_QC+"{sample}_fastqc.zip"
    wrapper:
        "0.34.0/bio/fastqc" # generates fastq qc statistics; quality control tool for high throughput sequence data

rule multiqc:
    input:
        expand(PATH_QC+"{sample}_fastqc.zip", sample=RNAIDs)
    output:
        PATH_OUT+'multiqc.html'
    log:
        PATH_LOG+'multiqc.log'
    wrapper:
        '0.35.0/bio/multiqc' # modular tool to aggregate results from bioinformatics analyses across many samples into a single report
