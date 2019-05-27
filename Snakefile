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
PATH_STAR = config['path']['star']

print(PATH_FASTQ)

Files = []
for p in PATH_FASTQ:
    Files = Files + glob.glob(path.join(p, '*.fastq.gz'))
    print(Files)
#   ignore = filter(lambda file: not file.endswith('_trimmed.fastq.gz'), p)
#   filter(lambda Files: not file.endswith('_trimmed.fastq.gz'), Files)


RNAIDs = [f.split('/')[-1].split('.')[0] for f in Files]
IDtoPath = dict()

for f in Files:
    IDtoPath[f.split('/')[-1].split('.')[0]] = '/'.join(f.split('/')[:-1]) + '/'

print(RNAIDs)

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
        PATH_QC + 'multiqc.html'
     
rule index:
        input:
            fa = config['fasta'], # provide your reference FASTA file
            gtf = config['gtf'] # provide your GTF file
        output: PATH_STAR # you can also rename the index folder
        threads: 20 # set the maximum number of available cores
        shell:"""
            STAR --runThreadN {threads} 
            --runMode genomeGenerate 
            --genomeDir {output} 
            --genomeFastaFiles {input.fa} 
            --sjdbGTFfile {input.gtf} 
            --sjdbOverhang 100
            # cleanup
            rm -rf _STARtmp
            mkdir -p logs/star
            mv Log.out logs/star/star_index.log
            """        

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
        annot = PATH_STAR + '{star}+_.gtf'
    shell:
        """
        featureCounts -p -g gene_id -a {params.annot} -o {output} -T {threads} {input} 2> {log}
        """

rule trimmomatic:
    input:
        '{path}/{fastq}.fastq.gz'
    output:
        temp('{path}/{fastq,[a-zA-Z0-9-_]+}_trimmed.fastq.gz')
    params:
       config['trim']['params']['trimmer']
    threads:
        config['trim']['threads']
    log:
        PATH_LOG + '{fastq}.log'
    wrapper:
        "0.34.0/bio/trimmomatic/se" # Trim single-end reads

rule star_alignment:
    input:
       fq1 = lambda wildcards: IDtoPath[wildcards.sample] + wildcards.sample + '_trimmed.fastq.gz'
    output:
        PATH_BAM + '{sample}/Aligned.out.bam'
    log:
        PATH_LOG + '{sample}_star.log'
    params:
        index = PATH_STAR 
    threads:
        config['star']['threads']
    wrapper:
        "0.34.0/bio/star/align" # Map SE reads with STAR

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
        lambda wildcards: IDtoPath[wildcards.sample] +wildcards.sample + '_trimmed.fastq.gz'
    output:
        html=PATH_QC+"{sample}_fastqc.html",
        zip=PATH_QC+"{sample}_fastqc.zip"
    wrapper:
        "0.34.0/bio/fastqc" # generates fastq qc statistics; quality control tool for high throughput sequence data

rule multiqc:
    input:
        expand(PATH_QC+"{sample}_fastqc.zip", sample=RNAIDs)
    output:
        PATH_QC+'multiqc.html'
    log:
        PATH_LOG+'multiqc.log'
    wrapper:
        '0.34.0/bio/multiqc' # modular tool to aggregate results from bioinformatics analyses across many samples into a single report
