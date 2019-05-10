import os.path as path
import pandas as pd
import os, sys, glob
import numpy as np

PATH_FASTQ = ['fastq/']

Files = []
for p in PATH_FASTQ:
    Files = Files + glob.glob(p + '*.fastq.gz')
#   ignore = filter(lambda file: not file.endswith('_trimmed.fastq.gz'), p)
#   filter(lambda Files: not file.endswith('_trimmed.fastq.gz'), Files)


RNAIDs = [f.split('/')[-1].split('.')[0] for f in Files]
IDtoPath = dict()
for f in Files:
    IDtoPath[f.split('/')[-1].split('.')[0]] = '/'.join(f.split('/')[:-1]) + '/'

PATH_LOG =  'log/'
PATH_QC = 'data/RNA_seq/qc/'
PATH_BAM = 'data/RNA_seq/bams/'
PATH_OUT = 'data/RNA_seq/output/'

print(RNAIDs)

rule all:
    input:
        PATH_OUT + 'featurecounts.log2.txt',
        PATH_QC + 'multiqc.html'

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
    threads:20
    params:
        annot = '/home/y.kim1/Resource/gtf/Homo_sapiens.GRCh37.75.gtf'
    shell:
        """
        featureCounts -p -g gene_id -a {params.annot} -o {output} -T {threads} {input} 2> {log}
        """

rule trimmomatic:
    input:
        '{path}/{fastq}.fastq.gz'
    output:
        '{path}/{fastq,[a-zA-Z0-9-_]+}_trimmed.fastq.gz',
    params:
        trimmer = ["TRAILING:3"] # Cut bases off the end of a read, if below a threshold quality, use LEADING for begining
    threads: 4
#    log:
#        PATH_LOG + '{fastq}.log'
    wrapper:
        "0.31.1/bio/trimmomatic/se" # Trim single-end reads

rule star_alignment:
    input:
       fq1 = lambda wildcards: IDtoPath[wildcards.sample] + wildcards.sample + '_trimmed.fastq.gz'
    output:
        PATH_BAM + '{sample}/Aligned.out.bam'
    log:
        PATH_LOG + '{sample}_star.log'
    params:
        index='/home/y.kim1/Resource/STAR_index/GRCh37.75/'
    threads:15
    wrapper:
        "0.34.0/bio/star/align" # v = 0.32.0 /

rule filter_bam:
    input:
        PATH_BAM+'{sample}/Aligned.out.bam'
    output:
        PATH_BAM+'{sample}.mq20.bam'
    params:
        '-b -h -q20' #-qINT Skip alignments with MAPQ smaller than INT [0].
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
