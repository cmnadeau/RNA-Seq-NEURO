import os.path as path
#ewimport pandas as pd
import os, sys, glob
#import numpy as np

configfile: 'config.yaml'

PATH_FASTQ = config['path']['fastq']
PATH_TRIMMED = config['path']['trimmed']
PATH_LOG = config['path']['log']
PATH_QC = config['path']['qc']
PATH_BAM = config['path']['bam']
PATH_HTSEQ = config['path']['htseq_counts']
PATH_OUT = config['path']['out']

PATH_HISAT2INDEX = config['hisat2']['index']
PATH_GTF = config['hisat2']['gtf']
PATH_HTSEQ_GTF = config['htseq']['gtf']

PLATFORM = config['platform']['SRorPE'] # SE or PE
PREFIX = config['platform']['prefix']

Files = []
RNAIDs = []
for p in PATH_FASTQ:
    R1 = glob.glob(path.join(p, '*'+PREFIX[0]+'.fastq.gz'))
    R2 = glob.glob(path.join(p, '*'+PREFIX[1]+'.fastq.gz'))
    ID_R1 = [f.split('/')[-1].split(PREFIX[0])[0] for f in R1]
    ID_R2 = [f.split('/')[-1].split(PREFIX[1])[0] for f in R2]
    if not set(ID_R1) == set(ID_R2):
        Mismatch = list(set(ID_R1)-set(ID_R2)) + list(set(ID_R2)-set(ID_R1))
        raise Exception("Missing pairs: " + ' '.join(Mismatch))
    RNAIDs = RNAIDs + ID_R1
    Files = Files + R1 + R2

RNAIDs = list(set(RNAIDs))
print(RNAIDs)
def ID2TrimmedFastq(ID, EXT):
   return [re.sub(".fastq.gz", "", file) + EXT for file in Files if ID in file]

def ID2FastqPath(ID):
    return '/'.join([s for s in Files if ID in s][0].split('/')[:-1])

#def normalize_counts(counts):
#    """Normalizes expression counts using DESeq's median-of-ratios approach."""
#
#    with np.errstate(divide="ignore"):
#        size_factors = estimate_size_factors(counts)
#        return counts / size_factors


#def estimate_size_factors(counts):
#    """Calculate size factors for DESeq's median-of-ratios normalization."""
#
#    def _estimate_size_factors_col(counts, log_geo_means):
#        log_counts = np.log(counts)
#        mask = np.isfinite(log_geo_means) & (counts > 0)
#        return np.exp(np.median((log_counts - log_geo_means)[mask]))
#
#    log_geo_means = np.mean(np.log(counts), axis=1)
#    size_factors = np.apply_along_axis(
#        _estimate_size_factors_col, axis=0,
#        arr=counts, log_geo_means=log_geo_means)
#
#    return size_factors

rule all:
    input:
        #PATH_OUT + 'featurecounts.log2.txt',
#        PATH_OUT + 'multiqc.html',
#        PATH_OUT + 'multiqc_raw.html',
        PATH_OUT + 'compress_fastq.zip',
        PATH_OUT + 'compress_bam.zip',
        PATH_OUT + 'full.counts.txt'

rule zip_fastq:
    input:
        Files
    output:
        PATH_OUT + 'compress_fastq.zip'
    shell:
        """
        zip -j {output} {input}
        """


rule zip_bams:
    input:
        expand(PATH_BAM + '{sample}.sorted.out.bam', sample=RNAIDs)
    output:
        PATH_OUT + 'compress_bam.zip'
    shell:
        """
        zip -j {output} {input}
        """
#rule hisat2_index:
#    input:
#        fasta=config['hisat2']['fasta_index']
#    output:
#        directory(config['hisat2']['index'])
#    threads:20
#    log:
#        PATH_LOG + "hisat2_index.log"
#    wrapper:
#        "0.65.0/bio/hisat2/index"
##rule index:
#        input:
#            fa = config['fa_hisat2'], # provide your reference FASTA file
#            gtf = config['gtf_hisat2'] # provide your GTF file
#        output: PATH_HISAT2INDEX
#        threads: 20
#        shell:"""
#            HISAT2 --runThreadN {threads}
#            --runMode genomeGenerate
#            --genomeDir {output}
#            --genomeFastaFiles {input.fa}
#            --sjdbGTFfile {input.gtf}
#            --sjdbOverhang 100
#            # cleanup
#            rm -rf _HISAT2tmp
#            mkdir -p log/hisat2
#            mv Log.out log/hisat2/hisat2_index.log
#            """

#rule normalize_counts:#
#    input:
#        PATH_OUT + "featurecounts.txt"
#    output:
#        PATH_OUT + "featurecounts.log2.txt"
#    run:
#        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)), skiprows=1)#
#        norm_counts = np.log2(normalize_counts(counts) + 1)
#        norm_counts.to_csv(output[0], sep="\t", index=True)

#def feature_counts_extra(wildcards):
#    extra = config["ftcount"]["option"]
#    if PLATFORM in ['PE', 'pe']:
#        extra += " -p"
#    return extra


#rule featurecount:
#    input:
#        expand(PATH_BAM + '{sample}.mq20.bam', sample=RNAIDs)
#    output:
#        PATH_OUT + 'featurecounts.txt'
#    log:
#        PATH_LOG + 'featurecount_log'
##    params:
#        annot = PATH_GTF,
#        others = feature_counts_extra
#    shell:
#        """
#        featureCounts {params.others} -a {params.annot} -o {output} -T {threads} {input} 2> {log}
#        """

rule trimmomatic:
    input:
        r1 = lambda wildcards: path.join(ID2FastqPath(wildcards.sample), wildcards.sample + PREFIX[0] + '.fastq.gz'),
        r2 = lambda wildcards: path.join(ID2FastqPath(wildcards.sample), wildcards.sample + PREFIX[1] + '.fastq.gz')
    output:
        r1= path.join(PATH_TRIMMED, '{sample}' + PREFIX[0] + '.trimmed.fastq.gz'),
        r2= path.join(PATH_TRIMMED, '{sample}' + PREFIX[1] + '.trimmed.fastq.gz'),
        r1_unpaired= path.join(PATH_TRIMMED, '{sample}' + PREFIX[0] + '.trimmed_unpaired.fastq.gz'),
        r2_unpaired= path.join(PATH_TRIMMED, '{sample}' + PREFIX[1] + '.trimmed_unpaired.fastq.gz')
    params:
        trimmer = config['trim']['params']['trimmer']
    threads:
        config['trim']['threads']
    log:
        PATH_LOG + 'trimmomatic/{sample}.trimmomatic.log'
    shell:
        """
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} \
        {output.r2} {output.r2_unpaired} {params.trimmer} -threads {threads} > {log}
        """
#    wrapper:
#        "v1.3.2/bio/trimmomatic/pe" # Trim paired-end reads
rule unzip_trim:
    input:
        r1 = path.join(PATH_TRIMMED, '{sample}' + PREFIX[0] + '.trimmed.fastq.gz'),
        r2 = path.join(PATH_TRIMMED, '{sample}' + PREFIX[1] + '.trimmed.fastq.gz')
    output:
        ur1= path.join(PATH_TRIMMED, '{sample}' + PREFIX[0] + '.trimmed.fastq'),
        ur2= path.join(PATH_TRIMMED, '{sample}' + PREFIX[1] + '.trimmed.fastq')
    shell:
        """
        gunzip -f {input.r1} {input.r2}
        """
rule hisat2_alignment:
    input:
       fq1 = path.join(PATH_TRIMMED, '{sample}' + PREFIX[0] + '.trimmed.fastq'),
       fq2 = path.join(PATH_TRIMMED, '{sample}' + PREFIX[1] + '.trimmed.fastq')
    output:
        bam = PATH_BAM + '{sample}.sorted.out.bam',
        sum = PATH_LOG + '{sample}_sum.txt',
        met = PATH_LOG + '{sample}_met.txt'
    log:
        PATH_LOG + 'align/{sample}_hisat2.log'
    params:
        extra = '',
        index = PATH_HISAT2INDEX
    threads:
        config['hisat2']['threads']
    shell:
        """
        hisat2 \
        --threads {threads} \
        -x {params.index} \
        --summary-file {output.sum} \
        --met-file {output.met} \
        -1 {input.fq1} \
        -2 {input.fq2} |\
        samtools sort -T {threads} -o {output.bam}
        """
#    wrapper:#
#        "v1.3.2/bio/hisat2/align" # Map PE reads with HISAT2

rule featureCounts:
    input:
        expand(PATH_BAM + '{sample}.sorted.out.bam', sample=RNAIDs)
    output:
        PATH_OUT + 'full.counts.txt'
    params:
        others = '-p -t exon -g gene_name',
        gtf = PATH_HTSEQ_GTF
    log:
        PATH_LOG + 'featurecount.log'
    threads:
        config['ftc']['threads']
    shell:
        """
        featureCounts {params.others} -T {threads} -a {params.gtf} -o {output} {input} > {log}
        """



rule fastqc_fastq:
    input:
        PATH_BAM + '{sample}.sorted.out.bam'
    output:
        html=PATH_QC+"{sample}_fastqc.html",
        zip=PATH_QC+"{sample}_fastqc.zip"
    log:
        path.join(PATH_QC, '{sample}.fastqc')
    wrapper:
        "0.65.0/bio/fastqc" # generates fastq qc statistics; quality control tool for high throughput sequence data

#rule multiqc:
#    input:
#        expand(PATH_QC+"{sample}_fastqc.zip", sample=RNAIDs)
#    output:
#        PATH_OUT+'multiqc.html'
#    log:
#        PATH_LOG+'multiqc.log'
#    wrapper:
#        '0.65.0/bio/multiqc' # modular tool to aggregate results from bioinformatics analyses across many samples into a single report

rule fastqc_raw:
    input:
        lambda wildcards: [file for file in Files if wildcards.sample in file]
    output:
        html = path.join(PATH_QC, 'raw', '{sample}.html'),
        zip = path.join(PATH_QC, 'raw', '{sample}_fastqc.zip')
    log:
        path.join(PATH_QC, '{sample}.fastqc')
    wrapper:
        'v1.3.2/bio/fastqc'

#rule multiqc_raw:
#    input:
#        lambda wildcards: [path.join(PATH_QC, 'raw', \
#            file.split('/')[-1].split('.fastq.gz')[0] + '_fastqc.zip') \
#            for file in Files]
#    output:
#        PATH_OUT+'multiqc_raw.html'
#    log:
#        PATH_LOG+'multiqc_raw.log'
#    wrapper:
#        '0.65.0/bio/multiqc'

ruleorder: fastqc_raw > hisat2_alignment
