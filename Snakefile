# The main script to manage the subworkflows of UMseqflow
import pandas as pd
import yaml

with open('configs/config_main.yaml') as yamlfile:
    config = yaml.safe_load(yamlfile)

# Load parameters to the workflow
THREADS = config["THREADS"]
GTF = config["GTF"]
HISAT2_INDEX = config['HISAT2_INDEX']

# read samples
SAMPLES = pd.read_csv(config["METAFILE"], sep = '\t', header = 0)['sample']
STRANDS = ["1", "2"]

rule all: 
    input:
        "multiqc/count_multiqc_report.html",
        "4-counts/all.id.txt"
        

rule qc:
	input:
		fq=expand("0-raw/{sample}_{strand}.fq.gz", sample=SAMPLES, strand=STRANDS)
	output:
		qc_zip=expand("1-fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS)
	threads: THREADS
	shell:
		"fastqc -t {threads} {input.fq} -o ./1-fastqc/"

rule trim:
    input:
        fq1 = "0-raw/{sample}_1.fq.gz",
        fq2 = "0-raw/{sample}_2.fq.gz",
    output:
        "2-clean/{sample}_1_val_1.fq.gz",
        "2-clean/{sample}_2_val_2.fq.gz",
        "2-clean/{sample}_1_val_1_fastqc.zip",
        "2-clean/{sample}_2_val_2_fastqc.zip"
    threads: THREADS
    shell:
        "trim_galore -j {threads} --phred33 -q 20 --length 20 --fastqc --paired -o ./2-clean/ {input.fq1} {input.fq2}"

rule map_sort:
    input:
        fq1 = "2-clean/{sample}_1_val_1.fq.gz",
        fq2 = "2-clean/{sample}_2_val_2.fq.gz"
    output:
        bam = "3-bam/{sample}_hisat_sorted.bam",
        bai = "3-bam/{sample}_hisat_sorted.bam.bai",
        depth = "3-bam/{sample}_hisat_sorted.bam.depth.txt",
        log = "3-bam/{sample}.log"
    params:
        index = HISAT2_INDEX
    threads: THREADS
    shell:
        "hisat2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} 2>{output.log} | samtools sort -@ {threads} -o {output.bam} -  && samtools index {output.bam} {output.bai} && samtools depth {output.bam} > {output.depth}"

rule featureCounts:
    input:
        bam = expand("3-bam/{sample}_hisat_sorted.bam",sample=SAMPLES)
    output:
        counts = "4-counts/all.id.txt",
        summary = "4-counts/all.id.txt.summary"
    params:
        gtf = GTF
    threads: THREADS
    shell:
    	"featureCounts -T {threads} -p -t exon -g gene_id -a {params.gtf} -o {output.counts} {input.bam}"

rule muitiqc_all:
	input:
       		raw_qc = expand("1-fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS),
		clean_qc = expand("2-clean/{sample}_{strand}_val_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS),
        	bam_qc = expand("3-bam/{sample}.log", sample=SAMPLES),
        	count_qc = "4-counts/all.id.txt.summary"
	output:
        	raw_report = "multiqc/raw_multiqc_report.html",
		clean_report = "multiqc/clean_multiqc_report.html",
        	bam_report = "multiqc/bam_multiqc_report.html",
        	count_report = "multiqc/count_multiqc_report.html"
	shell:
		'''
        	multiqc {input.raw_qc} --filename {output.raw_report}
        	multiqc {input.clean_qc} --filename {output.clean_report}
        	multiqc {input.bam_qc} --filename {output.bam_report}
        	multiqc {input.count_qc} --filename {output.count_report}
        	'''
