SAMPLES = []
STRANDS = ["1", "2"]
# REF = "/home/data/yzpeng/3-Data/Ref/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF = "/home/data/yzpeng/3-Data/Ref/human/Homo_sapiens.GRCh38.111.gtf"
# HISAT2_REF = "/home/data/yzpeng/3-Data/Ref/human/hisat_test/test_10M"
THREADS = 16

rule all: 
    input:
        "results/counts/all.id.txt",
        "results/fastqc/multiqc_report.html",
        "data/clean/multiqc_report.html"
        

rule qc:
	input:
		fq=expand("data/samples/{sample}_R{strand}.fastq.gz", sample=SAMPLES, strand=STRANDS)
	output:
            "results/fastqc/multiqc_report.html",
		qc_zip=expand("results/fastqc/{sample}_R{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS)
	threads: THREADS
	shell:
		"fastqc -t {threads} {input.fq} -o ./results/fastqc/ && multiqc {output.qc_zip} -o ./results/fastqc/"

rule trim:
    input:
        fq1 = "data/samples/{sample}_R1.fastq.gz",
        fq2 = "data/samples/{sample}_R2.fastq.gz",
    output:
        "data/clean/{sample}_R1_val_1.fq.gz",
        "data/clean/{sample}_R2_val_2.fq.gz",
        "data/clean/{sample}_R1_val_1_fastqc.zip",
        "data/clean/{sample}_R2_val_2_fastqc.zip"
    threads: THREADS
    shell:
        "trim_galore -j {threads} --phred33 -q 20 --length 20 --fastqc --paired -o ./data/clean/ {input.fq1} {input.fq2}"

rule qc_after_trim:
	input:
		expand("data/clean/{sample}_R{strand}_val_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS)
	output:
		"data/clean/multiqc_report.html"
	shell:
		"multiqc {input} -o ./data/clean/"

rule map_sort:
    input:
        fq1 = "data/clean/{sample}_R1_val_1.fq.gz",
        fq2 = "data/clean/{sample}_R2_val_2.fq.gz"
#        ref = HISAT2_REF

    output:
        bam = "results/sorted_reads/{sample}_hisat_sorted.bam",
        bai = "results/sorted_reads/{sample}_hisat_sorted.bam.bai",
        depth = "results/sorted_reads/{sample}_hisat_sorted.bam.depth.txt"
    threads: THREADS
    shell:
        "hisat2 -p {threads} -x /home/data/yzpeng/3-Data/Ref/human/hisat2/Homo_sapiens.GRCh38.dna.primary_assembly -1 {input.fq1} -2 {input.fq2} | samtools sort -@ {threads} -o {output.bam} -  && samtools index {output.bam} {output.bai} && samtools depth {output.bam} > {output.depth}"

rule featureCounts:
    input:
        bam = expand("results/sorted_reads/{sample}_hisat_sorted.bam",sample=SAMPLES),
        gtf = GTF
    output:
        counts = "results/counts/all.id.txt",
        summary = "results/counts/all.id.txt.summary"
    threads: THREADS
    shell:
    	"featureCounts -T {threads} -p -t exon -g gene_id -a {input.gtf} -o {output.counts} {input.bam}"