# HISAT2 mode rules
rule map_sort:
    input:
        fq1 = "clean/{sample}_1_val_" + STRANDS[0] + ".fq.gz",
        fq2 = "clean/{sample}_2_val_" + STRANDS[1] + ".fq.gz"
    output:
        bam = "bam/{sample}_hisat_sorted.bam",
        bai = "bam/{sample}_hisat_sorted.bam.bai"
    threads: THREADS
    params:
        index = HISAT2_INDEX
    log:
        hisat = "logs/hisat2/{sample}.log",
        samtools = "logs/samtools/{sample}.log"
    shell:
        """
        hisat2 -p {threads} -x {params.index} \
            -1 {input.fq1} -2 {input.fq2} 2> {log.hisat} | \
        samtools sort -@ {threads} -o {output.bam} - &> {log.samtools} && \
        samtools index {output.bam} {output.bai}
        """

rule featureCounts:
    input:
        bam = expand("bam/{sample}_hisat_sorted.bam",sample=SAMPLES)
    output:
        counts = "counts/all.id.txt",
        summary = "counts/all.id.txt.summary"
    threads: THREADS
    params:
        gtf = GTF
    log:
        "logs/featureCounts/counts.log"
    shell:
        "featureCounts -T {threads} -p -t exon -g gene_name "
        "-a {params.gtf} -o {output.counts} {input.bam} &> {log}" 