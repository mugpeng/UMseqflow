rule qc:
    input:
        fq=expand("raw/{sample}_{strand}.fq.gz", sample=SAMPLES, strand=STRANDS)
    output:
        qc_zip=expand("fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS)
    threads: THREADS
    log:
        "logs/fastqc/fastqc.log"
    shell:
        "fastqc -t {threads} {input.fq} -o ./fastqc/ &> {log}"

rule trim:
    input:
        fq1 = "raw/{sample}_" + STRANDS[0] + ".fq.gz",
        fq2 = "raw/{sample}_" + STRANDS[1] + ".fq.gz"
    output:
        fq1 = "clean/{sample}_1_val_" + STRANDS[0] + ".fq.gz",
        fq2 = "clean/{sample}_2_val_" + STRANDS[1] + ".fq.gz",
        qc1 = "clean/{sample}_1_val_" + STRANDS[0] + "_fastqc.zip",
        qc2 = "clean/{sample}_2_val_" + STRANDS[1] + "_fastqc.zip"
    threads: THREADS
    params:
        quality = 20,
        length = 20
    log:
        "logs/trim/{sample}.log"
    shell:
        "trim_galore -j {threads} --phred33 -q {params.quality} "
        "--length {params.length} --fastqc --paired "
        "-o ./clean/ {input.fq1} {input.fq2} &> {log}" 