# Salmon mode rules
rule salmon_quant:
    input:
        fq1 = "clean/{sample}_1_val_" + STRANDS[0] + ".fq.gz",
        fq2 = "clean/{sample}_2_val_" + STRANDS[1] + ".fq.gz"
    output:
        quant = directory("salmon/quant/{sample}.count")
    threads: THREADS
    params:
        index = SALMON_INDEX,
        gtf = GTF
    log:
        "logs/salmon/{sample}.log"
    shell:
        "salmon quant --gcBias -l A "
        "-1 {input.fq1} -2 {input.fq2} "
        "-i {params.index} "
        "-g {params.gtf} "
        "-o {output.quant} "
        "-p {threads} &> {log}"

rule salmon_merge:
    input:
        quants = expand("salmon/quant/{sample}.count", sample=SAMPLES)
    output:
        merged = "salmon/salmon.genes"
    shell:
        "salmon quantmerge --column NumReads --genes "
        "--quants {input.quants} "
        "-o {output.merged}" 