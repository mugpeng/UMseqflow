rule multiqc_all:
    input:
        raw_qc = expand("1-fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS) if qc == "yes" else [],
        clean_qc = expand("2-clean/{sample}_{strand}_val_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS) if trim == "yes" else [],
        hisat_counts = "4-counts/all.id.txt.summary" if HISAT2_INDEX else [],
        hisat_logs = expand("logs/hisat2/{sample}.log", sample=SAMPLES) if HISAT2_INDEX else [],
        salmon_outputs = expand("salmon/quant/{sample}.count/", sample=SAMPLES) if SALMON_INDEX else []
    output:
        raw_report = "multiqc/raw_multiqc_report.html" if qc == "yes" else touch("multiqc/.raw_dummy"),
        clean_report = "multiqc/clean_multiqc_report.html" if trim == "yes" else touch("multiqc/.clean_dummy"),
        count_report = "multiqc/count_multiqc_report.html" if HISAT2_INDEX else touch("multiqc/.count_dummy"),
        bam_report = "multiqc/bam_multiqc_report.html" if HISAT2_INDEX else touch("multiqc/.bam_dummy"),
        salmon_report = "multiqc/salmon_multiqc_report.html" if SALMON_INDEX else touch("multiqc/.salmon_dummy")
    shell:
        """
        if [ "{qc}" = "yes" ]; then
            multiqc {input.raw_qc} --filename {output.raw_report}
        fi
        
        if [ "{trim}" = "yes" ]; then
            multiqc {input.clean_qc} --filename {output.clean_report}
        fi
        
        if [ "{HISAT2_INDEX}" != "" ]; then
            multiqc {input.hisat_counts} --filename {output.count_report}
            multiqc {input.hisat_logs} --filename {output.bam_report}
        fi
        
        if [ "{SALMON_INDEX}" != "" ]; then
            multiqc {input.salmon_outputs} --filename {output.salmon_report}
        fi
        """ 