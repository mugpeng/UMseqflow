import pandas as pd
import yaml

with open('configs/config.yaml') as yamlfile:
    config = yaml.safe_load(yamlfile)

# Load parameters to the workflow
THREADS = config["THREADS"]
GTF = config["GTF"]
HISAT2_INDEX = config.get('HISAT2_INDEX', "")
SALMON_INDEX = config.get('SALMON_INDEX', "")

qc = config["QC"]
trim = config["TRIMMED"]

# read samples
SAMPLES = pd.read_csv(config["METAFILE"], sep = '\t', header = 0)['sample']
STRAND_PATTERN = ""  # Change this to "R" for "R1"/"R2" pattern or "" for "1"/"2" pattern
STRANDS = [f"{STRAND_PATTERN}1", f"{STRAND_PATTERN}2"]  # Will generate either ["R1", "R2"] or ["1", "2"]

def get_final_output():
    final_outputs = []
    
    # Add HISAT2 outputs if HISAT2_INDEX is provided
    if HISAT2_INDEX:
        final_outputs.extend([
            "multiqc/count_multiqc_report.html",
            "multiqc/bam_multiqc_report.html",
            "counts/all.id.txt"
        ])
    
    # Add Salmon outputs if SALMON_INDEX is provided
    if SALMON_INDEX:
        final_outputs.extend([
            "salmon/salmon.genes",
            "multiqc/salmon_multiqc_report.html"
        ])
    
    # Add QC outputs regardless of mode
    if qc == "yes":
        final_outputs.extend(expand("fastqc/{sample}_{strand}_fastqc.zip", sample=SAMPLES, strand=STRANDS))
    
    return final_outputs 