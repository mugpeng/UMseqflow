# Conda

```
conda install -y mamba
mamba create -n rnaseq -y python pandas snakemake pigz fastqc trim-galore multiqc samtools salmon hisat2 subread 
conda activate rnaseq
```



# Ref

I use ref files from gencode.



human ref47 files from gencode:

```
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz &
```



Mouse ref36:

```
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/GRCm39.primary_assembly.genome.fa.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.primary_assembly.annotation.gtf.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz &
```





# RNA

## salmon

Ref: [Selective Alignment](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)



make decoy seq:

```
mkdir salmon
grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > salmon/decoys.txt
cd salmon
sed -i.bak -e 's/>//g' decoys.txt
```



Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
```
cat ../gencode.v47.transcripts.fa.gz ../GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
```

```
# ./bin/salmon index -t transcripts.fa -i transcripts_index --decoys decoys.txt -k 31
nohup salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode -k 31 &
```



## Hisat2

```
nohup hisat2-build -p 32 ../GRCh38.primary_assembly.genome.fa genecode_v47.genome &
```



# WES

## bwa

Bwa bulid:

human:

```
pigz -d -p 16 GRCh38.primary_assembly.genome.fa.gz
nohup bwa index -a bwtsw -p bwa /home/data/yzpeng/3-Data/Ref/human/genome/gencode_v47/GRCh38.primary_assembly.genome.fa &
```

mouse:

```
pigz -d -p 16 GRCm39.primary_assembly.genome.fa.gz &
nohup bwa index -a bwtsw -p bwa/GRCm39 /home/data/yzpeng/3-Data/Ref/mouse/genome/gencode_v36/GRCm39.primary_assembly.genome.fa &
```



## GATK

GATK4 snp and indel:

```
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz &
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi &
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz &
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi &
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz &
nohup wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi &
```

[Resource bundle – GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)

[Google Cloud console](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/)

GATK4 build ref fa index:

```
samtools faidx GRCh38.primary_assembly.genome.fa
samtools dict GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.dict
```



bed files for mutect2:

```
zcat gencode.v47.primary_assembly.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sort -k1,1 -k2,2n| bedtools merge > hg38_gencode_gtf_exon.bed
```



## ANNOVAR

For annotation.



