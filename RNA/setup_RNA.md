# Conda

```
conda install -y mamba
mamba create -n rnaseq -y python pandas snakemake pigz fastqc trim-galore multiqc samtools salmon hisat2 subread 
conda activate rnaseq
```



# Ref

human ref47 files from gencode:

```
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.primary_assembly.annotation.gtf.gz &
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.transcripts.fa.gz &
```



# salmon

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



# Hisat2

```
nohup hisat2-build -p 32 ../GRCh38.primary_assembly.genome.fa genecode_v47.genome &
```





# run

```
nohup snakemake -p --cores 20 -s Snakefile_hisat &
# snakemake -np

snakemake -np -s Snakefile_rna_v02 

nohup snakemake -p --cores 8 -s Snakefile &
```



# Other Script

## DIY hisat2 index

## rename

```
mkdir 0-raw
```



```
ls | while read id
do
cd $id
mv *.fq.gz ~/2-Project/dongyang/0520-rnaseq/0-raw/
cd ..
done
```



## make metadata

```
ls *.fq.gz | while read id
do
base=${id%_*}
echo $base
done | uniq > ../configs/metadata.csv

sed -i '1s/^/sample\n/' configs/metadata.csv
```



# ref

[Fred-White94/snakemake_rnaseq: A Snakemake pipeline to go from fastq mRNA sequencing files to raw and normalised counts (usable for downstream EDA and differential analysis)](https://github.com/Fred-White94/snakemake_rnaseq)

[tjbencomo/ngs-pipeline: Pipeline for Somatic Variant Calling with WES and WGS data](https://github.com/tjbencomo/ngs-pipeline)

[zhxiaokang/RASflow: RNA-Seq analysis workflow](https://github.com/zhxiaokang/RASflow)

[基于GATK4标准找变异方法的自动化工作流程oVarFlow的使用-腾讯云开发者社区-腾讯云](https://cloud.tencent.com/developer/article/2032022)
