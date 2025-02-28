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
grep "^>" GRCh38.primary_assembly.genome.fa.gz | cut -d " " -f 1 > salmon/decoys.txt
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



mouse:

```
mkdir salmon
grep "^>" GRCm39.primary_assembly.genome.fa | cut -d " " -f 1 > salmon/decoys.txt
cd salmon
sed -i.bak -e 's/>//g' decoys.txt

zcat ../gencode.vM36.transcripts.fa.gz <(cat ../GRCh38.primary_assembly.genome.fa) | pigz -p 32 > gentrome.fa.gz
nohup salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode -k 31 &
```





## Hisat2

```
nohup hisat2-build -p 32 ../GRCh38.primary_assembly.genome.fa genecode_v47.genome &

nohup hisat2-build -p 32 ../GRCm39.primary_assembly.genome.fa genecode_v36.genome &
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

zcat gencode.vM36.primary_assembly.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sort -k1,1 -k2,2n| bedtools merge > GRCh39_gencode_gtf_exon.bed
```



## GATK build GRCm39 resource bundle

Ref:

[liftover,crossmap进行坐标转换时用到的chain文件介绍 - Zhongxu blog](https://www.zxzyl.com/archives/838/)

[为小鼠(GRCm38) 构建 GATK resource bundle - 简书](https://www.jianshu.com/p/dfe67e99151c)
[genomics/workflows/gatk-mouse-mm10.md at master · igordot/genomics](https://github.com/igordot/genomics/blob/master/workflows/gatk-mouse-mm10.md)



about vcf file: [30. 理解vcf 文件](https://www.yuque.com/mugpeng/sequence/bi0q4a)



download GRCm38 indel and snp:

```
# SNP
nohup wget -c https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz &

nohup wget -c https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz.tbi &

# INDEL
nohup wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz & 
nohup wget -c ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz.tbi &
```



SNP:

``` 
# from GRCm38 to mm10
zcat GRCm38/00-All.vcf.gz | sed 's/^\([0-9XY]\)/chr\1/' > mm10/00-All.vcf

# from mm10 to GRCm39
# nohup CrossMap vcf crossmap/mm10ToMm39.over.chain.gz mm10/00-All.vcf ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa GRCm39/00-All.vcf &

# crossmap it after sort
```



Indel:

```
# from GRCm38 to mm10
# adjust header
zcat GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> mm10/mgp.v5.indels.dbSNP142.vcf
# keep only passing and adjust chromosome name
zcat GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | cut -f 1-8 \
| sed 's/^\([0-9MXY]\)/chr\1/' \
>> mm10/mgp.v5.indels.dbSNP142.vcf
# filter pass
zcat GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | head -1000 | grep "^#" | cut -f 1-8 \
| grep -v "#contig" | grep -v "#source" \
> mm10/mgp.v5.indels.dbSNP142.pass.vcf
# keep only passing and adjust chromosome name
zcat GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz | grep -v "^#" | cut -f 1-8 \
| grep -w "PASS" | sed 's/^\([0-9MXY]\)/chr\1/' \
>> mm10/mgp.v5.indels.dbSNP142.pass.vcf &

# from mm10 to GRCm39
nohup CrossMap vcf crossmap/mm10ToMm39.over.chain.gz mm10/mgp.v5.indels.dbSNP142.pass.vcf ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa GRCm39/mgp.v5.indels.dbSNP142.pass.vcf &
```



sort vcf or it will produce errors:

```
> [E::hts_idx_push] Unsorted positions on sequence #3: 7719079 followed by 7719047
tbx_index_build3 failed: GRCm39/mgp.v5.indels.dbSNP142.pass.vcf.gz
```



```
samtools dict ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa > ../genome/gencode_v36/GRCm39.primary_assembly.genome.dict
# sort vcf
nohup gatk SortVcf  \
  -I GRCm39/mgp.v5.indels.dbSNP142.pass.vcf \
  -O GRCm39/mgp.v5.indels.dbSNP142.pass.sort.vcf \
  # -R ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa \
  -SD ../genome/gencode_v36/GRCm39.primary_assembly.genome.dict &


# FOR SNP 00-All.vcf file
nohup gatk UpdateVcfSequenceDictionary \
 I=mm10/00-All.vcf \
 O=mm10/00-All.update.vcf \
 SEQUENCE_DICTIONARY=../genome/gencode_v36/GRCm39.primary_assembly.genome.dict &

nohup gatk --java-options "-Xmx100G -Djava.io.tmpdir=./" SortVcf \
   -I mm10/00-All.update.vcf \
   -O mm10/00-All.sort.vcf \
   -SD ../genome/gencode_v36/GRCm39.primary_assembly.genome.dict &

# nohup gatk SortVcf  \
#    -I GRCm39/00-All.update.vcf \
#    -O GRCm39/00-All.sort.vcf \
#    -SD ../genome/gencode_v36/GRCm39.primary_assembly.genome.dict &

nohup CrossMap vcf crossmap/mm10ToMm39.over.chain.gz mm10/00-All.sort.vcf ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa GRCm39/00-All.sort.vcf &

# nohup bcftools reheader -f ../genome/gencode_v36/GRCm39.primary_assembly.genome.fa GRCm39/00-All.vcf > GRCm39/00-All.sort.vcf &
```

Use R parameter will produce: `Caused by: java.lang.AssertionError: SAM dictionaries are not the same`



compress vcf and make index:

```
# bgzip -@ 32 -d GRCm39/00-All.vcf.gz &
# bgzip -@ 32 -d GRCm39/mgp.v5.indels.dbSNP142.pass.vcf.gz &
bgzip -@ 32 GRCm39/00-All.vcf &
bgzip -@ 32 GRCm39/mgp.v5.indels.dbSNP142.pass.vcf &

tabix -p vcf GRCm39/00-All.vcf.gz &
tabix -p vcf GRCm39/mgp.v5.indels.dbSNP142.pass.vcf.gz &
```





## ANNOVAR

For annotation.



## Copy number

[Bioinformatics Pipeline: DNA-Seq Analysis: Whole Genome Sequencing Variant Calling - GDC Docs](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_WGS/)



### Gistic2

Ref: [GISTIC Documentation](https://broadinstitute.github.io/gistic2/)

fetch mouse gistic2 mat: [roland-rad-lab/MoCaSeq: Analysis pipelines for cancer genome sequencing in mice.](https://github.com/roland-rad-lab/MoCaSeq)

