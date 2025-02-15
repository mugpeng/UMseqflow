# UMseqflow

Welcome. UMseqflow is an easy and updated snakemake bioinformatics workflow for simplying the analysis workflow.



structure:

```
.
├── LICENSE
├── RNA
│   ├── Snakefile_rna_v04
│   ├── configs
│   │   ├── config.yaml
│   │   └── metadata.csv
│   └── rules
│       ├── common.smk
│       ├── hisat.smk
│       ├── multiqc.smk
│       ├── preprocessing.smk
│       └── salmon.smk
├── configs
│   ├── config.yaml
│   └── metadata.csv
├── readme.md
└── setup.md
```



Only provide RNA mode currently.

![](http://cos01.mugpeng.top/img/workflow.png)



# TODO

- Major

- [ ] WES mode
- [ ] package into docker
- [x]  Publish to GitHub



- Minor

- [ ] Revise trim mode
- [ ] Update Ref files into ali netdisk
- [ ] trim part qc is slow(run trim, run other part)



# Prepare

First prepare analysis environment follow `setup.md` include conda, ref files.



Set up the `configs/config.yaml`, `configs/metadata.csv`



First move or use script from `Other useful sript` by soft link to `./raw` folder:

```
$ ls raw
HB_1_1.fq.gz  HB_1_2.fq.gz  UH_1_1.fq.gz  UH_1_2.fq.gz
```



Then run for making  `metadata.csv` :

```
echo "sample" > metadata.csv
ls raw/*.fq.gz | sed 's/_[12].fq.gz$//' | sed 's|^raw/||' | sort | uniq >> metadata.csv
mv metadata.csv configs/metadata.csv
```



## RNA

You can decide call salmon, hisat2+featurecount or both by leaving the path for these index empty or not. 

![](http://cos01.mugpeng.top/img/workflow.png)



# Run

```
nohup snakemake -p --cores 42 -s Snakefile_rna_v04 &

snakemake -np -s Snakefile_rna_v04
# dry-mode for test
```



You can visualize the pipeline through `graphviz`:

```
snakemake --dag -s RNA/Snakefile_rna_v04 | dot -Tpdf > workflow.pdf

snakemake --dag -s Snakefile_rna_v04 | dot -Tpng > workflow.png
```





# Other useful sript

- soft link all fq into raw folder

```
find ../../X201SC24128617-Z01-F001/01.RawData -type f -name '*fq.gz' -exec ln -s {} . \;
```



# Milestones

## 250215

publish to github. Welcome!



# Contribution

You can follow other projects which I also referenced:

[Fred-White94/snakemake_rnaseq: A Snakemake pipeline to go from fastq mRNA sequencing files to raw and normalised counts (usable for downstream EDA and differential analysis)](https://github.com/Fred-White94/snakemake_rnaseq)

[tjbencomo/ngs-pipeline: Pipeline for Somatic Variant Calling with WES and WGS data](https://github.com/tjbencomo/ngs-pipeline)

[zhxiaokang/RASflow: RNA-Seq analysis workflow](https://github.com/zhxiaokang/RASflow)

[基于GATK4标准找变异方法的自动化工作流程oVarFlow的使用-腾讯云开发者社区-腾讯云](https://cloud.tencent.com/developer/article/2032022)



toturial:
[Snakemake for Biostatistics Quick-Start Tutorial](https://jean997.github.io/snakemake_tutorial/index.html)

https://carpentries-incubator.github.io/snakemake-novice-bioinformatics/

[#snakemake 北野茶缸子](https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzU5ODc3OTA0NQ==&action=getalbum&album_id=2426609164623003652&scene=173&sessionid=undefined&enterid=0&from_msgid=2247490534&from_itemidx=4&count=3&nolastread=1&uin=&key=&devicetype=iMac+Mac16%2C10+OSX+OSX+15.0+build(24A8332)&version=13080911&lang=en&nettype=WIFI&ascene=78&fontScale=100)



Or you can connect with me: mugpeng@outlook.com



Cursor, copilot, trae are good LLM IDE for coding your smk!



# Why I am choosing Smk instead of others

nextflow, galaxy...
