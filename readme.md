# Prepare

Set up the `configs/config_main.yaml`, `configs/metadata.csv`



# run

```
nohup snakemake -p --cores 20 &
snakemake -np
```



# Other Script

## DIY hisat2 index

``` 
hisat2-build Homo_sapiens.GRCh38_release95.genome.fa Homo_sapiens.GRCh38_release95.genome
```

fa need to be decompressed 



## rename

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

