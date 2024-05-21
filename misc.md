# rename

```
ls | while read id
do
cd $id
mv *.fq.gz ~/2-Project/dongyang/0520-rnaseq/0-raw/
cd ..
done
```



# make metadata

```
ls *.fq.gz | while read id
do
base=${id%_*}
echo $base
done | uniq > ../configs/metadata.csv

sed -i '1s/^/sample\n/' configs/metadata.csv
```



# run

```
nohup snakemake -p --cores 48 &
snakemake -np
```

