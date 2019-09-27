# Trichuris_transwells
Scripts for analysing Trichuris transwells data. 

Find files:

```
tail -n +2  meta_data.txt | cut -f 1 | while read -r sample; do 
	pf data -i ${sample} -t sample --filetype fastq) --symlink ./fastq --rename
	MY_PATHS=($(pf data -i ${sample} -t sample --filetype fastq))
	for MY_PATH in ${MY_PATHS[@]}; do
		echo ${sample}$'\t'${MY_PATH} >> data_locations.txt
	done
done
```

Generate STAR genome indices for mapping:

```
bsub -o index.o -e index.e -R'select[mem>=20000] rusage[mem=20000] span[hosts=1]' -M 20000 -n 8 STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm38.98.gtf
```

Map:


