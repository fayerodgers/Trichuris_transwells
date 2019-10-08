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
bsub -o index.o -e index.e -R'select[mem>=25000] rusage[mem=25000] span[hosts=1]' -M 25000 -n 8 STAR --runThreadN 8 --runMode genomeGenerate --genomeDir . --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm38.98.gtf
```

Map:

```
tail -n +2 meta_data.txt| cut -f 1 | while read -r sample; do 
	file=$(grep ${sample} data_locations.txt | grep -o [^/]*.fastq.gz | sed -e 's/#/_/g' | sed -e 's/_[12].fastq.gz//g' | sort -u )
	mkdir ./mapping/${sample}
	bsub -o ./mapping/${sample}/map.o -e ./mapping/${sample}/map.e -R'select[mem>=30000] rusage[mem=30000] span[hosts=1]' -M 30000 -n 8 STAR \
	--runThreadN 8 \
	--genomeDir ./genome_index \
	--readFilesIn ./fastq/${file}_1.fastq.gz ./fastq/${file}_2.fastq.gz \
	--readFilesCommand gunzip -c \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outFileNamePrefix ./mapping/${sample}/${sample} \
	--outSAMtype BAM Unsorted 
done 
```
Count features:

```
tail -n +2 meta_data.txt| cut -f 1 | while read -r sample; do 
	bsub -o ./mapping/${sample}/count.o -e ./mapping/${sample}/count.e -R 'select[mem>=200] rusage[mem=200]' -M 200 'htseq-count -f bam -s no -m union ./mapping/${sample}/${sample}Aligned.out.bam genome_index/Mus_musculus.GRCm38.98.gtf > ./mapping/${sample}/${sample}.count'
done
	
```

