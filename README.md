# bioinformatics 1

1. Скачал E. Coli SRR31294323: https://www.ncbi.nlm.nih.gov/sra/SRX26673257[accn]
! Отмечу, из пайплайна убрал сам процесс скачивания, так как столкнулся с неведомыми проблемами

```shell
wget -O SRR31294323.fastq.gz trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR31294323 && gunzip SRR31294323.fastq.gz
```

2. Скачивание референса
```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz && mv GCF_000005845.2_ASM584v2_genomic.fna.gz reference.fna.gz && gunzip reference.fna.gz
```

3. Установка conda, minimap, fastqc, samtools
fastqc:
```shell
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
```

samtools:
```shell
git clone https://github.com/samtools/samtools.git
cd samtools && make
```

4. Анализ fastqc
```shell
fastqc SRR31294323.fastq
```

5. Индексация референсного генома
```shell
minimap2 -d reference.mmi reference.fna
```

6. Картирование ридов на референсный геном
```shell
minimap2 -a reference.mmi SRR31294323.fastq > alignment.sam
```

7. Преобразование SAM в BAM
```shell
samtools view -Sb alignment.sam > alignment.bam
```

8. Сортировка BAM
```shell
samtools sort alignment.bam -o alignment.sorted.bam
```

9. Оценка качества картирования
```shell
samtools flagstat alignment.sorted.bam > flagstat_output.txt
```

10. Генерация визуализации пайплайна (DAG)
```shell
snakemake --dag | dot -Tpng > {pipeline.png}
```

11. Запуск пайплайна
```shell
snakemake --cores 2
```
