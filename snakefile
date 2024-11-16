# Главная цель — получить флагстат (анализ качества картирования)
rule all:
    input:
        "flagstat_output.txt"

# Шаг 1. Скачивание референсного генома
rule download_reference:
    output:
        "reference.fna"
    shell:
        "wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz && mv GCF_000005845.2_ASM584v2_genomic.fna.gz reference.fna.gz && gunzip reference.fna.gz"

# Шаг 2. Скачивание ридов
rule download_reads:
    input:
        "SRR31294323.fastq"
    output:
        "SRR31294323.html"
    shell:
        "fastqc {input}"

# Шаг 3. Индексация референсного генома
rule index_reference:
    input:
        "reference.fna"
    output:
        "reference.mmi"
    shell:
        "minimap2 -d {output} {input}"

# Шаг 4. Картирование ридов на референсный геном
rule map_reads:
    input:
        ref="reference.mmi",
        reads="SRR31294323.fastq"
    output:
        "alignment.sam"
    shell:
        "minimap2 -a {input.ref} {input.reads} > {output}"

# Шаг 5. Преобразование SAM в BAM
rule sam_to_bam:
    input:
        "alignment.sam"
    output:
        "alignment.bam"
    shell:
        "samtools view -Sb {input} > {output}"

# Шаг 6. Сортировка BAM
rule sort_bam:
    input:
        "alignment.bam"
    output:
        "alignment.sorted.bam"
    shell:
        "samtools sort {input} -o {output}"

# Шаг 7. Оценка качества картирования
rule flagstat:
    input:
        "alignment.sorted.bam"
    output:
        "flagstat_output.txt"
    shell:
        "samtools flagstat {input} > {output}"

# Шаг 8. Генерация визуализации пайплайна (DAG)
rule visualize_dag:
    output:
        "pipeline.png"
    shell:
        "snakemake --dag | dot -Tpng > {output}"

