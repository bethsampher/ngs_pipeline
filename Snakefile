rule all:
  input:
    "/home/jovyan/pipeline_data/SRR622461_1_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_2_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_sorted.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_sorted_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_refined.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_refined_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_CYP2C19.vcf"

rule fastp_filter:
  input:
    fwd="{filename}_1.fastq.gz",
    rev="{filename}_2.fastq.gz"
  output:
    fwd=temp("{filename}_1_fastp.fastq.gz"),
    rev=temp("{filename}_2_fastp.fastq.gz"),
    json="{filename}_fastp.json",
    html="{filename}_fastp.html",
    out="{filename}_fastp.txt"
  threads: 16
  shell:
    "fastp --thread 16 -i {input.fwd} -o {output.fwd} -I {input.rev} -O {output.rev} --disable_adapter_trimming --length_required 36 -3 --correction --json {output.json} --html {output.html} 2> {output.out}"

rule fastqc:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  threads: 16
  shell:
    "fastqc -t 16 {input}"

rule bwa_mem:
  input:
    ref="/home/jovyan/pipeline_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    fwd="{filename}_1_fastp.fastq.gz",
    rev="{filename}_2_fastp.fastq.gz"
  output:
    temp("{filename}.sam")
  threads: 16
  shell:
    "bwa mem -t 16 -R '@RG\\tID:1\\tLB:library\\tPL:Illumina\\tPU:lane1\\tSM:NA12878' {input.ref} {input.fwd} {input.rev} > {output}"

rule sam_to_bam:
  input:
    "{filename}.sam"
  output:
    temp("{filename}.bam")
  threads: 16
  shell:
    "samtools view -b {input} -o {output} -@ 16"

rule samtools_sort:
  input:
    "{filename}.bam"
  output:
    "{filename}_sorted.bam"
  threads: 16
  shell:
    "samtools sort {input} -o {output} -@ 16"

rule samtools_index:
  input:
    "{filename}.bam"
  output:
    "{filename}.bam.bai"
  threads: 16
  shell:
    "samtools index {input} -@ 16"

rule samtools_flagstat:
  input:
    "{filename}.bam"
  output:
    "{filename}_flagstats.txt"
  threads: 16
  shell:
    "samtools flagstat {input} -@ 16 > {output}"

rule picard_remove_duplicates:
  input:
    "{filename}_sorted.bam"
  output:
    bam="{filename}_refined.bam",
    metrics="{filename}_dupl_metrics.txt"
  shell:
    "picard MarkDuplicates I={input} O={output.bam} METRICS_FILE={output.metrics} REMOVE_DUPLICATES=true"

rule gatk_haplotype_caller:
  input:
    ref="/home/jovyan/pipeline_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    bam="{filename}_refined.bam"
  output:
    vcf="{filename}.vcf",
    idx="{filename}.vcf.idx"
  shell:
    "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} -mbq 20"

rule vcftools_filter:
  input:
    "{filename}.vcf"
  output:
    temp("{filename}.recode.vcf"),
    "{filename}.log"
  shell:
    "vcftools --vcf {input} --out {wildcards.filename} --minDP 3 --minQ 20 --recode --recode-INFO-all"

rule vcftools_exclude:
  input:
    "{filename}.recode.vcf"
  output:
    temp("{filename}_filtered.recode.vcf"),
    "{filename}_filtered.log"
  shell:
    "vcftools --vcf {input} --out {wildcards.filename} --max-missing 1 --recode --recode-INFO-all"

rule snpeff_annotate:
  input:
    "{filename}_filtered.recode.vcf"
  output:
    vcf="{filename}_annotated.vcf",
    html="{filename}_snpEff_summary.html",
    txt="{filename}_snpEff_summary.genes.txt"
  shell:
    "snpEff -Xmx4G GRCh38.104 {input} -stats {output.html} > {output.vcf}"

rule snpsift_filter:
  input:
    "{filename}_annotated.vcf"
  output:
    "{filename}_CYP2C19.vcf"
  shell:
    "cat {input} | SnpSift filter \"ANN[*].GENE = 'CYP2C19'\" > {output}"

