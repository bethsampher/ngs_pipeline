rule all:
  input:
    "/home/jovyan/pipeline_data/SRR622461_1_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_2_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_fastp.txt",
    "/home/jovyan/pipeline_data/SRR622461_1_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_2_fastp_fastqc.html",
    "/home/jovyan/pipeline_data/SRR622461_sorted.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_sorted_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_refined.bam.bai",
    "/home/jovyan/pipeline_data/SRR622461_refined_flagstats.txt",
    "/home/jovyan/pipeline_data/SRR622461_coverage.txt",
    "/home/jovyan/pipeline_data/SRR622461_CYP2C19.vcf"

rule fastqc:
  input:
    "{filename}.fastq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"

rule fastp:
  input:
    fwd="{filename}_1.fastq.gz",
    rev="{filename}_2.fastq.gz"
  output:
    fwd="{filename}_1_fastp.fastq.gz",
    rev="{filename}_2_fastp.fastq.gz",
    json="{filename}_fastp.json",
    html="{filename}_fastp.html",
    cmd="{filename}_fastp.txt"
  shell:
    "fastp --thread 8 -i {input.fwd} -o {output.fwd} -I {input.rev} -O {output.rev} --disable_adapter_trimming --length_required 36 -3 --correction --json {output.json} --html {output.html} > {output.cmd}"

rule bwa_mem:
  input:
    ref="/home/jovyan/pipeline_data/GCF_000001405.26_GRCh38_genomic.fna",
    fwd="{filename}_1_fastp.fastq.gz",
    rev="{filename}_2_fastp.fastq.gz"
  output:
    "{filename}.sam"
  shell:
    "bwa mem -t 8 {input.ref} {input.fwd} {input.rev} > {output}"

rule sam_to_bam:
  input:
    "{filename}.sam"
  output:
    "{filename}.bam"
  shell:
    "samtools view -b {input} -o {output} -@ 8"

rule sort_bam:
  input:
    "{filename}.bam"
  output:
    "{filename}_sorted.bam"
  shell:
    "samtools sort {input} -o {output} -@ 8"

rule index_bam:
  input:
    "{filename}.bam"
  output:
    "{filename}.bam.bai"
  shell:
    "samtools index {input} -@ 8"

rule samtools_flagstat:
  input:
    "{filename}.bam"
  output:
    "{filename}_flagstats.txt"
  shell:
    "samtools flagstat {input} -@ 8 > {output}"

rule mark_duplicates:
  input:
    "{filename}_sorted.bam"
  output:
    bam="{filename}_refined.bam",
    metrics="{filename}_dupl_metrics.txt"
  shell:
    "picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

rule samtools_depth:
  input:
    "{filename}_refined.bam"
  output:
    "{filename}_coverage.txt"
  shell:
    "samtools depth -a {input} > {output}"

rule call_variants:
  input:
    ref="/home/jovyan/pipeline_data/GCF_000001405.26_GRCh38_genomic.fna",
    bam="{filename}_refined.bam"
  output:
    "{filename}.vcf"
  shell:
    "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -L chr10 -mbq 20 --minimum-mapping-quality 50 -O {output}"

rule filter_variants:
  input:
    "{filename}.vcf"
  output:
    "{filename}_filtered.vcf",
    "{filename}_filtered.kept.sites",
    "{filename}_filtered.removed.sites"
  shell:
    "vcftools --vcf {input} --minDP 3 --minQ 20 --recode --recode-INFO-all --stdout | vcftools --max-missing 1 --out {output} --recode --recode-INFO-all --kept-sites --removed-sites"

rule annotate_variants:
  input:
    "{filename}_filtered.vcf"
  output:
    vcf="{filename}_annotated.vcf",
    html="{filename}_snpEff_summary.html",
    txt="{filename}_snpEff_genes.txt"
  shell:
    "snpEff GRCh38 {input} -stats '/home/jovyan/pipeline_data/' > {output.vcf}"

rule filter_annotations:
  input:
    "{filename}_annotated.vcf"
  output:
    "{filename}_CYP2C19.vcf"
  shell:
    "cat {input} | SnpSift filter (ANN[*].GENE = 'CYP2C19') > {output}"

