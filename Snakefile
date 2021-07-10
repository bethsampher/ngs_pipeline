rule_all:
  input:
    "/home/joyvan/pipeline_data/SRR622461_1_fastqc.html",
    "/home/joyvan/pipeline_data/SRR622461_2_fastqc.html"
    "/home/joyvan/pipeline_data/SRR622461_1_fastp_fastqc.html",
    "/home/joyvan/pipeline_data/SRR622461_2_fastp_fastqc.html",
    "/home/joyvan/pipeline_data/SRR622461_sorted.bam.bai",
    "/home/joyvan/pipeline_data/SRR622461_stats.txt",
    "/home/joyvan/pipeline_data/SRR622461_final.bam.bai"
    "/home/joyvan/pipeline_data/SRR622461_refined_stats.txt",
    "/home/joyvan/pipeline_data/SRR622461_coverage.txt"

rule fastqc_original:
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
    rev="{filename}_2_fastp.fastq.gz"
  shell:
    "fastp --thread 8 -i {input.fwd} -o {output.fwd} -I {input.rev} -O {output.rev} --disable_adapter_trimming --length_required 36 -3 --correction"

rule fastqc_fastp:
  input:
    "{filename}_fastp.fastq.gz"
  output:
    "{filename}_fastp_fastqc.html",
    "{filename}_fastp_fastqc.zip"
  shell:
    "fastqc {input}"

rule bwa_mem
  input:
    ref="",
    fwd="{filename}_1_fastp.fastq.gz",
    rev="{filename}_2_fastp.fastq.gz"
  output:
    "{filename}.sam"
  shell:
    "bwa mem -t 8 {input.ref} {input.fwd} {input.rev} > {output}"

rule sam_to_bam
  input:
    "{filename}.sam"
  output:
    "{filename}.bam"
  shell:
    "samtools view -b {input} -o {output} -@ 8"

rule sort_bam
  input:
    "{filename}.bam"
  output:
    "{filename}_sorted.bam"
  shell:
    "samtools sort {input} -o {output} -@ 8"

rule index_sorted_bam
  input:
    "{filename}_sorted.bam"
  output:
    "{filename}_sorted.bam.bai"
  shell:
    "samtools index {input} -@ 8"

rule samtools_flagstat
  input:
    "{filename}_sorted.bam"    
  output:
    "{filename}_stats.txt"
  shell:
    "samtools flagstat {input} -@ 8 > {output}" 

rule mark_duplicates
  input:
    "{filename}_sorted.bam"
  output:
    bam="{filename}_final.bam",
    metrics="dupl_metrics.txt"
  shell:
    "picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}"

rule index_final_bam
  input:
    "{filename}_final.bam"
  output:
    "{filename}_final.bam.bai"
  shell:
    "samtools index {input} -@ 8"

rule samtools_flagstat_refined
  input:
    "{filename}_final.bam"    
  output:
    "{filename}_refined_stats.txt"
  shell:
    "samtools flagstat {input} -@ 8 > {output}" 

rule samtools_depth
  input:
    "{filename}_final.bam"    
  output:
    "{filename}_coverage.txt"
  shell:
    "samtools depth -a {input} > {output}" 

rule call_variants
  input:
    ref="",
    bam="{filename}_final.bam"
  output:
    "{filename}.vcf"
  shell:
    "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -L chr10 -mbq 20 --minimum-mapping-quality 50 -O {output}"

