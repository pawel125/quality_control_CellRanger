import pandas as pd
import os


ref_genomes = dict(zip(samples.sample_id, samples.ref_genome))
ref_bed_files = {
    sample: config["rseqc_bed_file"][genome] 
    for sample, genome 
    in ref_genomes.items()}

# rule all:
    # input: "multiqc_report.html"

rule multiqc:
    input:
        fastqc=expand("qc/fastqc/{sample}_fastqc.html", sample = samples.sample_id),
        read_GC=expand("qc/rseqc/{sample}.read_GC.GC.xls", sample = samples.sample_id),
        bam_stat=expand("qc/rseqc/{sample}.bam_stat", sample = samples.sample_id),
        tin=expand("qc/rseqc/{sample}.tin.xls", sample = samples.sample_id),
        genebody_ceverage=expand("qc/rseqc/{sample}.geneBodyCoverage.txt", sample = samples.sample_id)
    output: "multiqc_report.html"
    shell: "multiqc qc"

rule mklink:
    input: 
        bam="single_cell_data/{sample}/outs/possorted_genome_bam.bam",
        bai="single_cell_data/{sample}/outs/possorted_genome_bam.bam.bai"
    output: 
        bam=temp("qc/{sample}.bam"),
        bai=temp("qc/{sample}.bam.bai")
    shell: 
        """
        mkdir -p qc
        ln {input.bam} {output.bam}
        ln {input.bai} {output.bai}
        """

rule fastqc:
    input: "qc/{sample}.bam"
    output: "qc/fastqc/{sample}_fastqc.html"
    params:
        outdir = "qc/fastqc"
    threads: 20
    shell: 
        """
        mkdir -p {params.outdir}
        srun fastqc -o {params.outdir} -t {threads} {input}
        """

rule read_GC:
    input: 
        bam="qc/{sample}.bam", 
        bai="qc/{sample}.bam.bai"
    output: "qc/rseqc/{sample}.read_GC.GC.xls"
    conda: "../envs/rseqc_env.yml"
    params:
        out_prefix="qc/rseqc/{sample}.read_GC"
    shell: "srun read_GC.py -i {input.bam} -o {params.out_prefix}"


rule bam_stat:
    input: 
        bam="qc/{sample}.bam", 
        bai="qc/{sample}.bam.bai"
    output: "qc/rseqc/{sample}.bam_stat"
    conda: "../envs/rseqc_env.yml"
    shell: "bam_stat.py  -i {input.bam} > {output}"

rule tin:
    input: ["qc/{sample}.bam", "qc/{sample}.bam.bai"]
    output: 
        summary="qc/rseqc/{sample}.summary.txt",
        tin="qc/rseqc/{sample}.tin.xls"
    conda: "../envs/rseqc_env.yml"
    params:
        bam="{sample}.bam",
        bed=lambda wildcards: ref_bed_files[wildcards.sample]
    shell: 
        """
        cd qc
        srun tin.py -i {params.bam} -r {params.bed}
        mv {wildcards.sample}.summary.txt rseqc/{wildcards.sample}.summary.txt
        mv {wildcards.sample}.tin.xls rseqc/{wildcards.sample}.tin.xls
        cd ..
        """
    
rule genebody_coverage:
    input: 
        bam="qc/{sample}.bam", 
        bai="qc/{sample}.bam.bai"
    output: "qc/rseqc/{sample}.geneBodyCoverage.txt"
    conda: "../envs/rseqc_env.yml"
    log: "outputs/genebody_coverage.{sample}.log"
    params:
        bam="../{sample}.bam",
        bed=lambda wildcards: ref_bed_files[wildcards.sample],
        tempdir = "qc/temp_genebody_coverage_{sample}"
    shell: 
        """
        mkdir -p {params.tempdir}
        cd {params.tempdir}
        srun geneBody_coverage.py -i {params.bam} -r {params.bed} -o {wildcards.sample}
        cp {wildcards.sample}.geneBodyCoverage.* ../rseqc/
        rm -rf {params.tempdir}
        """

