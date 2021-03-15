import os
from src.common import read_samples

container: "docker://continuumio/miniconda3:4.9.2"
configfile: "config/config.yml"
samples = read_samples(config["sample_list"])

localrules: all, link, multiqc

rule all:
    "Main rule for workflow"
    input:
        "results/multiqc/multiqc.html"

rule link:
    input:
        lambda wildcards: samples[wildcards.sample][wildcards.R]
    output:
        temp("results/intermediate/{sample}_{R}.fastq.gz")
    params:
        i = lambda wildcards, input: os.path.abspath(input[0]),
        o = lambda wildcards, output: os.path.abspath(output[0])
    shell:
        """
        ln -s {params.i} {params.o}
        """

rule cutadapt:
    input:
        R1 = "results/intermediate/{sample}_R1.fastq.gz",
        R2 = "results/intermediate/{sample}_R2.fastq.gz"
    output:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz"
    log:
        "results/logs/cutdapt/{sample}.log"
    threads: 4
    params:
        R1_adapter = config["cutadapt"]["R1_adapter"],
        R2_adapter = config["cutadapt"]["R2_adapter"]
    conda:
        "envs/cutadapt.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    shell:
        """
        cutadapt -p {threads} -a {params.R1_adapter} -A {params.R2_adapter} \
            -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log} 2>&1
        """

rule fastqc:
    input:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz"
    output:
        R1 = "results/fastqc/{sample}_R1_fastqc.zip",
        R2 = "results/fastqc/{sample}_R2_fastqc.zip"
    log:
        "results/logs/fastqc/{sample}.log"
    params:
        dir = lambda wildcards, output: os.path.dirname(output.R1)
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    conda:
        "envs/fastqc.yml"
    shell:
        """
        fastqc -q --noextract -o {params.dir} {input} >{log} 2>&1
        """

rule multiqc:
    input:
        expand("results/logs/cutdapt/{sample}.log",
            sample=samples.keys()),
        expand("results/fastqc/{sample}_R{i}_fastqc.zip",
            sample=samples.keys(), i=[1,2])
    output:
        "results/multiqc/multiqc.html"
    log:
        "results/logs/multiqc/multiqc.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        base = lambda wildcards, output: os.path.basename(output[0])
    conda:
        "envs/multiqc.yml"
    shell:
        """
        multiqc -o {params.outdir} -n {params.base} {input} > {log} 2>&1
        """