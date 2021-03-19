import os
from snakemake.utils import validate
from src.common import read_samples, dammit_input

container: "docker://continuumio/miniconda3:4.9.2"
configfile: "config/config.yml"

validate(config, schema="config/config_schema.yml", set_default=True)
samples = read_samples(config["sample_list"])

wildcard_constraints:
    assembler = "transabyss|trinity"

localrules: all, link, download_rna, multiqc

rule all:
    """Main rule for workflow"""
    input:
        "results/multiqc/multiqc.html",
        expand("results/transabyss/{sample}/merged.fa",
            sample = samples.keys()),
        expand("results/trinity/{sample}/Trinity.fasta",
            sample = samples.keys())

rule link:
    """Links input files so that naming is consistent within workflow"""
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
    """Runs cutadapt on raw input files"""
    input:
        R1 = ancient("results/intermediate/{sample}_R1.fastq.gz"),
        R2 = ancient("results/intermediate/{sample}_R2.fastq.gz")
    output:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz"
    log:
        "results/logs/cutadapt/{sample}.log"
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
        cutadapt -j {threads} -a {params.R1_adapter} -A {params.R2_adapter} \
            -o {output.R1} -p {output.R2} {input.R1} {input.R2} > {log} 2>&1
        """

rule download_rna:
    """Downloads RNA databases from SortMeRNA GitHub repo. 
    
    The db wildcard matches any *.fasta file, it gets expanded in the actual
    sortmerna rule below.
    """
    output:
        "resources/sortmerna/{db}.fasta"
    params:
        url_base = "https://raw.githubusercontent.com/biocore/sortmerna/master/data/rRNA_databases/",
        basename = lambda wildcards, output: os.path.basename(output[0]),
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    log:
        "resources/sortmerna/{db}.log"
    shell:
        """
        mkdir -p {params.outdir}
        curl -L -o {output[0]} {params.url_base}/{params.basename} > {log} 2>&1
        """

rule sortmerna:
    """
    Runs sortmerna on the output from cutadapt.
    
    Here the db wildcard is expanded using all RNA files set in the config
    for sortmerna. The rule produces four output files:
    - fwd and rev files for reads aligned to any of the rRNA databases (rRNA)
    - fwd and rev files for reads not aligned to any of the databases (mRNA)
    Note that although we name the non-aligned fraction 'mRNA' it is in fact 
    just 'non_rRNA'. The 'paired_out' strategy ensures that both reads for a 
    given pair end up in the mRNA fraction. 
    """
    input:
        R1 = "results/cutadapt/{sample}_R1.fastq.gz",
        R2 = "results/cutadapt/{sample}_R2.fastq.gz",
        db = expand("resources/sortmerna/{db}.fasta",
            db = config["sortmerna"]["dbs"])
    output:
        rR1 = "results/sortmerna/{sample}.rRNA_fwd.fastq.gz",
        rR2 = "results/sortmerna/{sample}.rRNA_rev.fastq.gz",
        mR1 = "results/sortmerna/{sample}.mRNA_fwd.fastq.gz",
        mR2 = "results/sortmerna/{sample}.mRNA_rev.fastq.gz"
    log:
        runlog="results/logs/sortmerna/{sample}.log",
        reportlog="results/sortmerna/{sample}.log"
    params:
        string = lambda wildcards, input: " ".join([f"--ref {x}" for x in input.db]),
        workdir = "$TMPDIR/sortmerna/{sample}.wd",
        outdir = lambda wildcards, output: os.path.dirname(output.rR1)
    threads: 10
    conda:
        "envs/sortmerna.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 48
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        rm -rf {params.workdir}
        mkdir -p {params.workdir}
        sortmerna --threads {threads} --workdir {params.workdir} --fastx \
            --reads {input.R1} --reads {input.R2} {params.string} --paired_out \
            --out2 --aligned {params.workdir}/{wildcards.sample}.rRNA \
            --other {params.workdir}/{wildcards.sample}.mRNA > {log.runlog} 2>&1
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.rRNA.log {log.reportlog}
        rm -rf {params.workdir}
        """

rule fastqc:
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq.gz"
    output:
        R1 = "results/fastqc/{sample}.{RNA}_fwd_fastqc.zip",
        R2 = "results/fastqc/{sample}.{RNA}_rev_fastqc.zip"
    log:
        "results/logs/fastqc/{sample}.{RNA}.log"
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
        expand("results/logs/cutadapt/{sample}.log",
            sample=samples.keys()),
        expand("results/sortmerna/{sample}.log",
            sample=samples.keys()),
        expand("results/fastqc/{sample}.{RNA}_{R}_fastqc.zip",
            sample=samples.keys(), R=["rev","fwd"], RNA=["mRNA","rRNA"])
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

rule transabyss:
    input:
        R1 = "results/sortmerna/{sample}.mRNA_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.mRNA_rev.fastq.gz"
    output:
        "results/transabyss/{sample}/{k}/{sample}.{k}-final.fa",
        "results/transabyss/{sample}/{k}/coverage.hist",
    log:
        "results/logs/transabyss/{sample}.{k}.log"
    conda:
        "envs/transabyss.yml"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = "$TMPDIR/{sample}.transabyss"
    threads: config["transabyss"]["threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 96
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        transabyss --pe {input.R1} {input.R2} -k {wildcards.k} \
            --outdir {params.tmpdir} --name {wildcards.sample}.{wildcards.k} \
            --threads {threads} >{log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        """

rule transabyss_merge:
    input:
        expand("results/transabyss/{{sample}}/{k}/{{sample}}.{k}-final.fa",
            k = config["transabyss"]["kmer"])
    output:
        "results/transabyss/{sample}/merged.fa"
    log:
        "results/logs/transabyss/{sample}.merge.log"
    conda:
        "envs/transabyss.yml"
    params:
        mink = min(config["transabyss"]["kmer"]),
        maxk = max(config["transabyss"]["kmer"]),
        i = lambda wildcards, input: sorted(input),
        prefix = [f"k{x}." for x in sorted(config["transabyss"]["kmer"])],
        tmpout = "$TMPDIR/{sample}.ta.merged.fa"
    threads: config["transabyss"]["threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 10
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        transabyss-merge {params.i} --mink {params.mink} --maxk {params.maxk} \
            --out {params.tmpout} --threads {threads} --prefix {params.prefix} > {log} 2>&1
        mv {params.tmpout} {output}
        """

rule trinity:
    input:
        R1="results/sortmerna/{sample}.mRNA_fwd.fastq.gz",
        R2="results/sortmerna/{sample}.mRNA_rev.fastq.gz"
    output:
        "results/trinity/{sample}/Trinity.fasta"
    log:
        "results/logs/trinity/{sample}.log"
    conda:
        "envs/trinity.yml"
    threads: config["trinity"]["threads"]
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        tmpdir = "$TMPDIR/{sample}.trinity",
        cpumem = config["mem_per_cpu"]
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 120
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        max_mem=$(({params.cpumem} * {threads}))
        Trinity --seqType fq --left {input.R1} --right {input.R2} --CPU {threads} --output {params.tmpdir} --max_memory ${{max_mem}}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        """

rule dammit_db:
    output:
        touch("resources/dammit/done")
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    conda:
        "envs/dammit.yml"
    shell:
        """
        dammit databases --database-dir {params.outdir} --install
        """

rule dammit:
    input:
        dammit_input
    output:
        touch("results/dammit/{sample}/{assembler}/done")
    log:
        "results/logs/dammit/{sample}.{assembler}.log"
    conda:
        "envs/dammit.yml"
    shell:
        """
        dammit -v > {log} 2>&1
        """
