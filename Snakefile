import os
from snakemake.utils import validate

include: "src/common.py"
container: "docker://continuumio/miniconda3:4.9.2"
configfile: prependWfd("config/config.yml")

validate(config, schema=prependWfd("config/config_schema.yml"), set_default=True)
samples = read_samples(prependWfd(config["sample_list"]))

wildcard_constraints:
    assembler = "transabyss|trinity"

localrules: all, link, download_rna, multiqc

def kallisto_output(samples):
    files = []
    for sample, vals in samples.items():
        if vals["type"] == "genome":
            files.append(f"results/kallisto/{sample}.mRNA.abundance.tsv")
    return files

rule all:
    """Main rule for workflow"""
    input:
        "results/multiqc/multiqc.html",
        expand("results/transabyss/{sample}/merged.fa",
            sample = samples.keys()),
        expand("results/trinity/{sample}/Trinity.fasta",
            sample = samples.keys()),
        kallisto_output(samples)

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
        R2_adapter = config["cutadapt"]["R2_adapter"],
        minlen = config["cutadapt"]["minlen"]
    conda:
        "envs/cutadapt.yml"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 2
    shell:
        """
        cutadapt -m {params.minlen} -j {threads} -a {params.R1_adapter} \
            -A {params.R2_adapter} -o {output.R1} -p {output.R2} \
            {input.R1} {input.R2} > {log} 2>&1
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
    just 'non_rRNA'. The 'paired_in' strategy means that if at least one read in
    a pair is flagged as rRNA, both are placed in the 'rRNA' fraction.
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
            --reads {input.R1} --reads {input.R2} {params.string} --paired_in \
            --out2 --aligned {params.workdir}/{wildcards.sample}.rRNA \
            --other {params.workdir}/{wildcards.sample}.mRNA > {log.runlog} 2>&1
        gzip {params.workdir}/*.fastq
        mv {params.workdir}/*.gz {params.outdir}
        mv {params.workdir}/{wildcards.sample}.rRNA.log {log.reportlog}
        rm -rf {params.workdir}
        """

rule extractTranscriptsFromGenome:
    input:
        fasta = lambda wc: config["genome"][wc.ref]["fasta"],
        gff = lambda wc: config["genome"][wc.ref]["gff"]
    output:
        fasta = "reference/{ref}_transcriptsFromGenome.fasta.gz"
    log:
        "reference/logs/{ref}_extractTranscriptsFromGenome.log"
    params:
        script = prependWfd("scripts/fixTranscriptId.py"),
        make_me_local = True
    conda:
        "envs/gffread.yaml"
    threads: 1
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    shell:
        """
        exec &> {log}        

        gffread {input.gff} -g {input.fasta} -w {output.fasta} -E -O
        # Additional mRNA-specific options -C -V -M 

        echo "Done!"
        """

rule kallisto_index:
    input:
        fasta = "reference/{ref}_transcriptsFromGenome.fasta.gz" \
            if samples["ref"]["type"] == "genome" else \
            "reference/{ref}_transcripts.fasta.gz"
    output:
        index = "reference/{ref}_transcripts.idx"
    log:
        "reference/logs/{ref}_kallisto_index.log"
    conda:
        "envs/kallisto.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    threads: 1
    shell:
        """
        exec &> {log}

        kallisto index \
        -i {output.index} \
        {input.fasta} 

        echo "Done!"
        """
        
rule kallisto_map:
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq.gz",
        index = lambda wc: expand("reference/{ref}_transcripts.idx",
                                  ref = samples[wc.sample]["reference"])
    output:
        tsv = "results/kallisto/{sample}.{RNA}.abundance.tsv",
        h5 = "results/kallisto/{sample}.{RNA}.abundance.h5",
        info = "results/kallisto/{sample}.{RNA}.run_info.json"
    log:
        "results/logs/{sample}.{RNA}_kallisto_map.log"
    params:
        out = "results/kallisto/"
    threads: 10
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    conda:
        "envs/kallisto.yml"
    shell:
        """
        exec &> {log}

        kallisto quant \
        -i {input.index} \
        -o {params.out} \
        -t {threads} \
        {input.R1} {input.R2}

        # Change to informative file names
        mv {params.out}/abundance.tsv {output.tsv}
        mv {params.out}/abundance.h5 {output.h5}
        mv {params.out}/run_info.json {output.info}

        echo "Done!"
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

rule dammit_busco:
    output:
        directory("resources/dammit/busco2db/{busco_group}_ensembl"),
        touch("resources/dammit/busco2db/download_and_untar:busco2db-{busco_group}.done")
    params:
        tmpdir = "$TMPDIR/dammit/{busco_group}"
    shell:
        """
        mkdir -p {params.tmpdir}
        curl -L https://busco.ezlab.org/v2/datasets/{wildcards.busco_group}_ensembl.tar.gz | tar -xz -C {params.tmpdir}
        mv {params.tmpdir}/* {output[0]} 
        """

rule dammit:
    input:
        dammit_input,
        lambda wildcards: f"resources/dammit/busco2db/download_and_untar:busco2db-{config['dammit']['busco_group'][wildcards.sample]}.done"
    output:
        touch("results/dammit/{sample}/{assembler}/done")
    log:
        "results/logs/dammit/{sample}.{assembler}.log"
    params:
        dbdir = "resources/dammit",
        outdir = lambda wildcards, output: os.path.dirname(output[0]),
        busco_group = lambda wildcards: config["dammit"]["busco_group"][wildcards.sample],
        tmpdir = "$TMPDIR/dammit/{sample}.{assembler}"
    conda:
        "envs/dammit.yml"
    threads: 10
    shell:
        """
        dammit annotate {input[0]} -n {wildcards.sample}.{wildcards.assembler} \
            --n_threads {threads} --database-dir {params.dbdir} \
            -o {params.tmpdir} --force --busco-group {params.busco_group} \
            --quick --verbosity 2 > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}
        """
