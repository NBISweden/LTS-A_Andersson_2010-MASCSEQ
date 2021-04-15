import os
from snakemake.utils import validate
from src.common import read_samples, assembly_input

include: "src/common.py"
container: "docker://continuumio/miniconda3:4.9.2"
configfile: prependWfd("config/config.yml")

validate(config, schema=prependWfd("config/config_schema.yml"), set_default=True)
samples = read_samples(prependWfd(config["sample_list"]))

wildcard_constraints:
    assembler = "transabyss|trinity"

localrules: all, link, download_rna, multiqc, extractTranscriptsFromGenome

def kallisto_output(samples, config):
    files = []
    for sample, vals in samples.items():
        t = samples[sample]["type"]
        if t == "genome":
            ref = samples[sample]["reference"]
            files.append(f"results/kallisto/{sample}/{ref}.mRNA.abundance.tsv")
        elif t == "transcriptome":
            files+=[f"results/kallisto/{sample}/{assembler}.mRNA.abundance.tsv" for assembler in config["assemblers"]]
    return files

rule all:
    """Main rule for workflow"""
    input:
        "results/multiqc/multiqc.html",
        kallisto_output(samples, config)

#####################
### PREPROCESSING ###
#####################

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

################
### ASSEMBLY ###
################

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
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
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
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
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
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        max_mem=$(({params.cpumem} * {threads}))
        Trinity --seqType fq --left {input.R1} --right {input.R2} --CPU {threads} --output {params.tmpdir} --max_memory ${{max_mem}}G > {log} 2>&1
        mv {params.tmpdir}/* {params.outdir}/
        """


###############
### MAPPING ###
###############

rule extractTranscriptsFromGenome:
    input:
        fasta = lambda wc: config["genome"][wc.ref]["fasta"],
        gff = lambda wc: config["genome"][wc.ref]["gff"]
    output:
        fasta = "reference/genome/{ref}_transcriptsFromGenome.fasta.gz"
    log:
        "reference/logs/genome/{ref}_extractTranscriptsFromGenome.log"
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

rule kallisto_index_asm:
    """
    Index an assembly for kallisto mapping
    """
    input:
        fasta = assembly_input
    output:
        index = "results/kallisto/{sample}/{assembler}_transcripts.idx"
    log:
        "results/logs/kallisto/{sample}/{assembler}.index.log"
    conda:
        "envs/kallisto.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    threads: 20
    shell:
        """
        exec &> {log}

        kallisto index \
        -i {output.index} \
        {input.fasta} 

        echo "Done!"
        """


rule kallisto_index_genome:
    input:
        fasta = "reference/genome/{ref}_transcriptsFromGenome.fasta.gz"
    output:
        index = "reference/genome/{ref}_transcripts.idx"
    log:
        "reference/logs/genome/{ref}_kallisto_index.log"
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

rule kallisto_map_asm:
    input:
        R1="results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2="results/sortmerna/{sample}.{RNA}_rev.fastq.gz",
        index="results/kallisto/{sample}/{assembler}_transcripts.idx"
    output:
        tsv = "results/kallisto/{sample}/{assembler}.{RNA}.abundance.tsv",
        h5 = "results/kallisto/{sample}/{assembler}.{RNA}.abundance.h5",
        info = "results/kallisto/{sample}/{assembler}.{RNA}.run_info.json"
    log:
        "results/logs/kallisto/{sample}/{assembler}.{RNA}.kallisto_map.log"
    params:
        out = "$TMPDIR/{sample}.{assembler}.{RNA}"
    threads: 10
    resources:
        runtime= lambda wildcards,attempt: attempt ** 2 * 60
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

rule kallisto_map_genome:
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq.gz",
        index = "reference/genome/{ref}_transcripts.idx"
    output:
        tsv = "results/kallisto/{sample}/{ref}.{RNA}.abundance.tsv",
        h5 = "results/kallisto/{sample}/{ref}.{RNA}.abundance.h5",
        info = "results/kallisto/{sample}/{ref}.{RNA}.run_info.json"
    log:
        "results/logs/kallisto/{sample}/{ref}.{RNA}.kallisto_map.log"
    params:
        out = "$TMPDIR/{sample}.{ref}.{RNA}"
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


rule star_index_genome:
    """
    Creates a STAR index file from a gzipped fasta file with either:
    - reference genome sequences, e.g., for chromosomes or contigs
    - transcript sequences from a de novo transcriptome assembly or 
      extracted from a reference genome (not yet implemented)
    """
    input:
        fasta = "results/{reftype}/reference//{ref}.fasta.gz"
    output:
        index = "results/{reftype}/reference/star/{ref}.idx"
    log:
        "results/logs/{reftype}/{ref}_star_index.log"
    params:
        genomedir = "results/{reftype}/reference/star/",
        sjdbOverhang = int(samples[sample]["read_length"]) - 1 # read length -1
    conda:
        "envs/star.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60 
    threads: 1
    shell:
        """
        exec &> {log}

        ## get some parameter values from params and input files and use these to set
        ## recommended optionparameter values for STAR options
        genomelength=$(cat {input.fasta} | awk 'BEGIN{ret=0} !/\>/ {ret=ret+length($0)} END{print ret}')
        nseqs=$(cat {input.fasta} | grep -c ">")
        
        # This should be min( 14, log2[ $genomelength ] / 2 - 1 )
        genomeSAindexNbases = echo $genomelength| awk '{s = int(log($0)/2log(2) - 1); if(s>14){s=14}; print s}'

        # This should be min( 18, log2[ max( $genomelength / $nseqs, {params.readlength} ) ] )
        genomeChrBinNbits = echo $genomelength $nseqs {params.readlength} | \ 
        awk '{ s = $1 / $2; if(s < $3){ s = $3 }; s = log(s) / log(2); if(s > 18){ s = 18 };print s }'

        ## Start STAR
        star \
        --runMode genomeGenerate \
        --genomeDir {output.index} \
        --genomeFastaFiles{input.fasta} \
        --runThreadN {threads} \
        # Splice junction option
        --sjdbGTFfile {input.gtf} \                        # Annotation file (for genome ref)
        --sjdbOverhang {params.sjdbOverhang} \             # nbases to use for splice junction signature
        --sjdbGTFfeatureExon exon \                        # What main feature to use from gff/gtf
        --sjdbGTFtagExonParentTranscript Parent \          # Parent main feature field name gff3
        #--sjdbFileChrStartEnd path/to/splicejunctionfile \ # Alternative Splice junction file (not used here)
        # Options reducing computational load
        --genomeSAindexNbases $genomeSAindexNbases \       # Very small genomes
        --genomeChrBinNbits  $genomeChrBinNbits            # Very many individual sequences (e.g.,transcriptome)

        echo "Done!"
        """

rule star_map:
    """
    Under construction
    """
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq.gz",
        index = "reference/{reftype}/star/{ref}_transcripts.idx"
    output:
        bam = "results/star/{sample}/{ref}.{RNA}.Aligned.out.bam",
        logout = "results/star/{sample}/{ref}.{RNA}.Log.out",
        logfinal = "results/star/{sample}/{ref}.{RNA}.Log.final.out",
    log:
        "results/logs/kallisto/{sample}/{ref}.{RNA}.star_map.log"
    params:
        genomedir = "reference/genome/star/",
        outprefix = "results/star/{sample}/{ref}.{RNA}"
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt ** 2 * 60
    conda:
        "envs/star.yml"
    shell:
        """
        exec &> {log}

        star \
        --genomeDir {params.genomedir} \
        --readFilesIn {input.R1},{input.R2} \
        --outFileNamePrefix {params.outprefix} \
        #Adjustable in st_pipeline
        --outFilterMultimapNmax 1 \   # option disable_multimap == True -> 1 (deafult 20) 
        --runThreadN {threads} \      # option threads (default 8)
        #--clip3pNbases 0 \            # option inverse_mapping_rv_trimming (default 0)
        #--clip5pNbases 0 \            # option mapping_rv_trimming (default 0)
        #--alignEndsType EndToEnd \    # option disable_clipping (default EndToEnd)
        #--alignIntronMin 1 \          # option min_intron_size (default 1)
        #--alignIntronMax 1 \          # option max_intron_size (default 20)
        #--outFilterMatchNmin 20 \     # option min_length_trimming (default 20)
        #--genomeLoad NoSharedMemory \ # option star_genome_loading (default NoSharedMemory)
        #--limitBAMsortRAM 0 \         # option star_sort_mem_limit (default 0)
        # Hardcoded by st_pipeline non-default STAR
        --outSAMtype BAM SortedByCoordinate \   # (STAR default: SAM)
        --outSAMmultNmax 1 \                    # (STAR default: -1)
        --outMultimapperOrder Random \          # (STAR default: Old 2.4)
        --outFilterMismatchNoverLmax 0.1 \      # (STAR default: 0.3)
        --readFilesType SAM SE \                # (STAR default: SAM)
        --readFilesCommand samtools view -h \   # (STAR default: Fastx)
        # Hardcoded by st_pipeline default STAR
        #--outFilterType Normal \
        #--outSAMorder Paired \
        #--outSAMprimaryFlag OneBestScore \
        #--readMatesLengthsIn NotEqual \


        echo "Done!"
        """
        
        
##################
### ANNOTATION ###
##################

rule transrate:
    input:
        R1="results/sortmerna/{sample}.mRNA_fwd.fastq.gz",
        R2="results/sortmerna/{sample}.mRNA_rev.fastq.gz",
        fa=assembler_input
    output:
        "results/transrate/{sample}/{assembler}.csv",
        directory("results/transrate/{sample}/{assembler}")
    params:
        outdir = "$TMPDIR/transrate-{sample}.{assembler}",
    threads: 10
    conda:
        "envs/transrate.yml"
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        transrate --assembly {input.fa} --left {input.R1} --right {input.R2} \
            --threads {threads} --output {params.outdir}
        mv {params.outdir}/assemblies.csv {output[0]}
        mv {params.outdir}/merged {output[1]}
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
        assembly_input,
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
