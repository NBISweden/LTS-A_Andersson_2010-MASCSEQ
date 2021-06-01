import os
from snakemake.utils import validate
from src.common import read_samples, assembly_input

include: "src/common.py"
container: "docker://continuumio/miniconda3:4.9.2"
configfile: prependWfd("config/config.yml")

validate(config, schema=prependWfd("config/config_schema.yml"), set_default=True)
samples = read_samples(config["sample_list"])

wildcard_constraints:
    assembler = "transabyss|trinity"

localrules:
    all,
    link,
    download_rna,
    multiqc,
    linkReferenceGenome,
    extractTranscriptsFromGenome,
    gunzipReads,
    busco_dl,
    busco

def busco_input(samples, config):
    files = []
    for sample, lineage in config["busco"]["sample_lineages"].items():
        if samples[sample]["type"] == "transcriptome":
            for assembler in config["assemblers"]:
                files.append(f"results/busco/{assembler}/{sample}/short_summary.specific.{lineage}.{sample}.txt")
    return files

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

def busco_input(samples, config):
    files = []
    for sample, lineage in config["busco"]["sample_lineages"].items():
        if samples[sample]["type"] == "transcriptome":
            for assembler in config["assemblers"]:
                files.append(f"results/busco/{assembler}/{sample}/short_summary.specific.{lineage}.{sample}.txt")
    return files

rule all:
    """Main rule for workflow"""
    input:
        "results/multiqc/multiqc.html",
        busco_input(samples, config)
        #kallisto_output(samples, config)

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
        tmpdir = "$TMPDIR/{sample}.transabyss",
        R1 = "$TMPDIR/{sample}.transabyss/R1",
        R2 = "$TMPDIR/{sample}.transabyss/R2"
    threads: config["transabyss"]["threads"]
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 150
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        transabyss --pe {params.R1} {params.R2} -k {wildcards.k} \
            --outdir {params.tmpdir} --name {wildcards.sample}.{wildcards.k} \
            --threads {threads} >{log} 2>&1
        rm {params.R1} {params.R2}
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

def trinity_strand_string(wildcards):
    try:
        strandness = samples[wildcards.sample]["strandness"]
    except KeyError:
        return ""
    if strandness == "sense":
        return "--SS_lib_type FR"
    elif strandness == "antisense":
        return "--SS_lib_type RF"

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
        ss_lib_type = trinity_strand_string,
        tmpdir = "$TMPDIR/{sample}.trinity",
        cpumem = config["mem_per_cpu"],
        R1 = "$TMPDIR/{sample}.trinity/R1",
        R2 = "$TMPDIR/{sample}.trinity/R2"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 48
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        gunzip -c {input.R1} > {params.R1}
        gunzip -c {input.R2} > {params.R2}
        max_mem=$(({params.cpumem} * {threads}))
        Trinity --seqType fq {params.ss_lib_type} --left {params.R1} --right {params.R2} --CPU {threads} --output {params.tmpdir} --max_memory ${{max_mem}}G > {log} 2>&1
        rm {params.R1} {params.R2}
        mv {params.tmpdir}/* {params.outdir}/
        """


###############
### MAPPING ###
###############

rule linkReferenceGenome:
    input:
        fasta = lambda wc: config["genome"][wc.ref]["fasta"],
        gff = lambda wc: config["genome"][wc.ref]["gff"]
    output:
        fasta = "results/genome/reference/{ref}.fasta.gz",
        gff = "results/genome/reference/{ref}.gff"
    log:
        "results/logs/genome/reference/{ref}_linkReferenceGenome.log"
    params:
        fasta = prependPwd("results/genome/reference/{ref}.fasta.gz"),
        gff = prependPwd("results/genome/reference/{ref}.gff")
    shell:
        """
        exec &> {log}        
        
        ln -s  {input.fasta} {params.fasta}
        ln -s  {input.gff} {params.gff}

        """
        
        
rule extractTranscriptsFromGenome:
    input:
        fasta = "results/genome/reference/{ref}.fasta.gz",
        gff = "results/genome/reference/{ref}.gff"
    output:
        fasta = "results/transcriptome/reference/{ref}_transcriptsFromGenome.fasta.gz"
    log:
        "results/logs/genome/reference/{ref}_extractTranscriptsFromGenome.log"
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

# Mapping with kallisto
#   - works only with transcriptome reference files
#   - handles multimappers
#   - uses hashed pseudo-alignment approach (fast)
#   - outputs read abundance/counts directly (no extra program)
###############################################################
rule kallisto_index:
    """ Create an index for kallisto from a transcriptome (fasta file
    with transcript sequences). """
    """
    Index an assembly for kallisto mapping
    """
    input:
        fasta = "resources/{reftype}/{ref}.fasta.gz"
    output:
        index = "resources/{reftype, transcriptome.*}/kallisto/{ref}.idx"
    log:
        "resources/logs/{reftype}/kallisto/{ref}.index.log"
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
    """ Map PE RNAseq reads to indexed transcriptome using kallisto.
    result in read abundance/counts output."""
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq.gz",
        index = "resources/{reftype}/kallisto/{ref}.idx"
    output:
        tsv = "results/{reftype, transcriptome.*}/kallisto/{ref}/{sample}.{RNA}.abundance.tsv",
        h5 = "results/{reftype}/kallisto/{ref}/{sample}.{RNA}.abundance.h5",
        info = "results/{reftype}/kallisto/{ref}/{sample}.{RNA}.run_info.json"
    log:
        "results/logs/{reftype}/kallisto/{ref}/{sample}.{RNA}.kallisto_map.log"
    params:
        out = "$TMPDIR/{ref}.{sample}.{RNA}"
    threads: 10
    resources:
        runtime = lambda wildcards,attempt: attempt ** 2 * 60
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

rule star_index_transcriptome:
    """Creates a STAR index file from a gzipped fasta file of transcript sequences
    from a de novo transcriptome assembly or extracted from a reference genome.
    STAR options are selected to align with those in ST-pipeline."""
    input:
        fasta = "resources/{reftype}/{ref}.fasta.gz"
    output:
        index = directory("resources/{reftype, transcriptome.*}/star/{ref}.idx")
    log:
        "resources/logs/{reftype}/star/{ref}_star_index.log"
    params:
        genomedir = "resources/{reftype}/star/",
        readlength = 100 # not solved yet: int(samples["\{sample\}"]["read_length"]) # read length
    conda:
        "envs/star.yml"
    shell:
        """
        exec &> {log}

        ## get some parameter values from params and input files and use these to set
        ## recommended optionparameter values for STAR options
        sjdbOverhang=100 # avoid error unbound var; actually set below
        let sjdbOverhang={params.readlength}-1

        genomelength=$( cat {input.fasta} | \
           awk 'BEGIN{{ret=0}} !/>/ {{ret=ret+length($0)}} END{{print ret}}' )

        nseqs=$(cat {input.fasta} | grep -c ">")

        # This should be min( 14, log2[ $genomelength ] / 2 - 1 )
        genomeSAindexNbases=$( echo $genomelength| \
           awk '{{ s = int(log($1)/2log(2) - 1); if(s>14){{ s=14 }}; print s }}' )

        # This should be min( 18, log2[ max( $genomelength / $nseqs, {params.readlength} ) ] )
        genomeChrBinNbits=$( echo $genomelength $nseqs {params.readlength} | awk \
                             '{{ s = $1 / $2; if(s < $3){{ s = $3 }}; s = log(s) / log(2); \
                              if(s > 18){{ s = 18 }};print s }}' )

        ## Start STAR, backticks are used to allow comments in multi-line bash command
        STAR \
        --runMode genomeGenerate \
        --genomeDir {output.index} \
        --genomeFastaFiles {input.fasta} \
        --runThreadN {threads} \
        `# Options reducing computational load` \
        --genomeSAindexNbases $genomeSAindexNbases        `# Very small genomes` \
        --genomeChrBinNbits  $genomeChrBinNbits            `# Very many individual sequences (e.g.,transcriptome)`

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
        fasta = "results/{reftype}/reference/{ref}.fasta.gz",
        gff = "results/{reftype}/reference/{ref}.gff"
    output:
        index = directory("results/{reftype}/reference/star/{ref}.idx")
    log:
        "results/logs/{reftype}/{ref}_star_index.log"
    params:
        genomedir = "results/{reftype}/reference/star/",
        readlength = 100 # not solved yet: int(samples["\{sample\}"]["read_length"]) # read length
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
        sjdbOverhang=100 # avoid error unbound var; actually set below
        let sjdbOverhang={params.readlength}-1

        genomelength=$( cat {input.fasta} | \
           awk 'BEGIN{{ret=0}} !/>/ {{ret=ret+length($0)}} END{{print ret}}' )

        nseqs=$(cat {input.fasta} | grep -c ">")
        
        # This should be min( 14, log2[ $genomelength ] / 2 - 1 )
        genomeSAindexNbases=$( echo $genomelength| \
           awk '{{ s = int(log($1)/2log(2) - 1); if(s>14){{ s=14 }}; print s }}' )

        # This should be min( 18, log2[ max( $genomelength / $nseqs, {params.readlength} ) ] )
        genomeChrBinNbits=$( echo $genomelength $nseqs {params.readlength} | awk \
                             '{{ s = $1 / $2; if(s < $3){{ s = $3 }}; s = log(s) / log(2); \
                              if(s > 18){{ s = 18 }};print s }}' )

        ## Start STAR, backticks are used to allow comments in multi-line bash command
        STAR \
        --runMode genomeGenerate \
        --genomeDir {output.index} \
        --genomeFastaFiles {input.fasta} \
        --runThreadN {threads} \
        `# Splice junction option` \
        --sjdbGTFfile {input.gff}                         `# Annotation file (for genome ref)` \
        --sjdbOverhang $sjdbOverhang                      `# nbases to use for splice junction signature` \
        --sjdbGTFfeatureExon exon                         `# What main feature to use from gff/gtf` \
        --sjdbGTFtagExonParentTranscript Parent           `# Parent main feature field name gff3` \
        `#--sjdbFileChrStartEnd path/to/splicejunctionfile  # Alternative Splice junction file (not used here)` \
        `# Options reducing computational load` \
        --genomeSAindexNbases $genomeSAindexNbases        `# Very small genomes` \
        --genomeChrBinNbits  $genomeChrBinNbits            `# Very many individual sequences (e.g.,transcriptome)`

        echo "Done!"
        """

rule gunzipReads:
    output:
        fastq = temp("results/{prefix}.fastq")
    input:
        fastq = "results/{prefix}.fastq.gz"
    shell:
        """
        gunzip -c {input.fastq} > {output.fastq}
        """

rule star_map:
    """ Map PE RNAseq reads to indexed transcriptome or genome using STAR.
    No estimation of read abundance/couots is deon in this step.
    STAR options are selected to align with those in ST-pipeline."""
    input:
        R1 = "results/sortmerna/{sample}.{RNA}_fwd.fastq",
        R2 = "results/sortmerna/{sample}.{RNA}_rev.fastq",
        index = "resources/{reftype}/star/{ref}.idx"
    output:
        bam = "results/{reftype}/star/{ref}/{sample}.{RNA}.Aligned.sortedByCoord.out.bam",
#        logout = "results/{reftype}/star/{ref}/{sample}.{RNA}.Log.out",
        logfinal = "results/{reftype}/star/{ref}/{sample}.{RNA}.Log.final.out",
        SJ = "results/{reftype}/star/{ref}/{sample}.{RNA}.SJ.out.tab",
    log:
        "results/logs/{reftype}/star/{ref}/{sample}.{RNA}.star_map.log"
    params:
        outprefix = "results/{reftype}/star/{ref}/{sample}.{RNA}."
    threads: 8
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 360
    conda:
        "envs/star.yml"
    shell:
        """
        exec &> {log}
        
        STAR \
        --genomeDir {input.index} \
        --readFilesIn {input.R1},{input.R2} \
        --outFileNamePrefix {params.outprefix} \
        --outFilterMultimapNmax 1 \
        --runThreadN {threads} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax -1 \
        --outMultimapperOrder Random \
        --outFilterMismatchNoverLmax 0.1 \
        --readFilesType Fastx \
        --limitBAMsortRAM 3826372101

        ##
        # STAR options explained\
        # --genomeDir {input.index} \
        # --readFilesIn {input.R1},{input.R2} \
        # --outFileNamePrefix {params.outprefix} \
        # `#Adjustable in st_pipeline` \
        # --outFilterMultimapNmax 1    `# option disable_multimap == True -> 1 (deafult 20)` \
        # --runThreadN {threads}       `# option threads (default 8)` \
        # `#--clip3pNbases 0             # option inverse_mapping_rv_trimming (default 0)` \
        # `#--clip5pNbases 0             # option mapping_rv_trimming (default 0)` \
        # `#--alignEndsType EndToEnd     # option disable_clipping (default EndToEnd)` \
        # `#--alignIntronMin 1           # option min_intron_size (default 1)` \
        # `#--alignIntronMax 1           # option max_intron_size (default 20)` \
        # `#--outFilterMatchNmin 20      # option min_length_trimming (default 20)` \
        # `#--genomeLoad NoSharedMemory  # option star_genome_loading (default NoSharedMemory)` \
        # `#--limitBAMsortRAM 0          # option star_sort_mem_limit (default 0)` \
        # `# Hardcoded by st_pipeline non-default STAR` \
        # --outSAMtype BAM SortedByCoordinate `# (STAR default: SAM)` \
        # --outSAMmultNmax 1                  `# (st default: 1 one multimapper alignment reported; STAR default: -1=all multimappers reported)` \
        # --outMultimapperOrder Random        `# (STAR default: Old 2.4)` pick best alignemnt by random\
        # --outFilterMismatchNoverLmax 0.1    `# (STAR default: 0.3)` \
        # --readFilesType Fastx               `# (st default: SAM SE) (STAR default: fastx)` \
        # `#--readFilesCommand -                # (st default: samtools view -h)` \
        # `# Limits recommended when running on UPPMAX (by STAR error message)` \
        # --limitBAMsortRAM 3826372101 `# max avail RAM for BAM sorting (STAR default 0)
        # #for debugging
        # # --readMapNumber 1 
        # # Hardcoded by st_pipeline default STAR
        # #--outFilterType Normal \
        # #--outSAMorder Paired \
        # #--outSAMprimaryFlag OneBestScore \
        # #--readMatesLengthsIn NotEqual \
        # #--limitBAMsortRAM 3826372101 `# max avail RAM for BAM sorting (STAR default 0)


        ls {params.outprefix}*
        rm -f *.sam
        rm -f *.progress.out
        echo "Done!"
        """

###################
### READ COUNTS ###
###################

rule htseq:
    input:
        bam = "results/{reftype}/star/{ref}/{sample}.{RNA}.Aligned.sortedByCoord.out.bam",
        gff = "resources/{reftype}/{ref}.gff"
    output:
        bam = "results/{reftype, genome.*}/star/{ref}/{sample}.{RNA}.annotated.bam",
        discarded = "results/{reftype}/star/{ref}/{sample}.{RNA}.annotated_discarded.bam",
    params:
        strandness = lambda wc: "yes" if samples[wc.sample]["strandness"] == "forward" else "reverse" # TODO: add option not stranded
    conda:
        "envs/htseq.yml"
    shell:
        """
        echo -e "
        from src.htseq import *
        annotateReads(
            {input.bam},
            {input.gff},
            {output.bam},
            {output.discarded},
            "intersection-nonempty", # st-option: htseq_mode (union/intersection-strict/intersection-nonempty*) [!= multimappers]
            {params.strandness},     # st-option: strandness (yes*[=forward]/no/reverse)
            False,                   # st-option: htseq_no_ambiguous (True/False*) discard htseqs ambiguous annotations
            False                    # st-option: include_non_annotated (True/False*) include unannotated reads as '_nofeature'
        )" | python3
        """
        

###################
### ASSEMBLY QC ###
###################
rule detonate:
    input:
        R1 = "results/sortmerna/{sample}.mRNA_fwd.fastq.gz",
        R2 = "results/sortmerna/{sample}.mRNA_rev.fastq.gz",
        fa = assembly_input
    output:
        "results/detonate/{assembler}/{sample}/{sample}.genes.results",
        "results/detonate/{assembler}/{sample}/{sample}.isoforms.results",
        "results/detonate/{assembler}/{sample}/{sample}.score",
        "results/detonate/{assembler}/{sample}/{sample}.score.genes.results",
        "results/detonate/{assembler}/{sample}/{sample}.score.isoforms.results",
        directory("results/detonate/{assembler}/{sample}/{sample}.stat")
    log:
        "results/logs/detonate/{sample}.{assembler}.log"
    params:
        outdir = lambda wildcards, output: os.path.dirname(output[1]),
        read_length = config["read_length"]
    conda:
        "envs/detonate.yml"
    threads: 10
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 10
    shell:
        """
        rsem-eval-calculate-score {input.R1},{input.R2} {input.fa} \
            {params.outdir}/{wildcards.sample} {params.read_length} -p {threads} >{log} 2>&1
        """


##################
### ANNOTATION ###
##################

def genetic_code(wildcards):
    try:
        return config["transdecode"]["genetic_codes"][wildcards.sample]
    except KeyError:
        return "Universal"

rule transdecoder_longorfs:
    input:
        assembly_input
    output:
        expand("results/transdecoder/{{assembler}}/{{sample}}/{f}",
            f = ["base_freqs.dat", "longest_orfs.cds", "longest_orfs.gff3", "longest_orfs.pep"]),
    params:
        gencode=lambda wildcards: genetic_code(wildcards),
        tmpdir="$TMPDIR/{assembler}.{sample}",
        ln="$TMPDIR/{assembler}.{sample}/{sample}",
        fa=lambda wildcards, input: os.path.abspath(input[0]),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    log: "results/logs/transdecoder/{assembler}.{sample}.longorfs.log"
    conda: "envs/transdecoder.yml"
    shadow: "full"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 48
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        mkdir -p {params.tmpdir}
        ln -s {params.fa} {params.ln}
        TransDecoder.LongOrfs -G {params.gencode} -O {params.outdir} -t {params.ln}
        rm -rf {params.tmpdir}
        """

rule transdecoder_predict:
    input:
        assembly_input,
        expand("results/transdecoder/{{assembler}}/{{sample}}/{f}",
            f=["base_freqs.dat", "longest_orfs.cds", "longest_orfs.gff3", "longest_orfs.pep"])
    output:
        expand("results/transdecoder/{{assembler}}/{{sample}}/{{sample}}.transdecoder.{suffix}",
            suffix=["bed", "cds", "gff3", "pep"])
    params:
        gencode=lambda wildcards: genetic_code(wildcards),
        tmpdir="$TMPDIR/{assembler}.{sample}",
        ln="$TMPDIR/{assembler}.{sample}/{sample}",
        fa=lambda wildcards, input: os.path.abspath(input[0]),
        outdir=lambda wildcards, output: os.path.dirname(output[0])
    log: "results/logs/transdecoder/{assembler}.{sample}.predict.log"
    conda: "envs/transdecoder.yml"
    shadow: "full"
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 48
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        mkdir -p {params.tmpdir}
        ln -s {params.fa} {params.ln}
        TransDecoder.Predict -t {params.ln} -O {params.outdir} -G {params.gencode}
        mv {wildcards.sample}.transdecoder* {params.outdir}
        rm -r {params.tmpdir}
        """

rule busco_dl:
    output:
        files = expand("resources/busco/lineages/{{lineage}}/{f}",
            f = ["ancestral", "ancestral_variants", "dataset.cfg",
                 "lengths_cutoff", "links_to_ODB10.txt", "refseq_db.faa.gz",
                 "scores_cutoff"]),
        flag = touch("resources/busco/{lineage}.done")
    params:
        url = lambda wildcards: config["busco"]["lineages"][wildcards.lineage]["url"],
        outdir = lambda wildcards, output: f"{os.path.dirname(output.flag)}/lineages",
        tmpfile = lambda wildcards: f"$TMPDIR/{wildcards.lineage}.tar.gz"
    log: "results/logs/busco/{lineage}.download.log"
    shell:
        """
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        mkdir -p {params.outdir}
        curl -L -o {params.tmpfile} {params.url} > {log} 2>&1
        tar -C {params.outdir} -xf {params.tmpfile}
        rm -r {params.tmpfile}
        """

rule busco:
    input:
        busco_input(samples, config)

rule run_busco:
    input:
        fa = "results/transdecoder/{assembler}/{sample}/{sample}.transdecoder.pep",
        flag = "resources/busco/{lineage}.done",
    output:
        "results/busco/{assembler}/{sample}/short_summary.specific.{lineage}.{sample}.txt"
    log: "results/logs/busco/{sample}.{assembler}.{lineage}.log"
    conda: "envs/busco.yml"
    params:
        dlpath = lambda wildcards, input: os.path.dirname(input.flag),
        tmpdir = "$TMPDIR/busco.{assembler}.{sample}",
        outdir = lambda wildcards, output: os.path.dirname(output[0])
    resources:
        runtime = lambda wildcards, attempt: attempt ** 2 * 60 * 24
    threads: 10
    shadow: "shallow"
    shell:
        """
        exec &> {log}
        if [ -z ${{TMPDIR+x}} ]; then TMPDIR=/scratch; fi
        mkdir -p {params.tmpdir}
        busco -f --offline --download_path {params.dlpath} -l {wildcards.lineage} \
            -i {input.fa} -m proteins --out_path {params.tmpdir} \
            -o {wildcards.sample} -c {threads} > {log} 2>&1
        rsync -azv {params.tmpdir}/{wildcards.sample}/* {params.outdir}/
        rm -rf {params.tmpdir}
        """