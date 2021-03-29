#!usr/bin/env bash -l
set -e

### Wrapper script for running snakemake with selected options, requires conda installed ###
# Usage: bash path/to/doSnakemake.sh [additional snakemake options] [optional snakemake output file]

# get the workflow base directory (=where this script reside)
DIR=${0%runSnakemake.sh}

# Default conda is currently a bit unstable, so we want to use the alternative
# mamba implementation instead
if ! [ -x "$(command -v mamba)" ]; then
    echo -e "NB! The default conda is currently a bit unstable, so we want to" \
         "use the alternative mamba implementation instead"
    echo -e "Please install mamba as follows:\n" \
         "   'conda install mamba'\n" \
	 "and rerun this script"
    exit -1
fi

# Set up a conda env for snakemake + mamba -- so snakemake also uses mamba
snakemake="masc-seq-asm"
if [[ "$(mamba info -e | awk -v var=$snakemake '{if($1==var)print "found"}')" != "found" ]]; then
    mamba env create --file $DIR/environment.yml
fi
snakemake=$(mamba info -e | awk -v var=$snakemake '{if($1==var)print $2}')
eval "$(conda shell.bash hook)"
conda activate $snakemake

# # Allow for either running on cluster or on local computer/laptop
if [ "$CLUSTER" = "rackham" ]; then # Change/add cluster name if needed
    echo "Running on rackham"
    time snakemake \
	 --conda-frontend mamba \
	 -j \
	 -k \
	 --snakefile $DIR/Snakefile \
	 --use-conda \
	 --profile s$DIR/slurm \
	 $@
	 # --use-singularity \
	 # --cluster-config $DIR/config/cluster.yml \
	 # --cluster " sbatch -J {cluster.name} -A {cluster.account} \
	 #             -p {cluster.partition} -n {cluster.n} -t {cluster.time} \
         #             {cluster.other} " \
    # Explanation of command and options:
    # `time`           just gives execution time report
    # --conda-frontend mamba 
    #                  tells snakemake to use the improved version of conda called `mamba`
    # -j               Tells snakemake all cores available for the job
    # -s               tells which Snakefile to use
    # --use-conda      tells snakemake to create and use rule-specific conda envs
    # --use-singularty together with a `container: <mycontainer>` in Snakefile,
    #                  this tells snakemake to run inside the container <mycontainer>
    # --cluster-config tells what cluster-config file to use
    # --cluster " sbatch -J {cluster.name} ..."
    #                  This is the sbatch command with variables
    #                  taken from the cluster-config
    # $@               pass additional arguments (requested outout file, options)
    #                  to snakemake, e.g. `doSlamdunk.sh -n myref.fasta` will add
    #                  `-n myref.fasta` directly to the end of the snakemake call
    # \                `line break`, marks that the command continues on the
    #                  next line (NB! no space or text after the `\`)
else
    # When run locally, we don't need --cluster-config or --cluster
    time snakemake \
	 --conda-frontend mamba \
	 -j \
	 --snakefile $DIR/Snakefile \
	 --use-conda \
	 $@
fi


# How snakemake handles the cluster calls:
# 1. Snakemake should be run on the login node in a virtual terminal (use tmux -- recommended -- or screen)
# 2. Snakemake first parses the Snakefile and create the workflow as a directed acyclic graph (DAG)
# 3. If needed it creates the requested conda environments
# 4. It parses the --cluster-config file and reads a python dictionary indexed by rule names and __default__
# 5. For each rule to be run, it the --cluster sbatch command and substitutes the {cluster.xxx}
#    variables for corresponding values from the --cluster-config for the rule
# 6. It sends the job to the slurm queue -- when possible it will run rules/jobs in parallel
# 7. When job is started, it first creates any conda-environment that is needed, and then runs the rule.
# 8. When the job returns, snakemake checks it it has failed or succeeded and then decide this
#    means that it should start another job/rule that depends on the output of the one just run.
