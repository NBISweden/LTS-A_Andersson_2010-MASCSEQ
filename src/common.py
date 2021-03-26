#!/usr/bin/env python

import pandas as pd

def read_samples(f):
    df = pd.read_csv(f, index_col=0)
    return df.to_dict(orient="index")


def dammit_input(wildcards):
    if wildcards.assembler == "transabyss":
        i = "merged.fa"
    elif wildcards.assembler == "trinity":
        i = "Trinity.fasta"
    return f"results/{wildcards.assembler}/{wildcards.sample}/{i}"

""" It is good keep the code and the results separate, so
    we run the workflow in a separate working directory. However,
    we sometimes need to access files in the code folder. This
    function prepends a file path with the path to the workflow
    base directory (kept in snakemake variable workflow.basedir) """
def prependWfd(path):
    return(os.path.normpath(os.path.join(workflow.basedir, path)))

""" Mostly, snakemake uses relative paths (to the working directory)
    when referring to, e.g., input and output files. However, when
    soft-linking (creating aliases) from external file, full paths
    must be used. This function prepends the path to the current
    working directory (accessed by the python function os.getcwd()). """
def prependPwd(path):
    return(os.path.normpath(os.path.join(os.getcwd(), path)))
