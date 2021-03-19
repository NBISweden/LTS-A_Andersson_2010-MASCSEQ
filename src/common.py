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
