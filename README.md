# A_Andersson_2010-MASCSEQ

## Introduction
This repository contains code for primarily RNA-seq analysis for the NBIS
long term support project A_Andersson_2010.

Two RNA-seq libraries will be analysed:

| NGI ID | User ID | Organism | Mreads | >=Q30(%) |
| ------ | ------- | -------- | ------ | -------- |
| P18363_101 | Ph6B | _Phaeodactylum tricornutum_ (diatom) | 207.86  | 95.02 |
| P18363_103 | He2B | _Heterocapsa sp._ (dinoflagellate) | 158.81 | 95.94 |

## Setup

Clone this repository to a convenient place on your computer:

```
git clone https://github.com/NBISweden/LTS-A_Andersson_2010-MASCSEQ.git
```

This creates a folder with a copy of all files in the repository; we will refer to this folder as the `gwd` (github working directory).

## Run

We suggest running the analyses in a dedicated folder, here called the `awd` (analysis working directory), that is different from the `gwd`.

For convenience, use the wrapper script `runSnakemake.sh` that encapsulates the snakemake command with necessary options:

```
bash path/to/runSnakemake.sh [additional options] [optional requested output file]
```

This sets up a conda environment for snakemake (if needed) and starts the snakemake forkflow. If no `optional requested output file]` is given, default outfiles will be produced (given by the rule `all` in Snakefile). It is possible to pass `[additional options]` to snakemake. One convenient option is the dry-run option `-n`; this will cause snakemake to show what will be done by the requested analysis, but it will not actually perform it -- very important safety check before running:)
