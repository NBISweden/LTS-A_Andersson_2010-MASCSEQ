FROM continuumio/miniconda3:4.9.2

MAINTAINER "John Sundh" john.sundh@nbis.se

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /project

# Set tmpdir
ENV TMPDIR="/scratch"

RUN mkdir $TMPDIR

# Install necessary tools
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 curl unzip wget && apt-get clean

RUN mkdir envs config src
# Add Conda environment file
COPY environment.yml .
# Add snakefile
COPY Snakefile .
# Add slurm profile
COPY slurm .
# Add src
COPY src/common.py ./src/common.py
# Add environments
COPY envs/*.yml ./envs/
# Add config
COPY config/*.yml ./config/

# Install the Conda environment
RUN conda install -c conda-forge mamba \
    && mamba env update -n base -f environment.yml \
    && mamba clean --all

ENTRYPOINT ["snakemake", "-rpk", "--use-conda", "--conda-frontend", "mamba"]
CMD ["-j", "4"]