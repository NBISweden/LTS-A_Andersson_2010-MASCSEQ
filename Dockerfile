FROM continuumio/miniconda3:4.9.2

MAINTAINER "John Sundh" john.sundh@nbis.se

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /project

# Set tmpdir
ENV TMPDIR="/scratch"

# Install necessary tools
RUN apt-get update && apt-get install -y --no-install-recommends bzip2 curl unzip wget && apt-get clean

# Add Conda environment file
COPY environment.yml ./

# Install the Conda environment
RUN conda install -c conda-forge mamba \
    && mamba env update -n base -f environment.yml \
    && mamba clean --all

# Start Bash shell by default
CMD /bin/bash