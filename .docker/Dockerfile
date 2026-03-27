FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bda57b6807c8922f85694577648e31c83b0ae0ace36d94f7f74c83c56e67c514"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.9.0/bio/cutadapt/se/environment.yaml
#   prefix: /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =5.2
RUN mkdir -p /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.9.0/bio/cutadapt/se/environment.yaml /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/971e30874532c7c91e6494d19e958a88
#   name: R
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - r-base=4.4
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.2.0
#     - r-ggseqlogo=0.2.2
#     - bioconda::bioconductor-deseq2=1.46.0
#     - bioconda::bioconductor-biostrings=2.74.0
#     - bioconda::bioconductor-rsamtools=2.22.0
#     - conda-forge::r-gridextra=2.3
RUN mkdir -p /conda-envs/971e30874532c7c91e6494d19e958a88
COPY workflow/envs/R.yaml /conda-envs/971e30874532c7c91e6494d19e958a88/environment.yaml

# Conda environment:
#   source: workflow/envs/pingpong.yaml
#   prefix: /conda-envs/c3303e4aa2c782fab6f9a32338f36644
#   name: pingpong
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bioconda::bowtie=1.3.1
#     - bioconda::fastx_toolkit=0.0.14
#     - conda-forge::python=3.12.12
#     - bioconda::pysam=0.23.3
#     - bioconda::samtools=1.22.1
#     - bioconda::seqtk=1.5
RUN mkdir -p /conda-envs/c3303e4aa2c782fab6f9a32338f36644
COPY workflow/envs/pingpong.yaml /conda-envs/c3303e4aa2c782fab6f9a32338f36644/environment.yaml

# Conda environment:
#   source: workflow/envs/process_reads.yaml
#   prefix: /conda-envs/01ea981e63e66e6a34e85f55e825d6be
#   name: process_reads
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bioconda::cutadapt=5.2
#     - bioconda::seqkit=2.12.0
#     - bioconda::bowtie2=2.5.4
RUN mkdir -p /conda-envs/01ea981e63e66e6a34e85f55e825d6be
COPY workflow/envs/process_reads.yaml /conda-envs/01ea981e63e66e6a34e85f55e825d6be/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32 --file /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml && \
    conda env create --prefix /conda-envs/971e30874532c7c91e6494d19e958a88 --file /conda-envs/971e30874532c7c91e6494d19e958a88/environment.yaml && \
    conda env create --prefix /conda-envs/c3303e4aa2c782fab6f9a32338f36644 --file /conda-envs/c3303e4aa2c782fab6f9a32338f36644/environment.yaml && \
    conda env create --prefix /conda-envs/01ea981e63e66e6a34e85f55e825d6be --file /conda-envs/01ea981e63e66e6a34e85f55e825d6be/environment.yaml && \
    conda clean --all -y