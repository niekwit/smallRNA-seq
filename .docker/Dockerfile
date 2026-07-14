FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7bfe7e929a740a55b3e850ec175cdb69f83fecbd71cbe4927498cbe92d75a309"

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
#   prefix: /conda-envs/7a22c6098e6e8d854820923f829b8498
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
#     - r-ggrepel=0.9.8
#     - bioconda::bioconductor-deseq2=1.46.0
#     - bioconda::bioconductor-biostrings=2.74.0
#     - bioconda::bioconductor-rsamtools=2.22.0
#     - conda-forge::r-gridextra=2.3
#     - conda-forge::jq
RUN mkdir -p /conda-envs/7a22c6098e6e8d854820923f829b8498
COPY workflow/envs/R.yaml /conda-envs/7a22c6098e6e8d854820923f829b8498/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/aaa4130c6d1032f188e6d44f6e17ca0d
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - deeptools=3.5.6
#     - pygenometracks=3.9
#     - pybigwig=0.3.25
#     - biopython=1.87
#     - pandas
RUN mkdir -p /conda-envs/aaa4130c6d1032f188e6d44f6e17ca0d
COPY workflow/envs/deeptools.yaml /conda-envs/aaa4130c6d1032f188e6d44f6e17ca0d/environment.yaml

# Conda environment:
#   source: workflow/envs/pingpong.yaml
#   prefix: /conda-envs/03d5424c620255c2c4b5db4e7d3a747f
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
#     - bioconda::seqkit=2.13.0
#     - conda-forge::wget=1.25.0
RUN mkdir -p /conda-envs/03d5424c620255c2c4b5db4e7d3a747f
COPY workflow/envs/pingpong.yaml /conda-envs/03d5424c620255c2c4b5db4e7d3a747f/environment.yaml

# Conda environment:
#   source: workflow/envs/te_small.yaml
#   prefix: /conda-envs/2b3c7df2b713c88df1c4a2c9ccbd1a0f
#   name: te_small
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bioconda::tesmall=2.1.0
RUN mkdir -p /conda-envs/2b3c7df2b713c88df1c4a2c9ccbd1a0f
COPY workflow/envs/te_small.yaml /conda-envs/2b3c7df2b713c88df1c4a2c9ccbd1a0f/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32 --file /conda-envs/24864fe2fb080d68e0e41c2f30a1cb32/environment.yaml && \
    conda env create --prefix /conda-envs/7a22c6098e6e8d854820923f829b8498 --file /conda-envs/7a22c6098e6e8d854820923f829b8498/environment.yaml && \
    conda env create --prefix /conda-envs/aaa4130c6d1032f188e6d44f6e17ca0d --file /conda-envs/aaa4130c6d1032f188e6d44f6e17ca0d/environment.yaml && \
    conda env create --prefix /conda-envs/03d5424c620255c2c4b5db4e7d3a747f --file /conda-envs/03d5424c620255c2c4b5db4e7d3a747f/environment.yaml && \
    conda env create --prefix /conda-envs/2b3c7df2b713c88df1c4a2c9ccbd1a0f --file /conda-envs/2b3c7df2b713c88df1c4a2c9ccbd1a0f/environment.yaml && \
    conda clean --all -y