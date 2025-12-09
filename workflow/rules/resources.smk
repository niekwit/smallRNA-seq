# Download genome index and annotation files from Hammell lab Dropbox
# -----------------------------------------
rule get_genome:
    output:
        resources=temp("resources/{genome}.tar.gz")
        genome_dir=directory("resources/{genome}")
    params:
        url=genome_url(GENOME)
    threads: 1
    resources:
        mem_mb=2000,
        time=30,
    shell:
        """
        wget -O {output.resources} "{params.url}"
        tar -xzf {output.resources} -C resources/
        """


# Build bowtie1 index for ping-pong analysis
# --------------------------------------------
rule build_bowtie1_index:
    input:
        fasta=config["pingpong"]["fasta"]
    output:
        index_files=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT)
    params:
        index_prefix=lambda wildcards, output: output.index_files[0].replace(".1.ebwt", "")
    threads: 1
    log:
        "logs/pingpong/build_bowtie1_index.log"
    conda:
        "../envs/bowtie1.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/pingpong {log}"
    