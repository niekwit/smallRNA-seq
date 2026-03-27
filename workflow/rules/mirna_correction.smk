rule download_mirna_fasta:
    output:
        fasta = "resources/mirhairpin.fa"
    params:
        url = "https://www.mirbase.org/download/hairpin.fa"
    log:
        "logs/download_mirna.log"
    shell:
        "wget -O {output.fasta} {params.url} > {log} 2>&1"


rule subset_mirna_fasta:
    input:
        fasta="resources/mirhairpin.fa"
    output:
        fasta=f"resources/mirhairpin_{MIRBASE_GENOME}.fa"
    params:
        genome=MIRBASE_GENOME.replace("_", " ")
    log:
        f"logs/seqkit_{MIRBASE_GENOME}.log"
    conda:
        "../envs/process_reads.yaml"
    shell:
        "seqkit grep -n -r -p '{params.genome}' {input.fasta} > {output.fasta} 2> {log}"


rule mirna_index:
    input:
        fasta=f"resources/mirhairpin_{MIRBASE_GENOME}.fa",
    output:
        index_files=expand("resources/bowtie1_index/mirna.{ext}", ext=EXT)
    log:
        "logs/bowtie/mirna.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/mirna {log}"


rule align_to_mirna:
    input:
        fasta="results/fasta/{sample}.collapsed.fasta",
        index=expand("resources/bowtie1_index/mirna.{ext}", ext=EXT)
    output:
        counts="results/mirna_correction/{sample}.txt"
    params:
        index_prefix="resources/bowtie1_index/mirna",
        mismatch=0
    threads: 4
    log:
        "logs/mirna_correction/{sample}_count.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        # mapping parameters according to:
        # https://www.nature.com/articles/s41467-023-42787-1
        "bowtie "
        "-n 2 -M 1 --best --strata --chunkmbs 1024 "
        "-f --sam "
        "--threads {threads} "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view "
        "-F 4 - | "
        "cut -f1 | cut -d'-' -f2 | "
        "awk '{{sum+=$1}} END {{print sum}}' > {output.counts}"