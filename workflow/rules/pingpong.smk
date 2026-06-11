# Collapse unique sequences and count occurrences
# ------------------------------------------------------------
rule collapse_sequences:
    input:
        fastq="results/trimmed/{sample}.fastq.gz",
    output:
        fasta="results/fasta/{sample}.collapsed.fasta",
    threads: 2
    log:
        "logs/pingpong/{sample}_collapse.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "zcat {input.fastq} | fastx_collapser -o {output.fasta} {log}"


# Download Ensembl ncRNA FASTA for the configured genome and release
# ------------------------------------------------------------
rule download_ncrna_fasta:
    output:
        fa=f"resources/ncrna/{GENOME}.ncrna.fa.gz",
    params:
        url=NCRNA_URL,
    retries: 3
    log:
        f"logs/ncrna/download.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "wget -q -O {output.fa} {params.url} 2> {log}"


# Extract excluded ncRNA biotypes to use as the filter reference.
# Reads mapping to these sequences will be removed before ping-pong
# analysis; lncRNAs are intentionally excluded from this filter so
# that reads from lncRNA-hosted piRNA clusters are retained.
# ------------------------------------------------------------
rule create_rrna_fasta:
    input:
        fa=f"resources/ncrna/{GENOME}.ncrna.fa.gz",
    output:
        fa=f"resources/ncrna/{GENOME}.rrna.fa",
    log:
        f"logs/ncrna/filter_fasta.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "seqkit grep "
        "-r "
        "-n "
        '-p "gene_biotype:(Mt_rRNA|rRNA)" '
        "{input.fa} "
        "-o {output.fa} "
        "2> {log}"


# Build bowtie1 index from the excluded ncRNA sequences
# ------------------------------------------------------------
rule bowtie_index_rrna:
    input:
        fa=f"resources/ncrna/{GENOME}.rrna.fa",
    output:
        index_files=expand(f"resources/bowtie1_index/rrna_{GENOME}.{{ext}}", ext=EXT),
    params:
        index_prefix=lambda wildcards, output: output.index_files[0].replace(
            ".1.ebwt", ""
        ),
    threads: 1
    log:
        f"logs/ncrna/bowtie_index.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fa} {params.index_prefix} > {log} 2>&1"


# Remove rRNA sequences
# Aligning reads are saved as a BAM for QC/troubleshooting.
# Alignment settings based on https://github.com/bowhan/piPipes
# ------------------------------------------------------------
rule remove_rrna_reads:
    input:
        fasta="results/fasta/{sample}.collapsed.fasta",
        index_files=expand(f"resources/bowtie1_index/rrna_{GENOME}.{{ext}}", ext=EXT),
    output:
        un_fasta="results/fasta/{sample}.rrna_removed.fasta",
        al_fasta="results/rrna/{sample}.fasta",
        bam="results/rrna/{sample}.bam",
    params:
        index_prefix=lambda wildcards, input: input.index_files[0].replace(
            ".1.ebwt", ""
        ),
    threads: 4
    log:
        "logs/rrna/{sample}_filter_reads.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie "
        "-f --sam "
        "--threads {threads} "
        "-v 2 "
        "-k 1 "
        "--best "
        "--un {output.un_fasta} "
        "--al {output.al_fasta} "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view -F 4 -bS - | "
        "samtools sort -o {output.bam}"


# Plot rRNA and miRNA read counts per sample
# ------------------------------------------------------------
rule plot_rna_counts:
    input:
        rrna=expand("results/rrna/{sample}.fasta", sample=SAMPLES),
        mirna=expand("results/mirna/{sample}.fasta", sample=SAMPLES),
    output:
        pdf="results/plots/rna_counts.pdf",
    log:
        "logs/plot_rna_counts.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_rna_counts.R"


rule download_mirna_fasta:
    output:
        fasta="resources/mirhairpin.fa",
    params:
        url="https://www.mirbase.org/download/hairpin.fa",
        url_ftp="ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz",
    log:
        "logs/download_mirna.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "timeout 15 wget -nv -O {output.fasta} {params.url} > {log} 2>&1 || "
        "(echo 'HTTPS failed, trying FTP fallback...' >> {log} && "
        "timeout 15 wget -nv -O - {params.url_ftp} 2>> {log} | gunzip > {output.fasta}) || "
        "(echo 'Both downloads failed. Please download hairpin.fa manually from https://mirbase.org or ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz and place it at {output.fasta}' >> {log} && exit 1)"


rule subset_mirna_fasta:
    input:
        fasta="resources/mirhairpin.fa",
    output:
        fasta=f"resources/mirhairpin_{MIRBASE_GENOME}.fa",
    params:
        genome=MIRBASE_GENOME.replace("_", " "),
    log:
        f"logs/seqkit_{MIRBASE_GENOME}.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "seqkit grep -n -r -p '{params.genome}' {input.fasta} > {output.fasta} 2> {log}"


rule mirna_index:
    input:
        fasta=f"resources/mirhairpin_{MIRBASE_GENOME}.fa",
    output:
        index_files=expand("resources/bowtie1_index/mirna.{ext}", ext=EXT),
    log:
        "logs/bowtie/mirna.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/mirna {log}"


rule remove_mirna_reads:
    input:
        fasta="results/fasta/{sample}.rrna_removed.fasta",
        index=expand("resources/bowtie1_index/mirna.{ext}", ext=EXT),
    output:
        mirna="results/mirna/{sample}.fasta",
        mirna_bam="results/mirna/{sample}.bam",
        filtered="results/fasta/{sample}.rrna_mirna_removed.fasta",
    params:
        index_prefix=lambda wildcards, input: input.index[0].replace(".1.ebwt", ""),
        mismatch=1,
    threads: 4
    log:
        "logs/mirna_correction/{sample}_count.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        # mapping parameters based on:
        # https://www.nature.com/articles/s41467-023-42787-1
        # https://github.com/bowhan/piPipes
        # --un:  truly unmapped reads (0 alignments)
        # --max: reads mapping to >M miRNA hairpins (multi-mappers)
        # Both are merged into {output.filtered} so that piRNA-like reads
        # sharing sequence with miRNA families are not silently discarded.
        "TMP_UN=$(mktemp -p {resources.tmpdir}) && "
        "TMP_MAX=$(mktemp -p {resources.tmpdir}) && "
        "bowtie "
        "-f --sam -M 1 --best --strata --chunkmbs 1024 "
        "--threads {threads} "
        "-x {params.index_prefix} "
        "--un $TMP_UN "
        "--max $TMP_MAX "
        "--al {output.mirna} "
        "{input.fasta} 2> {log} | "
        "samtools view -F 4 -bS - | "
        "samtools sort -o {output.mirna_bam} && "
        "cat $TMP_UN $TMP_MAX > {output.filtered} && "
        "rm $TMP_UN $TMP_MAX && "
        "touch {output.mirna}"


# Build bowtie1 index for ping-pong analysis
# --------------------------------------------
rule bowtie_index:
    input:
        fasta=config["pingpong"]["fasta"],
    output:
        index_files=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT),
    params:
        index_prefix=lambda wildcards, output: output.index_files[0].replace(
            ".1.ebwt", ""
        ),
    threads: 1
    log:
        "logs/bowtie/pingpong.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/pingpong {log}"


# Align collapsed sequences to bowtie1 index
# Bowtie parameters based on:
# https://github.com/rberrens/SPOCD1-piRNA_directed_DNA_met/blob/master/piRNA_analysis/snakemake/13_map_consensus_L1s_IAPs.py
# ------------------------------------------------------------
rule align:
    input:
        fasta="results/fasta/{sample}.rrna_mirna_removed.fasta",
        index=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT),
    output:
        bam="results/pingpong/{sample}.bam",
    params:
        index_prefix=lambda wildcards, input: input.index[0].replace(".1.ebwt", ""),
        mismatch=config["pingpong"]["mismatch"],
    threads: 4
    log:
        "logs/pingpong/{sample}_align.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie "
        "-f --sam "
        "--threads {threads} "
        "-v {params.mismatch} "
        "-a --best -k 1 "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view "
        "-F 4 -bS - | "
        "samtools sort -o {output.bam}"


# Perform ping-pong analysis
# ------------------------------------------------------------
rule pingpong_analysis:
    input:
        bam="results/pingpong/{sample}.bam",
    output:
        pingpong="results/pingpong/{sample}.csv",
        nt_bias="results/pingpong/{sample}_nt_bias.csv",
    params:
        window=config["pingpong"]["window"],
    threads: 2
    log:
        "logs/pingpong/{sample}_analysis.log",
    conda:
        "../envs/pingpong.yaml"
    script:
        "../scripts/pingpong.py"


# Plot ping-pong results
# ------------------------------------------------------------
rule plot_pingpong:
    input:
        pingpong=expand("results/pingpong/{sample}.csv", sample=SAMPLES),
    output:
        pdf="results/plots/pingpong.pdf",
        csv="results/plots/pingpong.csv",
    threads: 1
    log:
        "logs/pingpong/plot_pingpong.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_pingpong.R"


# Plot 1U/10A bias
# ------------------------------------------------------------
rule plot_sequence_bias_pirna:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
        nt_bias=expand("results/pingpong/{sample}_nt_bias.csv", sample=SAMPLES),
    output:
        logo="results/plots/pirna_sequence_bias/{te}_logo.pdf",
        bar="results/plots/pirna_sequence_bias/{te}_bargraph.pdf",
        csv="results/plots/pirna_sequence_bias/{te}_frequencies.csv",
    threads: 2
    log:
        "logs/pingpong/{te}_plot_sequence_bias.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/sequence_bias_pirna.R"


# Plot per-TE positional read coverage (forward/reverse strand)
# ------------------------------------------------------------
rule te_pos_coverage:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
    output:
        pdf="results/plots/te_pos_coverage/{te}.pdf",
    threads: 2
    log:
        "logs/pingpong/{te}_te_pos_coverage.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_pos_coverage.R"


"""
# Get length distribution of TE-aligned reads
# ------------------------------------------------------------
rule length_distribution_aligned_to_TE:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
    output:
        pdf="results/plots/length_distribution.pdf",
        csv="results/plots/length_distribution.csv"
    threads: 1
    log:
        "logs/pingpong/length_distribution.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/length_distribution.R"
"""
