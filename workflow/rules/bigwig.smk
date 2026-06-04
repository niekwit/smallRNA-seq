# Sort TEsmall BAM files by coordinate
# ------------------------------------------------------------
rule sort_te_small_bam:
    input:
        bam="results/te_small/{sample}/{sample}.3trf_free.bam",
    output:
        sorted_bam="results/te_small/{sample}/{sample}.3trf_free.sorted.bam",
    log:
        "logs/samtools_sort/{sample}.log",
    threads: 2
    conda:
        "../envs/pingpong.yaml"
    shell:
        "samtools sort {input.bam} -o {output.sorted_bam} 2> {log}"


# Index TEsmall BAM files
# ------------------------------------------------------------
rule index_te_small_bam:
    input:
        bam="results/te_small/{sample}/{sample}.3trf_free.sorted.bam",
    output:
        bai="results/te_small/{sample}/{sample}.3trf_free.sorted.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    threads: 1
    conda:
        "../envs/pingpong.yaml"
    shell:
        "samtools index {input.bam} 2> {log}"


# Intersect TEsmall genome-mapped reads with miRNA annotations to get miRNA-mapped reads
# ------------------------------------------------------------
rule mirna_te_small_bam:
    input:
        bam="results/te_small/{sample}/{sample}.3trf_free.sorted.bam",
        bai="results/te_small/{sample}/{sample}.3trf_free.sorted.bam.bai",
    output:
        bam="results/mirna/mapped_te_small_{sample}.bam",
    params:
        mirna_bed=f"{config["tesmall"]["dbfolder"]}/genomes/{GENOME}/annotation/miRNA.bed",
    log:
        "logs/mirna_te_small_bam/{sample}.log",
    threads: 1
    conda:
        "../envs/te_small.yaml"
    shell:
        "bedtools intersect -wa -u "
        "-a {input.bam} "
        "-b {params.mirna_bed} > {output.bam} 2> {log}"


# Index miRNA BAM files required by bamCoverage
# ------------------------------------------------------------
rule index_mirna_bam:
    input:
        bam="results/mirna/mapped_te_small_{sample}.bam",
    output:
        bai="results/mirna/mapped_te_small_{sample}.bam.bai",
    threads: 1
    log:
        "logs/samtools_index/{sample}.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "samtools index {input.bam} 2> {log}"


# Effective genome sizes from:
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
EFFECTIVE_GENOME_SIZES = {
    ("hg38", 50): 2_684_781_695,
    ("hg38", 75): 2_736_124_973,
    ("hg38", 100): 2_776_919_808,
    ("hg38", 150): 2_827_437_033,
    ("hg38", 200): 2_855_464_000,
    ("mm39", 50): 2_308_125_440,
    ("mm39", 75): 2_407_883_318,
    ("mm39", 100): 2_467_481_008,
    ("mm39", 150): 2_549_647_237,
    ("mm39", 200): 2_576_780_534,
}

# Map each condition to its sample names for use in expand() calls
condition_samples = {}
for sample, condition in config["samples"].items():
    condition_samples.setdefault(condition, []).append(sample)


# Convert miRNA BAM to bigWig using deeptools bamCoverage
# ------------------------------------------------------------
rule mirna_bigwig:
    input:
        bam="results/mirna/mapped_te_small_{sample}.bam",
        bai="results/mirna/mapped_te_small_{sample}.bam.bai",
    output:
        "results/bigwig/mirna_{sample}.bw",
    params:
        effective_genome_size=lambda wildcards: EFFECTIVE_GENOME_SIZES[
            (GENOME, config["bigwig"]["read_length"])
        ],
        normalize_using=config["bigwig"]["normalize_using"],
        extra=config["bigwig"].get("extra", ""),
    threads: 4
    log:
        "logs/bigwig/mirna_{sample}.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--effectiveGenomeSize {params.effective_genome_size} "
        "--normalizeUsing {params.normalize_using} "
        "-p {threads} "
        "{params.extra} "
        "2> {log}"


# Average per-sample bigWigs into a single bigWig per condition
# ------------------------------------------------------------
rule mirna_bigwig_average:
    input:
        lambda wildcards: expand(
            "results/bigwig/mirna_{sample}.bw",
            sample=condition_samples[wildcards.condition],
        ),
    output:
        "results/bigwig/mirna_{condition}_average.bw",
    threads: 12
    log:
        "logs/bigwig/mirna_{condition}_average.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bigwigAverage "
        "-b {input} "
        "-o {output} "
        "-p {threads} "
        "2> {log}"
