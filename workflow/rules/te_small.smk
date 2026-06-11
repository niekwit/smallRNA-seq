# Run TEsmall analysis
# --------------------------------------------------------
rule te_small:
    input:
        unpack(get_dependencies),
    output:
        txt="results/te_small/{sample}/count_summary.txt",
        report="results/te_small/{sample}/report.html",
        fastq="results/te_small/{sample}/{sample}.rm_rRNA.fastq",
        rlen="results/te_small/{sample}/{sample}.anno.rlen.info",
        bam="results/te_small/{sample}/{sample}.3trf_free.bam",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.txt),
        genome=GENOME,
        adapter=config["cutadapt"]["adapter"],
        dbfolder=config["tesmall"]["dbfolder"],
        minlen=config["tesmall"]["minlen"],
        maxlen=config["tesmall"]["maxlen"],
        maxaln=config["tesmall"]["maxaln"],
        mismatch=config["tesmall"]["mismatch"],
    threads: 4
    log:
        "logs/te_small/{sample}.log",
    conda:
        "../envs/te_small.yaml"
    script:
        "../scripts/te_small.py"


# Build DESeq2 object from all samples/conditions (used by PCA and per-comparison analysis)
# --------------------------------------------------------
rule deseq2_full:
    input:
        txt=expand("results/te_small/{sample}/count_summary.txt", sample=SAMPLES),
    output:
        rds="results/deseq2/dds_full.rds",
    log:
        "logs/te_small/deseq2_full.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/deseq2_full.R"


# Per-comparison differential expression analysis with DESeq2
# --------------------------------------------------------
rule deseq2:
    input:
        dds_full="results/deseq2/dds_full.rds",
    output:
        pdf="results/plots/{comparison}/volcano_plot.pdf",
        csv="results/plots/{comparison}/results.csv",
        dds="results/deseq2/{comparison}/dds.rds",
    log:
        "logs/te_small/deseq2_{comparison}.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/deseq2.R"


# Plot PCA of DESeq2 results (all conditions)
# --------------------------------------------------------
rule plot_pca:
    input:
        rds="results/deseq2/dds_full.rds",
    output:
        pdf="results/plots/pca_plot.pdf",
    log:
        "logs/te_small/plot_pca.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_pca.R"


# Pie chart of small RNA feature composition per condition
# --------------------------------------------------------
rule te_small_piechart:
    input:
        rlen=expand("results/te_small/{sample}/{sample}.anno.rlen.info", sample=SAMPLES),
    output:
        pdf="results/plots/te_small/piechart.pdf",
        csv="results/plots/te_small/piechart.csv",
    log:
        "logs/te_small/piechart.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_small_piechart.R"


# Bar graph of small RNA feature composition per read length per condition
# --------------------------------------------------------
rule te_small_bargraph:
    input:
        rlen=expand("results/te_small/{sample}/{sample}.anno.rlen.info", sample=SAMPLES),
    output:
        pdf="results/plots/te_small/bargraph.pdf",
    params:
        minlen=config["tesmall"]["minlen"],
        maxlen=config["tesmall"]["maxlen"],
    log:
        "logs/te_small/bargraph.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_small_bargraph.R"


# Bar graph of top TE subfamilies per class (All vs antisense piRNAs)
# --------------------------------------------------------
rule te_small_top_families:
    input:
        txt=expand("results/te_small/{sample}/count_summary.txt", sample=SAMPLES),
    output:
        pdf="results/plots/te_small/te_top_families.pdf",
    params:
        n_top=10,  # top N subfamilies per TE class
        te_classes=config["tesmall"]["te_classes"],
    log:
        "logs/te_small/te_top_families.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_small_top_families.R"


# Stacked bar of TE class composition of piRNAs per condition
# --------------------------------------------------------
rule te_small_stacked_bar_te:
    input:
        txt=expand("results/te_small/{sample}/count_summary.txt", sample=SAMPLES),
    output:
        pdf="results/plots/te_small/stacked_bar_te.pdf",
        csv="results/plots/te_small/stacked_bar_te.csv",
    log:
        "logs/te_small/stacked_bar_te.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_small_stacked_bar_te.R"
