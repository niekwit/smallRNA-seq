# Install TEsmall software from GitHub
# --------------------------------------------------------
rule install_te_small:
    output:
        directory("resources/te_small/"),
    params:
        url="https://github.com/mhammell-laboratory/TEsmall.git",
        version="-b 2.0.9",
    retries: 3
    log:
        "logs/install_te_small.log",
    conda:
        "../envs/te_small.yaml"
    shell:
        "git clone "
        "{params.url} "
        "{params.version} "
        "{output}; "
        "cd {output}; "
        "pip install --no-deps --prefix=. ."


# Run TEsmall analysis
# --------------------------------------------------------
rule run_te_small:
    input:
        unpack(get_dependencies),
    output:
        txt="results/te_small/{sample}/count_summary.txt",
        report="results/te_small/{sample}/report.html",
        fastq="results/te_small/{sample}/{sample}.rm_rRNA.fastq",
        rlen="results/te_small/{sample}/{sample}.anno.rlen.info",
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


# Differential expression analysis with DESeq2
# --------------------------------------------------------
rule deseq2:
    input:
        txt=expand("results/te_small/{sample}/count_summary.txt", sample=SAMPLES),
    output:
        pdf="results/plots/{comparison}/volcano_plot.pdf",
        csv="results/plots/{comparison}/results.csv",
        rds="results/deseq2/{comparison}/dds.rds",
    log:
        "logs/te_small/deseq2_{comparison}.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/deseq2.R"


# Plot PCA of DESeq2 results
# --------------------------------------------------------
rule plot_pca:
    input:
        rds=expand("results/deseq2/{comparison}/dds.rds", comparison=COMPARISONS[0]),
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
