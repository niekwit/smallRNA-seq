rule install_te_small:
    output:
        directory("resources/te_small/"),
    params:
        url="https://github.com/mhammell-laboratory/TEsmall.git",
        version="-b 2.0.9"
    retries: 3
    log:
        "logs/install_te_small.log"
    conda:
        "../envs/te_small.yaml"
    shell:
        "git clone "
        "{params.url} "
        "{params.version} "
        "{output}; "
        "cd {output}; "
        "python setup.py install"


rule run_te_small:
    input:
        unpack(get_dependencies)
    output:
        txt="results/te_small/{sample}/count_summary.txt",
        report="results/te_small/{sample}/report.html",
        fastq="results/te_small/{sample}/{sample}.rm_rRNA.fastq",
    params:
        outdir= lambda wildcards, output: os.path.dirname(output.txt),
        genome=GENOME,
        adapter=config["cutadapt"]["adapter"],
        dbfolder=config["tesmall"]["dbfolder"],
        minlen=config["tesmall"]["minlen"],
        maxlen=config["tesmall"]["maxlen"],
        maxaln=config["tesmall"]["maxaln"],
        mismatch=config["tesmall"]["mismatch"],
    threads: 8
    log:
        "logs/te_small/{sample}.log"
    conda:
        "../envs/te_small.yaml"
    script:
        "../scripts/te_small.py"


rule deseq2:
    input:
        txt=expand("results/te_small/{sample}/count_summary.txt", sample=SAMPLES),
    output:
        plot="results/plots/{comparison}/volcano_plot.pdf",
        csv="results/plots/{comparison}/results.csv",
        pca="results/plots/{comparison}/pca_plot.pdf",
    log:
        "logs/te_small/deseq2.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/deseq2.R"