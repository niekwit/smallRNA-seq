# Plot length distribution of trimmed sequences
# -----------------------------------------------------
rule plot_length_distribution:
    input:
        fasta=expand("results/fasta/{sample}_counts.fasta", sample=SAMPLES),
    output:
        pdf="results/plots/length_distribution.pdf",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_length_distribution.R"