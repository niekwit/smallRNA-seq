## Workflow overview

This workflow analyses small RNA-seq data, with a focus on TE-derived piRNAs.
It is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Trim 3' adapters (`cutadapt`) and plot read length distributions
2. Annotate small RNAs and count TE/miRNA reads (`TEsmall`)
3. Differential expression analysis, PCA and volcano plots (`DESeq2`)
4. Collapse reads (`fastx_collapser`) and remove rRNA (Ensembl ncRNA) and miRNA (miRBase hairpin) reads (`bowtie`)
5. Align the remaining reads to the TE consensus sequences (`bowtie`)
6. Ping-pong analysis (5' overlap distances between sense/antisense reads) and piRNA 1U/10A sequence bias
7. Strand-split TE coverage bigWigs (`deepTools`) and TE track plots (`pyGenomeTracks`)

See the main [README](../README.md) for the expected directory layout, the full list of output files and a detailed explanation of the bigWig normalisation options.

## Running the workflow

### Input data

Gzipped single-end FASTQ files, one per sample, at `reads/{sample}/{sample}.fastq.gz`.
Sample names and their conditions are set in `config/config.yaml` under `samples`.

A FASTA file with the TE consensus sequences of interest (`pingpong.fasta`) is needed for the ping-pong analysis and the TE coverage tracks.

### Parameters

This table lists all parameters that can be used to run the workflow.

| parameter                   | type        | details                                                                                                      | example                  |
| --------------------------- | ----------- | ------------------------------------------------------------------------------------------------------------ | ------------------------ |
| genome                      | str         | reference genome, `hg38` or `mm39`                                                                           | `mm39`                   |
| release                     | int         | Ensembl release used for reference data downloads                                                            | `115`                    |
| samples                     | map         | `sample: condition` pairs                                                                                    | `KO_1: KO`               |
| reference_condition         | str         | baseline condition for differential analysis and plot ordering                                               | `Control`                |
| **cutadapt**                |             |                                                                                                              |                          |
| adapter                     | str         | 3' adapter sequence to trim                                                                                  | `AGATCGGAAGAGCACACGTCT`  |
| extra                       | str         | additional `cutadapt` arguments                                                                              | `"-m 16 --trimmed-only -q 20"` |
| **tesmall**                 |             |                                                                                                              |                          |
| dbfolder                    | str         | path to the TEsmall database directory (built on the first run if absent)                                    | `/path/to/te_small_db`   |
| minlen                      | int         | minimum read length                                                                                          | `16`                     |
| maxlen                      | int         | maximum read length                                                                                          | `36`                     |
| maxaln                      | int         | maximum number of alignments per read                                                                        | `100`                    |
| mismatch                    | int         | number of mismatches allowed                                                                                 | `0`                      |
| extra_args                  | str         | additional TEsmall arguments                                                                                 | `""`                     |
| te_classes                  | str or list | TE classes for the top-families plot: `All` or a list (DNA, LINE, LTR, RC, Retroposon, Satellite, SINE, Unknown, Unspecified) | `[LINE, LTR]` |
| **length_distribution**     |             |                                                                                                              |                          |
| apply_mirna_correction      | bool        | normalise to miRNA reads; `False` normalises to total counts                                                 | `False`                  |
| mirna_correction_method     | str         | `length` (reads within the miRNA length window) or `align` (reads aligned to miRNA hairpins)                 | `length`                 |
| mirna_min                   | int         | minimum length of the miRNA window (`length` method only)                                                    | `19`                     |
| mirna_max                   | int         | maximum length of the miRNA window (`length` method only)                                                    | `25`                     |
| **bigwig**                  |             |                                                                                                              |                          |
| bin_size                    | int         | `bamCoverage` bin size (bp)                                                                                  | `10`                     |
| read_length                 | int         | read length used to select the effective genome size (50, 75, 100, 150 or 200)                               | `50`                     |
| normalize_using             | str         | `RPKM`, `CPM`, `BPM`, `RPGC`, `None`, `totalReads` or `miRNAReads`                                           | `CPM`                    |
| extra                       | str         | additional `bamCoverage` arguments                                                                           | `""`                     |
| **pingpong**                |             |                                                                                                              |                          |
| bowtie_params               | str         | `bowtie` alignment arguments                                                                                 | `"--best -k 1 -v 3"`     |
| fasta                       | str         | FASTA file with TE consensus sequences                                                                       | `resources/TE.fa`        |
| min_len                     | int         | minimum read length (inclusive) included in the ping-pong analysis and TE coverage tracks                    | `23`                     |
| max_len                     | int         | maximum read length (inclusive) included in the ping-pong analysis and TE coverage tracks                    | `32`                     |
| window                      | int         | maximum 5' overlap distance (nt) recorded between sense/antisense read pairs                                 | `30`                     |
| plot.type                   | str         | `line` or `bar`                                                                                              | `line`                   |
| plot.bar_position           | str         | `grouped` (dodged bars) or `overlap` (overlaid, semi-transparent bars); `bar` only                           | `grouped`                |
| **ngs_tracker** (optional)  |             |                                                                                                              |                          |
| enabled                     | bool        | register the run with [NGS Tracker](https://github.com/niekwit/ngs-tracker) on success/failure; see `config.yaml` for the other fields | `false` |
