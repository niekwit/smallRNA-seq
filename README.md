# Snakemake workflow: `smallRNA-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.25.5-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/smallRNA-seq/actions/workflows/main.yaml/badge.svg)](https://github.com/niekwit/smallRNA-seq/actions/workflows/main.yaml)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/niekwit/smallRNA-seq)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19258358.svg)](https://doi.org/10.5281/zenodo.19258358)

A Snakemake workflow for `smallRNA-seq` analysis.

Based on [Zoch et al.](https://www.nature.com/articles/s41586-020-2557-5), [Han et al.](https://academic.oup.com/bioinformatics/article/31/4/593/2748154?login=true), [Castañeda et al.](https://link.springer.com/article/10.15252/embj.201386855), and [Van Lopik et al.](https://www.nature.com/articles/s41467-023-42787-1).

## Usage

Prepare your analysis directory as follows:

```shell
.
├── config
│   ├── config.yaml
│   ├── README.md
│   └── schemas
│       └── config.schema.yaml
├── reads
│   ├── KO_1
│   │   └── KO_1.fastq.gz
│   ├── KO_2
│   │   └── KO_2.fastq.gz
│   ├── WT_1
│   │   └── WT_1.fastq.gz
│   └── WT_2
│       └── WT_2.fastq.gz
├── resources
│   └── TE.fa
└── workflow
    ├── envs
    │   ├── deeptools.yaml
    │   ├── pingpong.yaml
    │   ├── R.yaml
    │   └── te_small.yaml
    ├── rules
    │   ├── bigwig.smk
    │   ├── common.smk
    │   ├── length_distribution.smk
    │   ├── pingpong.smk
    │   └── te_small.smk
    ├── scripts
    │   ├── deseq2_full.R
    │   ├── deseq2.R
    │   ├── expand_collapsed_bam.py
    │   ├── pingpong.py
    │   ├── plot_length_distribution.R
    │   ├── plot_length_distribution_tesmall.R
    │   ├── plot_pca.R
    │   ├── plot_pingpong.R
    │   ├── plot_rna_counts.R
    │   ├── pygenometracks.py
    │   ├── sequence_bias_pirna.R
    │   ├── stranded_te_bigwig.py
    │   ├── strand_summary_table.py
    │   ├── te_pos_coverage.R
    │   ├── te_small_bargraph.R
    │   ├── te_small_cpm_te_classes.R
    │   ├── te_small.py
    │   ├── te_small_piechart.R
    │   ├── te_small_stacked_bar_te.R
    │   └── te_small_top_families.R
    └── Snakefile

12 directories, 38 files
```

The workflow and config directories can be copied from this repository.

The analysis settings are defined in `config/config.yaml`:

```yaml
genome: mm39 # hg38 or mm39
release: 115

samples:
  # sample: condition
  control_1: Control
  control_2: Control
  KO_1: KO
  KO_2: KO
reference_condition: Control

cutadapt:
  adapter: AGATCGGAAGAGCACACGTCT
  extra: "-m 16 --trimmed-only -q 20"

tesmall:
  run_tesmall: False
  dbfolder: /path/to/te_small_db
  minlen: 16
  maxlen: 36
  maxaln: 100
  mismatch: 0
  extra_args: ""
  # TE classes to include in top-families plot
  # Use "All" or a list: [DNA, LINE, LTR, RC, Retroposon, Satellite, SINE, Unknown, Unspecified]
  te_classes: All

length_distribution:
  apply_mirna_correction: TRUE
  mirna_correction_method: length # length or align

  # miRNA length range for correction
  # (for mirna_correction_method = length)
  mirna_min: 19
  mirna_max: 25

# Bigwig parameters
# See "Bigwig normalisation" below for a full explanation of normalize_using.
bigwig:
  bin_size: 10 # bamCoverage bin size (bp)
  read_length: 50 # used to select effective genome size for RPGC normalisation
  normalize_using: CPM # RPKM, CPM, BPM, RPGC, None, totalReads, miRNAReads
  extra: "" # additional bamCoverage arguments

pingpong:
  # bowtie alignment parameters for ping-pong analysis
  # Suggested parameters:
  # --best -k 1 -v 3: Maximizing alignment yield (forces a reported home for every alignable read).
  # -n 2 -M 1 --best --strata --nomaqround: High-stringency filtering (strictly grouping by mismatch count and catching multi-mappers).
  bowtie_params: "--best -k 1 -v 3"
  fasta: resources/TE.fa
  min_len: 23 # reads outside [min_len, max_len] are excluded from the ping-pong analysis
  max_len: 32
  window: 30 # maximum 5' overlap distance (nt) recorded between sense/antisense read pairs
  plot:
    type: line # line or bar
    bar_position: grouped # grouped or overlap (bar only)
```

`config.yaml` also has an optional `ngs_tracker` section for registering completed runs (and attaching key output files) with [NGS Tracker](https://github.com/niekwit/ngs-tracker). It's disabled by default (`ngs_tracker.enabled: false`); see the comments in `config/config.yaml` for the full set of fields (`base_url`, `project_id`, `files`, etc.). When enabled, the workflow registers the run automatically on both success and failure via `onsuccess`/`onerror` hooks.

The FASTA file under `config[pingpong][fasta]` should contain the TE consensus/repeat sequences of interest for the ping-pong analysis and the TE coverage tracks. Example contents:

<details>
  <summary>Example TE fasta file</summary>

```fasta
>L1MdTf_III
agcagcagcggtcgccatctkggttccgggactccgcgggacctaggaaattagtctgaa
caggttagagggtgcgccagagaaccggacagcttctgggacgggcggaagcacagagcc
gctgaggcagcacccttggcgggccgcagacagccggccaccatccggaccagaggacag
gtgtccgcctggcttgggaggcggcctcagcctcagcagcagcggtcgccatctkggttc
cgggactcmgcrgracytaggaaattagtctgaacaggttagagggtgcgccagagaacc
ggacagcttctgggacgggcagaagcacagagccgctgaggcagcascctkggcgggccg
cagacarccggccaccatccggaccagaggacaggtgtccgcctggcttgggaggcggcc
tcagcctcagcagcagcggtcgccatctkggttccgggactccgcggracytaggaaatt
agtctgaacaggttagagggtgcgccagagaaccggacagcttctgggacgggcagaagc
acagagccgctgaggcagcacccttggcgggccgcagacagccggccaccgtccggacca
gaggacaggtgtccgcctggctcgggaggcggcctcagcctcaggagcagcggtcrccat
cttggttccaggactccctggaacttaggawtttagtctgcacaggtgagagtctgcacc
acagaagctgacagcttctgggaactgccaaagcaacacagcttctgagagaggccctgt
tttgggccttcttcttcgaccaggaggaggtccaaaaacaagatatctgcgcaccttccc
tgtaagagagcttgccagcagagagtgctctgagcactgaaactcagaggagagaatctg
tctcccaggtctgctgagagacggtaacagaatcaccagaagaacaatctctaaacagag
tcaactataactactaactccagagattaccagatggcgaaaggtaaacgtaggaatcct
actaacaggaaccaagaccactcaccatcatcagaacccagcactcccacttcgcccagt
ccagggcaccccaacacacccgaaaacctagacctagatttaaaagcatatctcatgatg
```

</details>

To store frequently used `Snakemake` command line values, set up a `config.yaml` in `~/.config/snakemake/standard/`:

```yaml
cores: 32
latency-wait: 20
use-conda: True # recommended
rerun-incomplete: True
printshellcmds: True
cache: False
show-failed-logs: True
use-apptainer: True #recommended
keep-going: True
```

To run the workflow:

```shell
$ snakemake --profile /home/niek/.config/snakemake/standard/
```

### Bigwig normalisation

Before the TE-specific bigWigs are generated, the ping-pong BAM is expanded back to per-read alignments (`{sample}.expanded.bam`, undoing the earlier sequence collapsing) and then length-filtered to `pingpong.min_len`–`pingpong.max_len` (`{sample}.filtered.bam`), so the TE coverage tracks and ping-pong analysis are always restricted to the same read-length range.

`bigwig.normalize_using` controls how those coverage tracks are scaled. It applies to the TE-specific, strand-split bigWigs used for the ping-pong TE coverage tracks (`results/te_bigwig/`, generated by `workflow/scripts/stranded_te_bigwig.py`). Seven values are supported:

| Value        | Description                                                                                                                                                                                                                                                                                                             |
| ------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `RPKM`       | Reads Per Kilobase per Million mapped reads — deepTools' standard length- and depth-normalised coverage.                                                                                                                                                                                                                |
| `CPM`        | Counts Per Million mapped reads — depth normalisation only, no length term.                                                                                                                                                                                                                                             |
| `BPM`        | Bins Per Million mapped reads (TPM-like) — normalises so bin values sum to the same total across samples.                                                                                                                                                                                                               |
| `RPGC`       | Reads Per Genomic Content (1x normalisation) — scales coverage to an average of 1x across the effective genome size. Requires `bigwig.read_length` to pick the correct precomputed effective genome size (50/75/100/150/200) for the configured `genome`.                                                               |
| `None`       | No normalisation; raw read counts per bin.                                                                                                                                                                                                                                                                              |
| `totalReads` | Custom method (not a deepTools built-in). Computes a scale factor as `1 / total_reads * 1e6` from the TEsmall `count_summary.txt` for that sample (i.e. reads per million of the _entire_ small-RNA library, not just TE-mapped reads), then applies it via `bamCoverage --normalizeUsing None --scaleFactor <factor>`. |
| `miRNAReads` | Same approach as `totalReads`, but the scale factor is `1 / total_miRNA_reads * 1e6`, using only reads TEsmall annotated as `miRNA`. Useful for normalising against a stable miRNA reference pool when total small-RNA composition varies a lot between samples/conditions.                                             |

The reverse-strand track's scale factor is negated (multiplied by -1) so forward- and reverse-strand coverage render on opposite sides of zero when viewed together (e.g. in the `te_tracks` pyGenomeTracks plots or a genome browser). The resolved scale factor is also written out per sample/strand for record-keeping.

> **Note:** the miRNA bigWigs under `results/bigwig/` (`mirna_bigwig` rule) currently only support the five standard deepTools values (`RPKM`/`CPM`/`BPM`/`RPGC`/`None`) — `totalReads`/`miRNAReads` are not yet implemented for that rule and will cause `bamCoverage` to error if selected there.

> **Note:** the miRNA bigWig targets (`results/bigwig/mirna_*.bw`) are currently commented out of the default `rule all` targets in the `Snakefile`, so they are not generated by a normal run even though the rules producing them still exist in `bigwig.smk`.

`workflow/scripts/strand_summary_table.py` (`strand_summary_table` rule) collects, per sample and per TE, the total mapped small-RNA read count (from TEsmall's `count_summary.txt`), the raw and normalised forward/reverse ping-pong read counts, and the number of ping-pong pairs analysed, into a single table at `results/pingpong/strand_summary.csv`. The normalisation reuses the same scale factor written out by `stranded_te_bigwig.py`, so it's only meaningful when `bigwig.normalize_using` is `totalReads` or `miRNAReads` — for the standard deepTools methods that scale factor is just `±1`, so the "normalised" columns equal the raw counts in that case.

## Output

```shell
results/
├── bigwig                              # only generated if the mirna_bigwig targets are re-enabled in the Snakefile
│   ├── mirna_KO_1.bw
│   ├── mirna_KO_2.bw
│   ├── mirna_KO_average.bw
│   ├── mirna_WT_1.bw
│   ├── mirna_WT_2.bw
│   └── mirna_WT_average.bw
├── deseq2
│   └── dds_full.rds
├── fasta
│   ├── KO_1.collapsed.fasta
│   ├── KO_1.rrna_mirna_removed.fasta
│   ├── KO_1.rrna_removed.fasta
│   ├── WT_1.collapsed.fasta
│   ├── WT_1.rrna_mirna_removed.fasta
│   └── WT_1.rrna_removed.fasta
├── length_distribution
│   ├── KO_1_length_distribution.txt
│   └── WT_1_length_distribution.txt
├── mirna
│   ├── KO_1.bam
│   ├── KO_1.fasta
│   ├── WT_1.bam
│   └── WT_1.fasta
├── pingpong
│   ├── KO_1.bam
│   ├── KO_1.expanded.bam
│   ├── KO_1.filtered.bam
│   ├── KO_1.csv
│   ├── KO_1_nt_bias.csv
│   ├── WT_1.bam
│   ├── WT_1.expanded.bam
│   ├── WT_1.filtered.bam
│   ├── WT_1.csv
│   ├── WT_1_nt_bias.csv
│   └── strand_summary.csv
├── plots
│   ├── length_distribution_total_count_normalisation.csv
│   ├── length_distribution_total_count_normalisation.pdf
│   ├── length_distribution_tesmall.csv
│   ├── length_distribution_tesmall.pdf
│   ├── pca_plot.pdf
│   ├── pingpong.csv
│   ├── pingpong.pdf
│   ├── pirna_sequence_bias
│   │   ├── L1MdTf_III_bargraph.pdf
│   │   ├── L1MdTf_III_frequencies.csv
│   │   └── L1MdTf_III_logo.pdf
│   ├── rna_counts.pdf
│   ├── te_pos_coverage
│   │   └── L1MdTf_III.pdf
│   ├── te_small
│   │   ├── bargraph.pdf
│   │   ├── piechart.csv
│   │   ├── piechart.pdf
│   │   ├── stacked_bar_te.csv
│   │   ├── stacked_bar_te.pdf
│   │   ├── te_cpm.csv
│   │   ├── te_cpm.pdf
│   │   └── te_top_families.pdf
│   ├── te_tracks
│   │   ├── L1MdTf_III.ini
│   │   └── L1MdTf_III.pdf
│   └── KO_vs_Control
│       ├── results.csv
│       └── volcano_plot.pdf
├── rrna
│   ├── KO_1.bam
│   ├── KO_1.fasta
│   ├── WT_1.bam
│   └── WT_1.fasta
├── te_bigwig
│   ├── KO_1_fwd.bw
│   ├── KO_1_rev.bw
│   ├── KO_1_scale_factor_fwd.txt
│   ├── KO_1_scale_factor_rev.txt
│   ├── KO_average_fwd.bw
│   ├── KO_average_rev.bw
│   ├── WT_1_fwd.bw
│   ├── WT_1_rev.bw
│   ├── WT_1_scale_factor_fwd.txt
│   ├── WT_1_scale_factor_rev.txt
│   ├── WT_average_fwd.bw
│   └── WT_average_rev.bw
├── te_small
│   └── KO_1
│       ├── count_summary.txt
│       ├── KO_1.anno.rlen.info
│       └── report.html
└── trimmed
    ├── KO_1_qc.txt
    └── WT_1_qc.txt
```

## Citation

If you find this workflow useful, please consider citing it:

> Niek Wit. (2026). niekwit/smallRNA-seq: v0.4.0 (v0.4.0). Zenodo. https://doi.org/10.5281/zenodo.19258358

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
