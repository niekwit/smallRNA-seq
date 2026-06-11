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
│       ├── config.schema.yaml
│       └── samples.schema.yaml
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
    │   ├── process_reads.yaml
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
    │   ├── pingpong.py
    │   ├── plot_length_distribution.R
    │   ├── plot_pca.R
    │   ├── plot_pingpong.R
    │   ├── plot_rna_counts.R
    │   ├── sequence_bias_pirna.R
    │   ├── te_pos_coverage.R
    │   ├── te_small_bargraph.R
    │   ├── te_small.py
    │   ├── te_small_piechart.R
    │   ├── te_small_stacked_bar_te.R
    │   └── te_small_top_families.R
    └── Snakefile

12 directories, 33 files
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

bigwig:
  read_length: 50 # used to select effective genome size for bamCoverage
  normalize_using: CPM # RPKM, CPM, BPM, RPGC, or None
  extra: "" # additional bamCoverage arguments

pingpong:
  mismatch: 3
  fasta: resources/L1mTf.fa
  window: 30
```

The FASTA file under config[pingpong][fasta] should contain the TE sequences of interest for ping-pong analysis. Example contents:

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
atggtagaggacatcaagaaggactttaataaatcacttaaagaaatacaggagaacact
gctaaagagttacaagtccttaaagaaaaacaggaaaacacaatcaaacaggtagaagtc
cttacagaaaaagaggaaaaaacatacaaacaggtgatggaaatgaacaaaaccatacta
gacctaaaaagggaagtagaaacaataaagaaaactcaaagtgaggcaacactggagata
gaaaccctaggaaagaaatctggaaccatagatttgagcatcagcaacagaatacaagag
atggaagagagaatctcaggtgcagaagattccatagagaacatcggcacaacaatcaaa
gaaaatggaaaatgcaaaaagatcctaactcaaaatatccaggaaatccaggacacaatg
agaagaccaaacctacggataataggagtggatgagaatgaagattttcaactcaaagga
ccagcaaacatcttcaacaaaattattgaagaaaacttcccaaatataaagaaagagata
cctatgaacatacaagaagcctacagaactccaaatagactggaccagaaaagaaattcc
tcccgacacataataatcagaacaacaaatgcactaaataaagatagaatactaaaagca
gtaagggaaaaaggtcaagtaacatacaaaggcaagcctatcagaattacaccagatttt
tcaccagagactatgaaagccagaagagcctggacagatgttatacagacactaagagaa
cacaaattccagcccaggctactatacccagccaaactctcaattaccatagatggagaa
accaaagtattccacgacaaaaccaaattcacacattatctctccacgaatccagccctt
caaaggttaataacagaaaaaaaccaatacaagaacgggaacaatgccctagaaaaaaca
agaaggtaatccctcaacaaacctaaaagaagacagccacaagaacagaatgccaacttt
aacaacaaaaataacaggaagcaacaattacttttccttaatatctcttaacatcaatgg
tctcaactccccaataaaaagacatagactaacaaactggctacacaaacaagacccaac
attttgctgtttacaggagacacatctcagagaaaaagatagacactacctcagaataaa
aggctggaaaacaattttccaagcaaatggtatgaagaaacaagctggagtagccatcct
aatatctgataagattgacttccaacccaaagtcatcaaaaaagacaaggaggggcactt
ygttctcatcaaaggtaaaatcctccaagaggaactctcaattctgaatatctatgctcc
aaatacaagggcagccacattcattaaagaaactttagtaaagctcaaagcacacattgc
acctcacacaataatagtgggagacttcaacacaccactttcaccaatggacagatcatg
gaaacagaaactaaacagggacacactgaaactaacagaagtgatgaaacaaatggatct
gacagatatctacagaacattttatcctaaaacaaaaggatataccttcttctcagcacc
tcatggtaccttctccaaaattgaccacataataggtcacaaaacaggcctcaacagatt
caaaaatattgaaattgtcccatgtatcctatcagatcaccatgcactaaggctgatctt
caataacaaaaaaaataacagaaagccaacactcacgtggaaactgaacaacactcttct
caatgataccttggtcaaggaaggaataaagaaagaaattaaagactttttagagtttaa
tgaaaatgaagccacaacgtacccaaacctttgggacacaatgaaagcatttctaagagg
gaaactcatagctctgagtgcctccatgaagaaacgggagagagcacatactagcagctt
gacaacacatctaaaagctctagaaaaaaaggaagcaaattcacccaagaggagtagacg
gcaggaaataatcaaactcaggggtgaaatcaaccaagtggaaacaagaagaactattca
aagaattaaccaaacgaggagttggttctttgagaaaatcaacaagatagataaaccctt
agctagactcactagagggcacagagacaaaatcctaattaacaaaatcagaactgaaaa
gggagacataacaacagatcctgaagaaatccaaaacaccatcagatccttctacaaaag
gctatactcaacaaaactggaaaacctggacgaaatggacaaatttctggacagatacca
ggtaccaaagttgaatcaggatcaagttgaccttctaaacagtcccatatcccctaaaga
aatagaagcagttataaatagtctcccagccaaaaaaagcccaggaccagacgggtttag
tgcagagttctatcagaccttcaaagaagatctaattccagttctgcacaaactttttca
caagatagaagtagaaagtactctacccaactcattttatgaagccactattactctgat
acctaaaccacagaaagatccaacaaagatagagaacttcagaccaatttctcttatgaa
tatcgatgcaaaaatcctcaataaaattctcgctaaccgaatccaagaacacattaaagc
aatcatccatcctgaccaagtaggttttattccaggratgcagggatggtttaatatacg
aaaatccatcaatgtaatccactatataaacaaactcaaagacaaaaaccacatgatcat
ctcgttrgatgcagaaaaagcatttgacaagatccaacacccattcatgataaaagttct
ggaaagatcaggaattcaaggcccatacctaaacatgataaaagcaatctacagcaaacc
agtagccaacatcaaagtaaatggagagaagctggaagcaatcccactaaaatcagggac
tagacaaggctgcccactttctccctaccttttcaacatagtacttgaagtattagccag
agcaattcgacaacaaaaggagatcaaggggatacaaattggaaaggaggaagtcaaaat
atcactttttgcagatgatatgatagtatatataagtgaccctaaaaattccaccagaga
actcctaaacctgataaacagcttcggtgaagtagctggatataaaattaactcaaacaa
gtcaatggcctttctctacacaaagaataaacaggctgagaaagaaattagggaaacaac
acccttctcaatagtcacaaataatataaaatatctcggcgtgactctaactaaggaagt
gaaagatctgtatgataaaaacttcaagtctctgaagaaagaaattaaagaagatctcag
aagatggaaagatctcccatgctcatggattggcaggatcaayattgtaaaaatggctat
cttgccaaaagcaatctacagattcaatgcaatccccatcaaaattccaactcaattctt
caacgaattagaaggagcaatttgcaaattcatctggaataacaaaaaacctaggatagc
aaaaactcttctcaaggataaaagaacctctggtggaatcaccatgcctgacctaaagct
ttactacagagcaattgtggtaaaaactgcatggtactggtatagagacagacaagtaga
ccaatggaatagaattgaagacccagaaatgaacccacacacctatggtcacttgatctt
cgacaagggagctaaaaccatccagtggaagaaagacagcattttcaacaaatggtgctg
gcacaactggttgttatcatgtagaagaatgcgaatcgatccatacttatctccttgtac
taaggtcaaatctaaatggatcaaagaacttcacataaaaccagagacactgaaacttat
agaggagaaagtggggaaaagccttgaagatatgggcacaggggaaaaattcctgaacag
aacagcaatggcttgtgctgtaagatcgagaattgacaaatgggacctaatgaaactcca
aagtttctgcaaggcaaaagacaccgtcaataagacaaagagaccaccaacagattggga
aaggatctttacctatcctaaatcagataggggactaatatccaacatatataaagaact
caagaaggtggacttcagaaaatcaaayaaccccattaaaaaatggggctcagaactgaa
caaagaattctcacctgaggaataccgaatggcagagaagcacctgaaaaaatgttcaac
atccttaatcatcagggaaatgcaaatcaaaacaaccctgagattccacctcacaccagt
cagaatgtctaagatcaaaaattcaggtgacagcagatgctggcgaggatgtggagaaag
aggaacactcctccattgttggtgggattgcaggcttgtacaaccactctggaaatccgt
ctggcggttcctcagaaaattggacatagtactaccggaggatccagcaatacctctcct
gggcatatatccagaagatgccccaactggtaagaaggacacatgctccactatgttcat
agcagccttatttataatagccagaagctggaaagaacccagatgcccctcaacagagga
atggatacagaaaatgtggtacatctacacaatggagtactactcagctattaaaaagaa
tgaatttatgaaattcctagccaaatggatggacctggagggcatcatcctgagtgaggt
aacacattcacaaaggaactcacacaatatgtactcactgataagtggatattagcccaa
aacctaggatacccaagatataagatacaatttcctaaacacatgaaactcaagaaaaat
gaagactgaagtgtggacactatgcccctccttagaagtgggaacaaaacacccttggaa
ggagttacagagacaaagtttggagctgagatgaaaggatggaccatgtagagactgcct
tatccagggatccaccccataatcagcatccaaacgctgacaccattgcatacactagca
agattttatcgaaaggacccagatgtagctgtctcttgtgagactatgccggggcctagc
aaacacagaagtggatgcccacagtcagctaatggatggatcacagggctcccaatggag
gagctagagaaagtacccaaggagctaaagggatctgcaaccctataggtggatcaacat
tatgaactaaccagtaccccggagctcttgactctagctgcatatgtatcaaaagatggc
ctagtcggccatcactggaaagagaggcccattggacacacaaactttatatgccccaga
acaggggaacgccagggccaaaaagggggagtgggcgggtaggggagtgggggtgggtgg
gtatgggggacttttggtatagcattggaaatgtaaatgagctaaatacctaataaaaaa
tggaaagaaa
```

</details>

To store frequently used `Snakemake` command line values, set up a `config.yaml` in ~/.config/snakemake/standard/:

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

## Output

```shell
results/
├── bigwig
│   ├── mirna_KO_1.bw
│   ├── mirna_KO_2.bw
│   ├── mirna_KO_average.bw
│   ├── mirna_WT_1.bw
│   ├── mirna_WT_2.bw
│   └── mirna_WT_average.bw
├── fasta
│   ├── KO_1.collapsed.fasta
│   ├── KO_1.rrna_mirna_removed.fasta
│   ├── KO_1.rrna_removed.fasta
│   ├── KO_2.collapsed.fasta
│   ├── KO_2.rrna_mirna_removed.fasta
│   ├── KO_2.rrna_removed.fasta
│   ├── WT_1.collapsed.fasta
│   ├── WT_1.rrna_mirna_removed.fasta
│   ├── WT_1.rrna_removed.fasta
│   ├── WT_2.collapsed.fasta
│   ├── WT_2.rrna_mirna_removed.fasta
│   └── WT_2.rrna_removed.fasta
├── length_distribution
│   ├── KO_1_length_distribution.txt
│   ├── KO_2_length_distribution.txt
│   ├── WT_1_length_distribution.txt
│   └── WT_2_length_distribution.txt
├── mirna
│   ├── KO_1.bam
│   ├── KO_1.fasta
│   ├── KO_2.bam
│   ├── KO_2.fasta
│   ├── WT_1.bam
│   ├── WT_1.fasta
│   ├── WT_2.bam
│   └── WT_2.fasta
├── pingpong
│   ├── KO_1.bam
│   ├── KO_1.csv
│   ├── KO_1_nt_bias.csv
│   ├── KO_2.bam
│   ├── KO_2.csv
│   ├── KO_2_nt_bias.csv
│   ├── WT_1.bam
│   ├── WT_1.csv
│   ├── WT_1_nt_bias.csv
│   ├── WT_2.bam
│   ├── WT_2.csv
│   └── WT_2_nt_bias.csv
├── plots
│   ├── length_distribution_total_count_normalisation.csv
│   ├── length_distribution_total_count_normalisation.pdf
│   ├── pingpong.csv
│   ├── pingpong.pdf
│   ├── pirna_sequence_bias
│   │   ├── L1MdTf_III_bargraph.pdf
│   │   ├── L1MdTf_III_frequencies.csv
│   │   └── L1MdTf_III_logo.pdf
│   ├── rna_counts.pdf
│   └── te_pos_coverage
│       └── L1MdTf_III.pdf
├── rrna
│   ├── KO_1.bam
│   ├── KO_1.fasta
│   ├── KO_2.bam
│   ├── KO_2.fasta
│   ├── WT_1.bam
│   ├── WT_1.fasta
│   ├── WT_2.bam
│   └── WT_2.fasta
├── seqs
│   ├── KO_1_seqs.txt
│   ├── KO_2_seqs.txt
│   ├── WT_1_seqs.txt
│   └── WT_2_seqs.txt
└── trimmed
    ├── KO_1_qc.txt
    ├── KO_2_qc.txt
    ├── WT_1_qc.txt
    └── WT_2_qc.txt
```

## Citation

If you find this workflow useful, please consider citing it:

> Niek Wit. (2026). niekwit/smallRNA-seq: v0.4.0 (v0.4.0). Zenodo. https://doi.org/10.5281/zenodo.19258358

## References

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, 10, 33, **2021**. https://doi.org/10.12688/f1000research.29032.2.
