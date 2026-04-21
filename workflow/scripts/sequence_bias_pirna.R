# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(ggseqlogo)
library(Rsamtools)
library(cowplot)

# Get data
bam_files <- snakemake@input[["bam"]]
nt_bias_files <- snakemake@input[["nt_bias"]]
comparison <- snakemake@wildcards[["comparison"]]
reference_condition <- str_split(comparison, "_vs_")[[1]][2]
test_condition <- str_split(comparison, "_vs_")[[1]][1]
te <- snakemake@wildcards[["te"]]
all_samples <- snakemake@config$samples

sequence_lists <- list()

# Read sequences for sequence logo
for (condition in c(reference_condition, test_condition)) {
  # Get sample names for this condition from Snakemake config
  samples <- names(all_samples)[all_samples == condition]

  # Get BAM files that contain these samples
  condition_bam_files <- bam_files[
    sapply(samples, function(s) {
      grepl(s, bam_files)
    }) %>%
      rowSums() >
      0
  ]

  all_seqs <- c()

  for (bam in condition_bam_files) {
    # Load BAM file
    bam_data <- BamFile(bam)

    # Get the sequences from the BAM file as named character vectors
    data <- scanBam(bam_data)[[1]]

    seqs <- as.character(data$seq)
    names(seqs) <- as.character(data$rname)

    # Keep only sequences that map to the specified TE
    seqs <- seqs[names(seqs) == te]

    # Remove sequences that are not 24-32 nt long
    seqs <- seqs[nchar(seqs) >= 24 & nchar(seqs) <= 32]

    # Remove sequences that contain any number of Ns
    seqs <- seqs[!grepl("N", seqs)]

    # Trim sequences to 10 nt
    seqs <- substr(seqs, 1, 10)

    # Convert T to U for RNA representation
    seqs <- str_replace_all(seqs, "T", "U")

    # Get the names of all sequences
    # As sequences were collapsed in the BAM, we need to replicate
    # each sequence according to its count
    seq_names <- data$qname[which(names(seqs) == te)]

    # Sequence names are in the format: seqname-count
    counts <- as.integer(str_split_fixed(seq_names, "-", 2)[, 2])

    # Replicate sequences according to their counts
    seqs <- rep(seqs, times = counts)

    all_seqs <- c(all_seqs, seqs)
  }
  sequence_lists[[condition]] <- all_seqs
}

# Generate sequence logos
p1 <- ggseqlogo(
  sequence_lists[[reference_condition]],
  method = "bits",
  seq_type = "rna"
)

p2 <- ggseqlogo(
  sequence_lists[[test_condition]],
  method = "bits",
  seq_type = "rna"
)

pdf(file = snakemake@output[["logo"]])
gridExtra::grid.arrange(
  p1 + ggtitle(reference_condition),
  p2 + ggtitle(test_condition),
  ncol = 1
)
dev.off()

# Build a sample -> condition lookup from config
sample_to_condition <- setNames(
  unlist(all_samples),
  names(all_samples)
)

# Read nucleotide bias CSVs (ping-pong pairs at distance 10 only) and
# compute per-condition nucleotide frequencies at:
#   position 1  -> antisense read position 1
#   position 10 -> sense read position 10
freq_df <- lapply(nt_bias_files, read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  filter(repeat_id == te) |>
  mutate(
    Condition = sample_to_condition[sample],
    Nucleotide = str_replace_all(nucleotide, "T", "U"),
    Position = ifelse(strand == "antisense", "1", "10")
  ) |>
  filter(Condition %in% c(reference_condition, test_condition)) |>
  group_by(Condition, Position, Nucleotide) |>
  summarise(count = sum(count), .groups = "drop") |>
  group_by(Condition, Position) |>
  mutate(Frequency = count / sum(count)) |>
  ungroup() |>
  mutate(
    Condition = factor(
      Condition,
      levels = c(reference_condition, test_condition)
    ),
    Nucleotide = factor(Nucleotide, levels = c("C", "A", "G", "U")),
    Position = factor(Position, levels = c("1", "10")),
    te = te
  )

# Plot as stacked bar graph
bar_pdf <- snakemake@output[["bar"]]
p_bar <- ggplot(
  freq_df,
  aes(x = Condition, y = Frequency, fill = Nucleotide)
) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = NULL,
    y = "Nucleotide frequency",
    title = paste0(
      "Nucleotide frequency at\nping-pong positions for ",
      te
    )
  ) +
  theme_cowplot(16) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = c(
      "A" = "#109648",
      "U" = "#d62839",
      "C" = "#255c99",
      "G" = "#f7b32b"
    )
  ) +
  facet_wrap(
    ~Position,
    ncol = 1,
    labeller = labeller(
      Position = c(
        "1" = "Position 1 (antisense)",
        "10" = "Position 10 (sense)"
      )
    )
  )

ggsave(
  filename = bar_pdf,
  plot = p_bar,
  width = 4,
  height = 8
)

# Save frequency data to CSV
write_csv(
  freq_df,
  snakemake@output[["csv"]]
)
