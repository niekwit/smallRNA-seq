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
te <- snakemake@wildcards[["te"]]
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition
other_conditions <- setdiff(unique(unlist(all_samples)), reference_condition)
conditions <- c(reference_condition, other_conditions)

sequence_lists <- list()

# Read sequences for sequence logo
for (condition in conditions) {
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

  all_sense_seqs <- c()
  all_antisense_seqs <- c()

  for (bam in condition_bam_files) {
    bam_data <- BamFile(bam)
    data <- scanBam(bam_data)[[1]]

    # Identify reads mapping to the specified TE and extract all fields
    # together so indices stay in sync through all filtering steps.
    te_idx <- which(as.character(data$rname) == te)
    seqs <- as.character(data$seq[te_idx])
    strands <- as.character(data$strand[te_idx])
    seq_names <- as.character(data$qname[te_idx])

    # Keep only piRNA-length reads (24-32 nt)
    len_ok <- nchar(seqs) >= 24 & nchar(seqs) <= 32
    seqs <- seqs[len_ok]
    strands <- strands[len_ok]
    seq_names <- seq_names[len_ok]

    # Remove reads containing Ns
    no_n <- !grepl("N", seqs)
    seqs <- seqs[no_n]
    strands <- strands[no_n]
    seq_names <- seq_names[no_n]

    # BAM stores minus-strand reads as the reverse complement of the
    # original piRNA (same issue as the Python ping-pong script). RC
    # them back so all sequences run 5'->3' in their original orientation,
    # ensuring position 1 carries U and position 10 carries A.
    minus_idx <- which(strands == "-")
    if (length(minus_idx) > 0) {
      seqs[minus_idx] <- as.character(
        reverseComplement(DNAStringSet(seqs[minus_idx]))
      )
    }

    # Trim to first 10 nt and convert T->U for RNA representation
    seqs <- substr(seqs, 1, 10)
    seqs <- str_replace_all(seqs, "T", "U")

    # Replicate each sequence by its original read count, encoded in
    # the collapsed FASTA qname as "readname-count". Split by strand
    # before replication so sense and antisense stay separate.
    counts <- as.integer(str_split_fixed(seq_names, "-", 2)[, 2])

    sense_mask <- strands == "+"
    antisense_mask <- strands == "-"

    all_sense_seqs <- c(
      all_sense_seqs,
      rep(seqs[sense_mask], times = counts[sense_mask])
    )
    all_antisense_seqs <- c(
      all_antisense_seqs,
      rep(seqs[antisense_mask], times = counts[antisense_mask])
    )
  }
  sequence_lists[[condition]] <- list(
    sense = all_sense_seqs,
    antisense = all_antisense_seqs
  )
}

# Generate sequence logos: N conditions x 2 strands arranged in a grid
# (columns = conditions, rows = sense / antisense).
make_logo <- function(seqs, title) {
  ggseqlogo(seqs, method = "bits", seq_type = "rna") +
    ggtitle(title)
}

n_conditions <- length(conditions)
sense_logos <- lapply(conditions, function(cond) {
  make_logo(sequence_lists[[cond]]$sense, paste0(cond, "\nsense"))
})
antisense_logos <- lapply(conditions, function(cond) {
  make_logo(sequence_lists[[cond]]$antisense, paste0(cond, "\nantisense"))
})
logo_plots <- c(sense_logos, antisense_logos)

pdf(file = snakemake@output[["logo"]], width = 3 * n_conditions, height = 6)
gridExtra::grid.arrange(grobs = logo_plots, ncol = n_conditions)
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
freq_df <- lapply(
  nt_bias_files,
  read_csv,
  show_col_types = FALSE,
  col_types = cols(count = col_double())
) |>
  bind_rows() |>
  filter(repeat_id == te) |>
  mutate(
    Condition = sample_to_condition[sample],
    Nucleotide = str_replace_all(nucleotide, "T", "U"),
    Position = ifelse(strand == "antisense", "1", "10")
  ) |>
  group_by(Condition, Position, Nucleotide) |>
  summarise(count = sum(count), .groups = "drop") |>
  group_by(Condition, Position) |>
  mutate(Frequency = count / sum(count)) |>
  ungroup() |>
  mutate(
    Condition = factor(Condition, levels = conditions),
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
    ncol = 2,
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
  width = 8,
  height = 8
)

# Save frequency data to CSV
write_csv(
  freq_df,
  snakemake@output[["csv"]]
)
