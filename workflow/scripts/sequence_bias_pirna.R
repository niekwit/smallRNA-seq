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
comparison <- snakemake@wildcards[["comparison"]]
reference_condition <- str_split(comparison, "_vs_")[[1]][2]
test_condition <- str_split(comparison, "_vs_")[[1]][1]
te <- snakemake@wildcards[["te"]]
all_samples <- snakemake@config$samples

sequence_lists <- list()

# Read sequences
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

# Calculate the frequency of the first and tenth nucleotide
freq_df <- data.frame(
  Condition = character(),
  Nucleotide = character(),
  Frequency = numeric(),
  Position = character()
)

for (condition in c(reference_condition, test_condition)) {
  seqs <- sequence_lists[[condition]]
  first_nt <- substr(seqs, 1, 1)
  first_nt_table <- table(first_nt)
  first_nt_freq <- first_nt_table / sum(first_nt_table)

  temp_df <- data.frame(
    Condition = condition,
    Nucleotide = names(first_nt_freq),
    Frequency = as.numeric(first_nt_freq),
    Position = "1"
  )

  freq_df <- rbind(freq_df, temp_df)

  tenth_nt <- substr(seqs, 10, 10)
  tenth_nt_table <- table(tenth_nt)
  tenth_nt_freq <- tenth_nt_table / sum(tenth_nt_table)
  temp_df <- data.frame(
    Condition = condition,
    Nucleotide = names(tenth_nt_freq),
    Frequency = as.numeric(tenth_nt_freq),
    Position = "10"
  )
  freq_df <- rbind(freq_df, temp_df)
}

# Set factor levels for consistent plotting order
freq_df$Condition <- factor(
  freq_df$Condition,
  levels = c(reference_condition, test_condition)
)
freq_df$Nucleotide <- factor(
  freq_df$Nucleotide,
  levels = c("C", "A", "G", "U")
)
freq_df$Position <- factor(
  freq_df$Position,
  levels = c("1", "10")
)

freq_df$te <- te

# Plot as stacked bar graph
bar_pdf <- snakemake@output[["bar"]]
p_bar <- ggplot(freq_df, aes(x = Condition, y = Frequency, fill = Nucleotide)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = NULL,
    y = "Nucleotide frequency",
    title = paste0("Nucleotide frequency at\npositions 1 and 10 for ", te)
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
    labeller = labeller(Position = c("1" = "Position 1", "10" = "Position 10"))
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
