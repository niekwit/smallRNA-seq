library(tidyverse)
library(cowplot)

# Get fasta files
fasta_files <- snakemake@input[["fasta"]]

# For each fasta extract length and counts
length_distributions <- lapply(fasta_files, function(fasta_file) {
  condition_name <- str_extract(basename(fasta_file), "^[^_]+")
  sample_name <- str_replace(basename(fasta_file), "_counts.fasta$", "")
  seqs <- read_lines(fasta_file)
  headers <- seqs[seq(1, length(seqs), by = 2)]

  # Split header with : as delimiter into 3 parts: sequence, length, count
  split_headers <- str_split_fixed(headers, ":", 3)

  # Store lengths and counts in vectors
  lengths <- as.integer(split_headers[, 2])
  counts <- as.integer(split_headers[, 3])

  # Create data frame with lengths and counts
  tmp <- data.frame(
    length = lengths,
    count = counts,
    sample = sample_name,
    condition = condition_name,
    stringsAsFactors = FALSE
  ) %>%
    group_by(length) %>%
    mutate(count = sum(count)) %>%
    unique() %>%
    arrange(length) %>%
    ungroup()
})

# Combine all data frames into one
df <- bind_rows(length_distributions) %>%
  # Only keep lengths from 18 to 32
  filter(length >= 18 & length <= 32) %>%
  # Sum of counts per sample
  group_by(sample) %>%
  mutate(sample_sum = sum(count)) %>%
  # Calculate frequency and SEM for each condition and length
  group_by(length) %>%
  mutate(sample_frequency = count / sample_sum) %>%
  ungroup() %>%
  # Calculate SEM per length across samples of the same condition
  group_by(condition, length) %>%
  mutate(
    condition_frequency = mean(sample_frequency),
    sem = sd(sample_frequency) / sqrt(n())
  ) %>%
  ungroup()

# Create line plot with error bars
p <- ggplot(df, aes(x = length, y = condition_frequency, color = condition)) +
  geom_line(aes(group = sample)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = condition_frequency - sem, ymax = condition_frequency + sem),
    width = 0.2,
    color = "black"
  ) +
  labs(
    x = "Length (nt)",
    y = "Frequency"
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  scale_x_continuous(
    breaks = seq(18, 32, by = 2),
    guide = guide_axis(minor.ticks = TRUE)
  )

# Save plot
ggsave(filename = snakemake@output[["pdf"]], plot = p, width = 6, height = 4)

# Save data as CSV
# Only keep and save unique condition, length, condition_frequency, sem
df_output <- df %>%
  # Add columns with sums per condition and length, and total sum
  group_by(condition) %>%
  mutate(
    condition_sum = sum(count)
  ) %>%
  group_by(length) %>%
  mutate(
    length_sum = sum(count)
  ) %>%
  select(
    condition,
    length,
    length_sum,
    condition_sum,
    condition_frequency,
    sem
  ) %>%
  unique() %>%
  arrange(condition, length)

# Save to CSV
write_csv(df_output, snakemake@output[["csv"]])
