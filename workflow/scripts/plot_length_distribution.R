# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


library(tidyverse)
library(cowplot)

# Get fasta files
count_files <- snakemake@input[["counts"]]

# For each fasta extract length and counts
length_distributions <- lapply(count_files, function(count_file) {
  sample_name <- str_replace(
    basename(count_file),
    "_length_distribution.txt$",
    ""
  )

  condition_name <- str_extract(sample_name, "^[^_]+")

  read_delim(
    count_file,
    col_names = FALSE,
    show_col_types = FALSE
  ) %>%
    transmute(
      sample = X3,
      length = X2,
      count = X1
    ) %>%
    group_by(length) %>%
    mutate(count = sum(count)) %>%
    unique() %>%
    arrange(length) %>%
    ungroup() %>%
    mutate(condition = condition_name)
})

# Combine into single data frame and calculate frequencies and SD
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
    sd = sd(sample_frequency)
  ) %>%
  ungroup()

# Create line plot with error bars
p <- ggplot(df, aes(x = length, y = condition_frequency, color = condition)) +
  geom_line(aes(group = sample)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = condition_frequency - sd, ymax = condition_frequency + sd),
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
# Only keep and save unique condition, length, condition_frequency, sd
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
    sd
  ) %>%
  unique() %>%
  arrange(condition, length)

# Save to CSV
write_csv(df_output, snakemake@output[["csv"]])
