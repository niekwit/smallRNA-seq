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

  condition_name <- snakemake@config$samples[[sample_name]]
  read_delim(
    count_file,
    col_names = FALSE,
    show_col_types = FALSE
  ) %>%
    transmute(
      sample = sample_name,
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

# Combine into single data frame and filter for lengths 18-34
df <- bind_rows(length_distributions) %>%
  filter(length >= 18 & length <= 34)

if (snakemake@config$length_distribution$apply_mirna_correction) {
  # Apply miRNA-based normalization
  print("Applying miRNA-based normalization")
  # Based on:
  # https://link.springer.com/article/10.15252/embj.201386855
  mirna_min <- snakemake@config$length_distribution$mirna_min
  mirna_max <- snakemake@config$length_distribution$mirna_max

  df <- df %>%
    group_by(sample) %>%
    mutate(
      mirna_sum = sum(count[length >= mirna_min & length <= mirna_max])
    ) %>%
    # Normalize every length's count by the miRNA sum of that sample
    mutate(
      mi_norm_frequency = count / mirna_sum
    ) %>%
    # Calculate mean and SD based on miRNA-normalized values
    group_by(condition, length) %>%
    mutate(
      mi_norm_condition_avg = mean(mi_norm_frequency),
      mi_norm_sd = sd(mi_norm_frequency)
    ) %>%
    ungroup()

  # Set plotting and output columns
  y_value <- "mi_norm_condition_avg"
  y_sd <- "mi_norm_sd"
  y_label <- "Fraction (miRNA-Normalized)"
} else {
  print("No miRNA-based normalization applied")
  df <- df %>%
    # Sum of counts per sample
    group_by(sample) %>%
    mutate(sample_sum = sum(count)) %>%
    # Calculate frequency and SD for each condition and length
    group_by(length) %>%
    mutate(sample_frequency = count / sample_sum) %>%
    ungroup() %>%
    # Calculate SD per length across samples of the same condition
    group_by(condition, length) %>%
    mutate(
      condition_frequency = mean(sample_frequency),
      sd = sd(sample_frequency)
    ) %>%
    ungroup()

  # Set plotting and output columns
  y_value <- "condition_frequency"
  y_sd <- "sd"
  y_label <- "Fraction (Total Count Normalized)"
}

# Set colour palette
conditions <- unique(df$condition)
reference_condition <- snakemake@config$reference_condition
other_conditions <- setdiff(conditions, reference_condition)
new_levels <- c(reference_condition, other_conditions)
df$condition <- factor(df$condition, levels = new_levels)
if (length(conditions) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- RColorBrewer::brewer.pal(n = length(conditions), name = "Set3")
}

# Create bar plot with error bars
p <- ggplot(
  df,
  aes(x = factor(length), y = .data[[y_value]], fill = condition)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = .data[[y_value]] - .data[[y_sd]],
      ymax = .data[[y_value]] + .data[[y_sd]]
    ),
    width = 0.2,
    color = "black",
    position = position_dodge(width = 0.9)
  ) +
  labs(
    x = "Length (nt)",
    y = y_label
  ) +
  theme_cowplot() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save plot
ggsave(filename = snakemake@output[["pdf"]], plot = p, width = 6, height = 4)

# Select relevant columns for the final CSV
df_output <- df %>%
  select(
    condition,
    length,
    any_of(c(
      "condition_frequency",
      "sd",
      "mi_norm_condition_avg",
      "mi_norm_sd"
    ))
  ) %>%
  distinct() %>%
  arrange(condition, length)

# Save to CSV
write_csv(df_output, snakemake@output[["csv"]])
