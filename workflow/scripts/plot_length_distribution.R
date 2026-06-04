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
  method <- snakemake@config$length_distribution$mirna_correction_method

  if (method == "align") {
    # Normalise by total aligned miRNA counts (reads per million miRNA).
    # miRNA FASTA headers encode read count as ">readID-count"; sum all
    # counts across sequences to get the total miRNA library size per sample.
    print("Applying miRNA-based normalization (align method)")
    mirna_files <- snakemake@input[["mirna_fasta"]]

    mirna_totals <- lapply(mirna_files, function(f) {
      sample_name <- str_replace(basename(f), "\\.fasta$", "")
      headers <- grep("^>", readLines(f), value = TRUE)
      counts <- as.integer(
        str_split_fixed(str_remove(headers, "^>"), "-", 2)[, 2]
      )
      data.frame(sample = sample_name, mirna_total = sum(counts, na.rm = TRUE))
    })
    mirna_totals_df <- bind_rows(mirna_totals)

    df <- df %>%
      left_join(mirna_totals_df, by = "sample") %>%
      group_by(sample) %>%
      mutate(rpm_mirna = count / mirna_total * 1e6) %>%
      group_by(condition, length) %>%
      mutate(
        rpm_mirna_avg = mean(rpm_mirna),
        rpm_mirna_sd = sd(rpm_mirna)
      ) %>%
      ungroup()

    y_value <- "rpm_mirna_avg"
    y_sd <- "rpm_mirna_sd"
    y_label <- "Reads per million miRNA"
  } else {
    # Length-based miRNA normalisation.
    # Based on:
    # https://link.springer.com/article/10.15252/embj.201386855
    print("Applying miRNA-based normalization (length method)")
    mirna_min <- snakemake@config$length_distribution$mirna_min
    mirna_max <- snakemake@config$length_distribution$mirna_max

    df <- df %>%
      group_by(sample) %>%
      mutate(
        mirna_sum = sum(count[length >= mirna_min & length <= mirna_max])
      ) %>%
      mutate(mi_norm_frequency = count / mirna_sum) %>%
      group_by(condition, length) %>%
      mutate(
        mi_norm_condition_avg = mean(mi_norm_frequency),
        mi_norm_sd = sd(mi_norm_frequency)
      ) %>%
      ungroup()

    y_value <- "mi_norm_condition_avg"
    y_sd <- "mi_norm_sd"
    y_label <- "Fraction (miRNA-Normalized)"
  }
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
