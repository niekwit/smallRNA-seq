# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Count files
count_files <- snakemake@input[["pingpong"]]

# Create empty df with columns repeat_id, distance, count, sample
df <- data.frame(
  repeat_id = character(),
  distance = integer(),
  count = integer(),
  sample = character(),
  condition = character(),
  total_reads_per_te = integer(),
  fraction_of_te = double(),
  avg_fraction_condition = double(),
  sd_fraction_condition = double()
)

# Load each count file and rbind to df
for (count_file in count_files) {
  sample <- basename(count_file) %>% str_replace(".csv$", "")
  condition <- snakemake@config$samples[[sample]]

  temp_df <- read_delim(count_file) %>%
    group_by(repeat_id, sample) %>%
    mutate(
      condition = condition,
      total_reads_per_te = sum(count),
      fraction_of_te = count / total_reads_per_te
    ) %>%
    ungroup() %>%
    # Calculate average fraction of all TEs per distance per sample
    group_by(sample, distance) %>%
    mutate(avg_fraction = mean(fraction_of_te)) %>%
    ungroup() %>%
    # Calculate average fraction of per distance per condition
    group_by(condition, distance) %>%
    mutate(
      avg_fraction_condition = mean(fraction_of_te),
      sd_fraction_condition = sd(fraction_of_te)
    ) %>%
    ungroup()

  df <- bind_rows(df, temp_df)
}

# Prepare colours for plotting
conditions <- unique(df$condition)
if (length(conditions) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- RColorBrewer::brewer.pal(n = length(conditions), name = "Set3")
}

# Set factor levels for conditions
reference_condition <- snakemake@config[["reference_condition"]]
other_conditions <- setdiff(conditions, reference_condition)
new_levels <- c(reference_condition, other_conditions)
df$condition <- factor(df$condition, levels = new_levels)

# Plot
p <- ggplot(
  df,
  aes(x = distance, y = avg_fraction_condition, color = condition)
) +
  geom_line(aes(group = condition)) +
  geom_errorbar(
    aes(
      ymin = avg_fraction_condition - sd_fraction_condition,
      ymax = avg_fraction_condition + sd_fraction_condition
    ),
    width = 0.4,
    color = "black",
    linewidth = 0.25
  ) +
  labs(
    title = snakemake@config[["pingpong"]][["name"]],
    x = "Distance between 5' ends (nt)",
    y = "Frequency"
  ) +
  theme_cowplot(12) +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours)

# Save plot
ggsave(filename = snakemake@output[["pdf"]], plot = p, width = 6, height = 4)

# Save data as CSV
write_csv(
  df,
  file = snakemake@output[["csv"]]
)
