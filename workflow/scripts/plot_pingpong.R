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
  fraction_of_te = double()
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
    ungroup()
  df <- bind_rows(df, temp_df)
}

# Create summary statistics per distance and condition
df_summary <- df %>%
  group_by(repeat_id, distance, condition) %>%
  summarise(
    fraction_per_condition = mean(fraction_of_te, na.rm = TRUE),
    sd_fraction_condition = sd(fraction_of_te, na.rm = TRUE),
    .groups = "drop"
  )

# Add summary stats to main df
df <- df %>%
  left_join(
    df_summary,
    by = c("repeat_id", "distance", "condition")
  )

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

# Set factor levels for repeat_id
repeat_ids <- unique(df$repeat_id)
df$repeat_id <- factor(df$repeat_id, levels = repeat_ids)

# Plot
p <- ggplot(
  df,
  aes(x = distance, y = fraction_per_condition, color = condition)
) +
  geom_line(aes(group = condition)) +
  geom_errorbar(
    aes(
      ymin = fraction_per_condition - sd_fraction_condition,
      ymax = fraction_per_condition + sd_fraction_condition
    ),
    width = 0.4,
    color = "black",
    linewidth = 0.25
  ) +
  theme_cowplot(12) +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours)

# Check if multiple repeat_ids are present
# If so, facet by repeat_id
if (length(unique(df$repeat_id)) > 1) {
  p <- p +
    facet_wrap(~repeat_id, ncol = 2) +
    labs(
      x = "Distance between 5' ends (nt)",
      y = "Frequency"
    )
} else {
  p <- p +
    labs(
      title = unique(df$repeat_id),
      x = "Distance between 5' ends (nt)",
      y = "Frequency"
    )
}

# Save plot
width <- if (length(unique(df$repeat_id)) > 1) {
  8
} else {
  5
}
height <- if (length(unique(df$repeat_id)) > 1) {
  ceiling(length(unique(df$repeat_id)) / 2) * 3
} else {
  4
}
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = width,
  height = height
)

# Save data as CSV
write_csv(
  df,
  file = snakemake@output[["csv"]]
)
