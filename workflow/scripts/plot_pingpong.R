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
  temp_df <- read_delim(
    count_file,
    show_col_types = FALSE,
    col_types = cols(distance = col_double(), count = col_double())
  ) %>%
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
  colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = length(conditions) - 1, name = "Set1")
  )
}

# Set factor levels for conditions
reference_condition <- snakemake@config[["reference_condition"]]
other_conditions <- setdiff(conditions, reference_condition)
new_levels <- c(reference_condition, other_conditions)
df$condition <- factor(df$condition, levels = new_levels)
df_summary$condition <- factor(df_summary$condition, levels = new_levels)

# Set factor levels for repeat_id
repeat_ids <- unique(df$repeat_id)
df$repeat_id <- factor(df$repeat_id, levels = repeat_ids)
df_summary$repeat_id <- factor(df_summary$repeat_id, levels = repeat_ids)

# Plot
plot_type <- snakemake@config$pingpong$plot$type
bar_pos <- if (plot_type == "bar") {
  snakemake@config$pingpong$plot$bar_position
} else {
  NULL
}
square_panels <- plot_type == "line" || identical(bar_pos, "overlap")

if (plot_type == "line") {
  p <- ggplot(
    df_summary,
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
    scale_color_manual(values = colours)
} else {
  bar_pos <- snakemake@config$pingpong$plot$bar_position
  if (bar_pos == "grouped") {
    position_spec <- position_dodge(width = 0.9)
    alpha_val <- 1
  } else {
    position_spec <- "identity"
    alpha_val <- 0.6
  }
  p <- ggplot(
    df_summary,
    aes(x = factor(distance), y = fraction_per_condition, fill = condition)
  ) +
    geom_bar(
      stat = "identity",
      position = position_spec,
      color = "black",
      linewidth = 0.25,
      alpha = alpha_val
    ) +
    geom_errorbar(
      aes(
        ymin = fraction_per_condition - sd_fraction_condition,
        ymax = fraction_per_condition + sd_fraction_condition
      ),
      width = 0.2,
      color = "black",
      position = position_spec
    ) +
    theme_cowplot(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = colours) +
    scale_x_discrete(breaks = function(x) x[as.integer(x) %% 5 == 0]) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
}

# Check if multiple repeat_ids are present
# If so, facet by repeat_id
if (length(unique(df$repeat_id)) > 1) {
  p <- p +
    facet_wrap(~repeat_id, ncol = 3, scales = "free") +
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

if (square_panels) {
  p <- p + theme(aspect.ratio = 1)
}

# Save plot
n_repeats <- length(unique(df$repeat_id))
n_rows <- ceiling(n_repeats / 3)
height <- if (n_repeats > 1) n_rows * 4 else 4
width <- if (n_repeats > 1) {
  if (square_panels) 3 * 4 else 20
} else {
  if (square_panels) 4 else 5
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
