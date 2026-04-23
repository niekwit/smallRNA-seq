log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

rlen_files <- snakemake@input[["rlen"]]
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition
minlen <- snakemake@params[["minlen"]]
maxlen <- snakemake@params[["maxlen"]]

sample_to_condition <- setNames(unlist(all_samples), names(all_samples))

df <- lapply(rlen_files, function(f) {
  sample_name <- basename(dirname(f))
  read_delim(f, show_col_types = FALSE) |>
    mutate(
      sample = sample_name,
      condition = sample_to_condition[sample_name]
    )
}) |>
  bind_rows() |>
  pivot_longer(
    cols = -c(rlen, sample, condition),
    names_to = "feature",
    values_to = "count"
  )

df_summary <- df %>%
  group_by(condition, feature, rlen) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  group_by(condition, rlen) %>%
  mutate(perc_count = (total_count / sum(total_count)) * 100) %>%
  ungroup() %>%
  filter(rlen >= minlen & rlen <= maxlen)

other_conditions <- setdiff(unique(df_summary$condition), reference_condition)
df_summary$condition <- factor(
  df_summary$condition,
  levels = c(reference_condition, other_conditions)
)

p <- ggplot(
  df_summary,
  aes(x = rlen, y = perc_count, fill = feature)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    colour = "black",
    linewidth = 0.25
  ) +
  facet_wrap(~condition) +
  labs(
    x = "Read length (nt)",
    y = "Percentage of total reads (%)",
    fill = "Feature"
  ) +
  theme_cowplot() +
  theme(
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(minlen, maxlen, by = 2),
    guide = guide_axis(minor.ticks = TRUE),
    expand = c(0, 0)
  )

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 10,
  height = 5
)
