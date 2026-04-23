log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

rlen_files <- snakemake@input[["rlen"]]
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition

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

df_summary <- df |>
  group_by(condition, feature) |>
  summarise(total_count = sum(count), .groups = "drop") |>
  group_by(condition) |>
  mutate(perc_count = (total_count / sum(total_count)) * 100) |>
  ungroup()

other_conditions <- setdiff(unique(df_summary$condition), reference_condition)
df_summary$condition <- factor(
  df_summary$condition,
  levels = c(reference_condition, other_conditions)
)

p <- ggplot(
  df_summary,
  aes(x = "", y = perc_count, fill = feature)
) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = pi / 2) +
  facet_wrap(~condition) +
  labs(title = NULL, fill = "Feature") +
  theme_void() +
  theme(
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  geom_text(
    aes(label = paste0(round(perc_count, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_brewer(palette = "Set3")

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 8,
  height = 5
)

write_csv(df_summary, snakemake@output[["csv"]])
