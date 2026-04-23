log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

txt_files <- snakemake@input[["txt"]]
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition

sample_to_condition <- setNames(unlist(all_samples), names(all_samples))

df <- lapply(txt_files, function(f) {
  sample_name <- basename(dirname(f))
  read_delim(
    f,
    col_names = c("fid", "ftype", "counts"),
    show_col_types = FALSE
  ) |>
    filter(ftype %in% c("sense_TE", "anti_TE")) |>
    separate(
      fid,
      into = c("Class", "Family", "Subfamily", "Instance"),
      sep = ":"
    ) |>
    mutate(
      Class = case_when(
        Class == "LTR" & grepl("IAP", Subfamily) ~ "IAP",
        TRUE ~ Class
      ),
      counts = as.numeric(counts),
      total_count = sum(counts, na.rm = TRUE),
    ) |>
    group_by(Class) |>
    summarise(
      Class_Count = sum(counts, na.rm = TRUE),
      Total_Reads = first(total_count),
      .groups = "drop"
    ) |>
    mutate(
      Fraction_of_Total = Class_Count / Total_Reads,
      sample = sample_name,
      condition = sample_to_condition[sample_name]
    )
}) |>
  bind_rows()

df_summary <- df |>
  group_by(condition, Class) |>
  summarise(avg_percentage = mean(Fraction_of_Total) * 100, .groups = "drop")

write_csv(df_summary, snakemake@output[["csv"]])

# Collapse minor classes into Other
df_summary <- df_summary |>
  mutate(
    Class = case_when(
      Class %in% c("RC", "Retroposon", "Unspecified", "Unknown", "Satellite") ~
        "Other",
      TRUE ~ Class
    )
  ) |>
  group_by(condition, Class) |>
  summarise(avg_percentage = sum(avg_percentage), .groups = "drop")

other_conditions <- setdiff(unique(df_summary$condition), reference_condition)
df_summary$condition <- factor(
  df_summary$condition,
  levels = c(reference_condition, other_conditions)
)

df_summary$Class <- factor(
  df_summary$Class,
  levels = c("LINE", "SINE", "LTR", "IAP", "DNA", "Other")
)

p <- ggplot(
  df_summary,
  aes(x = condition, y = avg_percentage, fill = Class)
) +
  geom_bar(stat = "identity", colour = "black") +
  labs(
    x = NULL,
    y = "Mapped TE-derived piRNAs (%)",
    fill = NULL
  ) +
  theme_cowplot(14) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(
    values = c(
      "LINE" = "white",
      "SINE" = "#a6cee3",
      "LTR" = "#1f78b4",
      "IAP" = "#174891",
      "DNA" = "#103061",
      "Other" = "grey"
    )
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101))

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 3.5,
  height = 6
)
