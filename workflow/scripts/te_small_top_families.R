log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

txt_files <- snakemake@input[["txt"]]
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition
n_top <- snakemake@params[["n_top"]]
te_classes <- snakemake@params[["te_classes"]]

sample_to_condition <- setNames(unlist(all_samples), names(all_samples))
other_conditions <- setdiff(unique(unlist(all_samples)), reference_condition)
condition_levels <- c(reference_condition, other_conditions)

df_list <- list()

for (f in txt_files) {
  sample_name <- basename(dirname(f))

  raw <- read_tsv(
    f,
    col_names = c("fid", "ftype", "counts"),
    show_col_types = FALSE
  ) %>%
    filter(ftype %in% c("sense_TE", "anti_TE")) %>%
    mutate(counts = as.numeric(counts))

  total_reads <- sum(raw$counts, na.rm = TRUE)

  aggregate_counts <- function(data, mode_label) {
    data %>%
      separate(
        fid,
        into = c("Class", "Family", "Subfamily", "Instance"),
        sep = ":",
        extra = "drop",
        fill = "right"
      ) %>%
      mutate(
        Subfamily = str_remove(Subfamily, "_copy\\d+.*$|_dup\\d+.*$"),
        Subfamily = str_remove(Subfamily, "_[IVX]+$")
      ) %>%
      group_by(Class, Subfamily) %>%
      summarise(count = sum(counts, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        fraction = count / total_reads,
        sample = sample_name,
        condition = sample_to_condition[[sample_name]],
        mode = mode_label
      )
  }

  df_list[[sample_name]] <- bind_rows(
    aggregate_counts(raw, "All piRNAs"),
    aggregate_counts(filter(raw, ftype == "anti_TE"), "Antisense piRNAs")
  )
}

df <- bind_rows(df_list)

# Filter to user-specified TE classes (or keep all)
if (!(length(te_classes) == 1 && te_classes == "All")) {
  df <- df %>% filter(Class %in% te_classes)
}
top_class_names <- unique(df$Class)

# Top N subfamilies per class (All piRNAs)
top_families <- df %>%
  filter(mode == "All piRNAs") %>%
  group_by(Class, Subfamily) %>%
  summarise(mean_fraction = mean(fraction), .groups = "drop") %>%
  group_by(Class) %>%
  slice_max(mean_fraction, n = n_top) %>%
  select(Class, Subfamily)

df_filtered <- df %>%
  semi_join(top_families, by = c("Class", "Subfamily")) %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    mode = factor(mode, levels = c("All piRNAs", "Antisense piRNAs"))
  )

# Inline reorder_within: encodes Class into the factor level so that
# ggplot2 free_x scales sort each panel independently by abundance.
reorder_within <- function(x, by, within) {
  stats::reorder(paste(x, within, sep = "___"), by, mean)
}
scale_x_reordered <- function(...) {
  ggplot2::scale_x_discrete(labels = function(x) sub("___.*$", "", x), ...)
}

df_summary <- df_filtered %>%
  group_by(condition, Class, Subfamily, mode) %>%
  summarise(
    mean_fraction = mean(fraction),
    sd_fraction = sd(fraction),
    .groups = "drop"
  )

# Ordering key: mean fraction in "All piRNAs" mode, averaged across conditions.
# Used for both panels of each class so All/Antisense share the same x order.
order_df <- df_summary %>%
  filter(mode == "All piRNAs") %>%
  group_by(Class, Subfamily) %>%
  summarise(order_val = mean(mean_fraction), .groups = "drop")

df_summary <- df_summary %>%
  left_join(order_df, by = c("Class", "Subfamily")) %>%
  mutate(Subfamily = reorder_within(Subfamily, -order_val, Class))

df_filtered <- df_filtered %>%
  left_join(order_df, by = c("Class", "Subfamily")) %>%
  mutate(Subfamily = reorder_within(Subfamily, -order_val, Class))

if (length(condition_levels) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = length(condition_levels) - 1, name = "Set1")
  )
}

if (length(top_class_names) <= 2) {
  facet_spec <- facet_grid(mode ~ Class, scales = "free_x")
} else {
  facet_spec <- facet_wrap(
    vars(Class, mode),
    ncol = min(length(top_class_names), 2L) * 2L,
    scales = "free_x"
  )
}

p <- ggplot(
  df_summary,
  aes(x = Subfamily, y = mean_fraction, fill = condition)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_errorbar(
    aes(
      ymin = mean_fraction - sd_fraction,
      ymax = mean_fraction + sd_fraction
    ),
    position = position_dodge(0.8),
    width = 0.25,
    na.rm = TRUE
  ) +
  geom_point(
    data = df_filtered,
    aes(x = Subfamily, y = fraction),
    position = position_dodge(0.8),
    colour = "black",
    size = 1.5,
    show.legend = FALSE
  ) +
  facet_spec +
  scale_x_reordered() +
  scale_fill_manual(values = colours, name = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of TE-mapped small RNAs") +
  theme_cowplot(10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

if (length(top_class_names) <= 2) {
  plot_width <- max(6, length(top_class_names) * n_top * 0.5)
  plot_height <- 4
} else {
  n_class_rows <- ceiling(length(top_class_names) / 2)
  plot_width <- max(6, min(length(top_class_names), 2) * n_top * 0.5)
  plot_height <- 4 * n_class_rows
}
ggsave(snakemake@output[["pdf"]], p, width = plot_width, height = plot_height)
