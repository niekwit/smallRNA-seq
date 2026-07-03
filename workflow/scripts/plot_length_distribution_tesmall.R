log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

ref_condition <- snakemake@config[["reference_condition"]]

# TEsmall's *.anno.rlen.info classes (columns) that go into this plot, their
# display labels, and which facet ("TE" vs "miRNA") they belong to.
plot_classes <- c("sense_TE", "anti_TE", "miRNA")
class_labels <- c(sense_TE = "Sense", anti_TE = "Antisense", miRNA = "miRNA")
feature_group_map <- c(sense_TE = "TE", anti_TE = "TE", miRNA = "miRNA")

# Get all .anno.rlen.info files
rlen_files <- snakemake@input[["rlen_info"]]

# Load data from each file and add sample and condition information.
# Each file is a wide table: one row per read length, one column per TEsmall
# feature class (e.g. rlen, anti_TE, exon, ..., sense_TE), counts as values.
# Pivot to long format (rlen, class, count) to match the rest of the script.
rlen <- lapply(rlen_files, function(rlen_file) {
  sample_name <- str_replace(
    basename(rlen_file),
    "\\.anno\\.rlen\\.info$",
    ""
  )

  condition_name <- snakemake@config$samples[[sample_name]]
  read_delim(
    rlen_file,
    show_col_types = FALSE
  ) %>%
    pivot_longer(-rlen, names_to = "class", values_to = "count") %>%
    mutate(
      sample = sample_name,
      condition = condition_name
    )
})

# --- Calculate per-sample CPM (counts per million reads), then mean/SD ---
# CPM is normalised against each sample's full read total (all classes, all
# lengths), so it reflects the actual proportion of the library, not just of
# the TE-mapped reads plotted below.

df <- bind_rows(rlen)
conditions <- unique(df$condition)
condition_levels <- c(ref_condition, conditions[conditions != ref_condition])

df <- df %>%
  group_by(sample) %>%
  mutate(sample_sum = sum(count)) %>%
  ungroup() %>%
  mutate(cpm = count / sample_sum * 1e6) %>%
  group_by(condition, rlen, class) %>%
  mutate(condition_cpm = mean(cpm)) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels = condition_levels))

# Per-sample total CPM within each facet (anti_TE + sense_TE summed for the
# "TE" facet; miRNA alone for the "miRNA" facet), used for the overall bar
# height and its error bar.
df_total <- df %>%
  filter(class %in% plot_classes) %>%
  mutate(feature_group = feature_group_map[class]) %>%
  group_by(sample, condition, rlen, feature_group) %>%
  summarise(total_cpm = sum(cpm), .groups = "drop") %>%
  group_by(condition, rlen, feature_group) %>%
  mutate(
    condition_total_cpm = mean(total_cpm),
    total_sd = sd(total_cpm)
  ) %>%
  ungroup()

# --- Stacked bar graph (TE: anti_TE + sense_TE; miRNA: single class),
# faceted by feature group, with error bars ---
# Manually offset x-positions so bars can be dodged by condition while still
# letting ggplot stack sub-classes within each condition's slot.
# Conditions are laid out symmetrically around each rlen tick, splitting a
# fixed total width evenly so any number of conditions fits without overlap.
n_conditions <- length(condition_levels)
total_width <- 0.9
bar_width <- total_width / n_conditions
nudge_by_condition <- setNames(
  (seq_len(n_conditions) - (n_conditions + 1) / 2) * bar_width,
  condition_levels
)
nudge <- function(condition) nudge_by_condition[as.character(condition)]

df_summary <- df %>%
  filter(class %in% plot_classes) %>%
  select(condition, rlen, class, condition_cpm) %>%
  distinct() %>%
  mutate(
    feature_group = factor(feature_group_map[class], levels = c("TE", "miRNA")),
    x_pos = rlen + nudge(condition),
    fill_group = paste0(condition, " (", class_labels[class], ")")
  )

df_err <- df_total %>%
  select(condition, rlen, feature_group, condition_total_cpm, total_sd) %>%
  distinct() %>%
  mutate(
    feature_group = factor(feature_group, levels = c("TE", "miRNA")),
    x_pos = rlen + nudge(condition)
  )

# Class order within each condition's block follows plot_classes (Sense,
# Antisense, miRNA); condition blocks follow condition_levels (reference
# condition first).
class_order <- unname(class_labels[plot_classes])
fill_levels <- as.vector(outer(
  class_order,
  condition_levels,
  function(cl, cond) paste0(cond, " (", cl, ")")
))

df_summary$fill_group <- factor(df_summary$fill_group, levels = fill_levels)

# One base colour per condition (reference condition is grey, others follow
# Set1), with a darkened variant used for the "Antisense" class so each
# condition's bars stay visually grouped in the legend.
if (n_conditions == 2) {
  base_colours <- c("#cccccc", "#dd3b3b")
} else {
  base_colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = n_conditions - 1, name = "Set1")
  )
}
names(base_colours) <- condition_levels

darken <- function(hex, factor = 0.55) {
  rgb_mat <- grDevices::col2rgb(hex) * factor
  grDevices::rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ], maxColorValue = 255)
}

fill_colours <- setNames(
  as.vector(outer(
    class_order,
    condition_levels,
    function(cl, cond) {
      base <- base_colours[cond]
      ifelse(cl == "Antisense", darken(base), base)
    }
  )),
  fill_levels
)

p <- ggplot() +
  geom_col(
    data = df_summary,
    aes(x = x_pos, y = condition_cpm, fill = fill_group, group = x_pos),
    position = position_stack(reverse = TRUE),
    width = bar_width,
    colour = "black"
  ) +
  geom_errorbar(
    data = df_err,
    aes(
      x = x_pos,
      ymin = condition_total_cpm - total_sd,
      ymax = condition_total_cpm + total_sd
    ),
    width = 0.15
  ) +
  scale_fill_manual(
    values = fill_colours,
    name = NULL
  ) +
  guides(fill = guide_legend(nrow = length(class_order))) +
  facet_wrap(~feature_group, nrow = 1, scales = "free_y") +
  scale_x_continuous(
    breaks = sort(unique(df_summary$rlen)),
    labels = sort(unique(df_summary$rlen))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Length (nt)", y = "Reads per million (CPM)") +
  theme_cowplot(14) +
  theme(legend.position = "bottom")


ggsave(
  file.path(snakemake@output[["pdf"]]),
  p,
  width = 11,
  height = 6
)

# Save data as CSV
write_csv(
  df,
  file = snakemake@output[["csv"]]
)
