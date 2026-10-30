log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

ref_condition <- snakemake@config[["reference_condition"]]

# TEsmall's *.anno.rlen.info classes (columns) that go into this plot and
# which panel ("TE" vs "miRNA") they belong to.
plot_classes <- c("sense_TE", "anti_TE", "miRNA")
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

# --- Bar graphs: TE (sense up, antisense down) and miRNA, with error bars
# and individual sample values ---
# Manually offset x-positions so bars can be dodged by condition.
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

# Reference condition is grey, others follow Set1
if (n_conditions == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = n_conditions - 1, name = "Set1")
  )
}
names(colours) <- condition_levels

# TE panel: antisense CPM is plotted as negative values to show the
# sense/antisense split more clearly
te_classes <- c("sense_TE", "anti_TE")

df_te_points <- df %>%
  filter(class %in% te_classes) %>%
  mutate(
    cpm_signed = cpm * if_else(class == "anti_TE", -1, 1),
    x_pos = rlen + nudge(condition)
  )

df_te <- df_te_points %>%
  group_by(condition, rlen, class, x_pos) %>%
  summarise(
    condition_cpm_signed = mean(cpm_signed),
    sd = sd(cpm_signed),
    .groups = "drop"
  )

# miRNA panel
df_mirna_points <- df_total %>%
  filter(feature_group == "miRNA") %>%
  mutate(x_pos = rlen + nudge(condition))

df_mirna <- df_mirna_points %>%
  select(
    feature_group,
    condition,
    rlen,
    x_pos,
    condition_total_cpm,
    total_sd
  ) %>%
  distinct()

# Shared styling for both panels
rlen_breaks <- sort(unique(df$rlen))
point_layer <- function(data, y) {
  geom_point(
    data = data,
    aes(x = x_pos, y = .data[[y]]),
    position = position_jitter(width = 0.05, height = 0, seed = 1),
    size = 1.5,
    shape = 16,
    colour = "black",
    show.legend = FALSE
  )
}
common_layers <- list(
  scale_fill_manual(values = colours, name = NULL),
  scale_x_continuous(breaks = rlen_breaks, labels = rlen_breaks),
  theme_cowplot(14),
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
)

p_te <- ggplot(
  df_te,
  aes(x = x_pos, y = condition_cpm_signed, fill = condition)
) +
  geom_col(width = bar_width, colour = "black") +
  geom_errorbar(
    aes(
      ymin = condition_cpm_signed - sd,
      ymax = condition_cpm_signed + sd
    ),
    width = 0.15
  ) +
  point_layer(df_te_points, "cpm_signed") +
  geom_hline(yintercept = 0, colour = "black") +
  facet_wrap(~"TE") +
  common_layers +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    x = "Length (nt)",
    y = "Reads per million (CPM)\nAntisense (-) / Sense (+)"
  )

p_mirna <- ggplot(
  df_mirna,
  aes(x = x_pos, y = condition_total_cpm, fill = condition)
) +
  geom_col(width = bar_width, colour = "black") +
  geom_errorbar(
    aes(
      ymin = condition_total_cpm - total_sd,
      ymax = condition_total_cpm + total_sd
    ),
    width = 0.15
  ) +
  point_layer(df_mirna_points, "total_cpm") +
  facet_wrap(~feature_group) +
  common_layers +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Length (nt)", y = "Reads per million (CPM)")

# Combine panels with one shared legend
legend <- get_plot_component(p_te, "guide-box-bottom")
p <- plot_grid(
  plot_grid(
    p_te + theme(legend.position = "none"),
    p_mirna + theme(legend.position = "none"),
    nrow = 1,
    align = "h",
    axis = "tb"
  ),
  legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

ggsave(
  file.path(snakemake@output[["pdf"]]),
  p,
  width = 13.2,
  height = 6
)

# Save data as CSV
write_csv(
  df,
  file = snakemake@output[["csv"]]
)
