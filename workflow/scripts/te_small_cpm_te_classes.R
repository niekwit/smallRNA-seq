# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

### Plot CPM of TE classes (sense/antisense), per condition ###
# Same input/classification as te_small_stacked_bar_te.R, but reporting
# CPM (per million total reads in the sample) instead of percentage of
# TE-mapped reads, with conditions plotted side by side per Class.
count_files <- snakemake@input[["txt"]]

# Empty list to store results
df_list <- list()

for (i in seq_along(count_files)) {
  sample_name <- basename(dirname(count_files[i]))
  condition_name <- snakemake@config$samples[[sample_name]]
  # Open count summary file
  raw <- read_delim(
    count_files[i],
    show_col_types = FALSE,
    col_names = c(
      "fid",
      "ftype",
      "counts"
    )
  ) %>%
    mutate(counts = as.numeric(counts)) %>%
    filter(!is.na(counts))

  # CPM denominator is the sample's full library (all ftypes), not just
  # the TE-mapped reads plotted below.
  total_reads <- sum(raw$counts, na.rm = TRUE)

  counts <- raw %>%
    filter(ftype %in% c("sense_TE", "anti_TE")) %>%
    separate(
      fid,
      into = c("Class", "Family", "Subfamily", "Instance"),
      sep = ":"
    ) %>%
    mutate(
      Class = case_when(
        # Split LTR class into LTR and IAP
        Class == "LTR" & grepl("IAP", Subfamily) ~ "IAP",
        TRUE ~ Class
      ),
      # Convert minor classes to Other
      Class = case_when(
        Class %in%
          c(
            "DNA",
            "SINE",
            "RC",
            "Retroposon",
            "Unspecified",
            "Unknown",
            "Satellite"
          ) ~ "Other",
        TRUE ~ Class
      )
    ) %>%
    group_by(Class, ftype) %>%
    summarise(
      Class_Count = sum(counts, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cpm = Class_Count / total_reads * 1e6,
      sample = sample_name,
      condition = condition_name
    )

  df_list[[i]] <- counts
}

# Combine all samples into single data frame
df <- bind_rows(df_list)

# Calculate mean/SD CPM per condition, Class and strand across replicates
df_summary <- df %>%
  group_by(condition, Class, ftype) %>%
  summarise(
    mean_cpm = mean(cpm),
    sd_cpm = sd(cpm),
    .groups = "drop"
  )

# Save to file
write_csv(df_summary, snakemake@output[["csv"]])

# Convert condition to factor with desired order (reference condition first)
reference_condition <- snakemake@config$reference_condition
conditions <- unique(df$condition)
condition_levels <- c(
  reference_condition,
  setdiff(conditions, reference_condition)
)
df$condition <- factor(df$condition, levels = condition_levels)
df_summary$condition <- factor(df_summary$condition, levels = condition_levels)

# Set colour palette (reference condition grey, others from Set1)
if (length(condition_levels) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = length(condition_levels) - 1, name = "Set1")
  )
}

# Convert Class to factor with desired order
class_levels <- c("LINE", "LTR", "IAP", "Other")
df$Class <- factor(df$Class, levels = class_levels)
df_summary$Class <- factor(df_summary$Class, levels = class_levels)

# Convert ftype to factor with desired order and readable facet labels
ftype_levels <- c("anti_TE", "sense_TE")
ftype_labels <- c("Antisense", "Sense")
df$ftype <- factor(df$ftype, levels = ftype_levels, labels = ftype_labels)
df_summary$ftype <- factor(
  df_summary$ftype,
  levels = ftype_levels,
  labels = ftype_labels
)

# Create grouped bar plot (conditions dodged per Class), faceted by strand
p <- ggplot(
  df_summary,
  aes(x = Class, y = mean_cpm, fill = condition)
) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_errorbar(
    aes(ymin = mean_cpm - sd_cpm, ymax = mean_cpm + sd_cpm),
    position = position_dodge(0.8),
    width = 0.25
  ) +
  geom_point(
    data = df,
    aes(x = Class, y = cpm),
    position = position_dodge(0.8),
    colour = "black",
    size = 1.5,
    show.legend = FALSE
  ) +
  facet_wrap(~ftype) +
  labs(
    x = NULL,
    y = "TE-derived piRNAs (CPM)",
    fill = NULL
  ) +
  theme_cowplot(14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# Save plot to file
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 6,
  height = 6
)
