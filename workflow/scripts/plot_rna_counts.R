# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

rrna_files <- snakemake@input[["rrna"]]
mirna_files <- snakemake@input[["mirna"]]

# Sum counts encoded in collapsed FASTA headers (">readID-count")
fasta_total_count <- function(fasta_file) {
  headers <- grep("^>", readLines(fasta_file), value = TRUE)
  if (length(headers) == 0) {
    return(0L)
  }
  counts <- as.integer(str_split_fixed(str_remove(headers, "^>"), "-", 2)[, 2])
  sum(counts, na.rm = TRUE)
}

df <- bind_rows(
  data.frame(
    sample = str_remove(basename(rrna_files), "\\.fasta$"),
    count = sapply(rrna_files, fasta_total_count),
    type = "rRNA"
  ),
  data.frame(
    sample = str_remove(basename(mirna_files), "\\.fasta$"),
    count = sapply(mirna_files, fasta_total_count),
    type = "miRNA"
  )
)

# Order samples and types consistently
all_samples <- snakemake@config$samples
reference_condition <- snakemake@config$reference_condition
sample_to_condition <- setNames(unlist(all_samples), names(all_samples))

df <- df %>%
  mutate(
    condition = sample_to_condition[sample],
    type = factor(type, levels = c("rRNA", "miRNA"))
  )

# Reference condition samples first, then others
ref_samples <- names(all_samples)[all_samples == reference_condition]
other_samples <- names(all_samples)[all_samples != reference_condition]
df$sample <- factor(df$sample, levels = c(ref_samples, other_samples))

p <- ggplot(df, aes(x = sample, y = count, fill = type)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::scientific
  ) +
  scale_fill_manual(values = c(rRNA = "#4e79a7", miRNA = "#f28e2b")) +
  labs(
    x = NULL,
    y = "Read count",
    fill = NULL,
    title = "rRNA/miRNA read counts removed with Bowtie"
  ) +
  theme_cowplot() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = max(4, length(unique(df$sample)) * 0.8),
  height = 5
)
