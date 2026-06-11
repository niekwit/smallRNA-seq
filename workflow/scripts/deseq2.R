log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)
library(DESeq2)

# Get reference and test condition names from wildcard
comparison <- snakemake@wildcards[["comparison"]]
reference_condition <- snakemake@config$reference_condition
test_condition <- str_split(comparison, "_vs_")[[1]][1]

# Load full multi-condition dds built by deseq2_full rule
load(snakemake@input[["dds_full"]])

# Get results for this specific pairwise contrast
res <- results(
  dds,
  contrast = c("condition", test_condition, reference_condition)
)

# Subset dds to the two comparison conditions for normalised counts
keep_samples <- rownames(colData(dds))[
  colData(dds)$condition %in% c(reference_condition, test_condition)
]
dds <- dds[, keep_samples]
dds$condition <- droplevels(dds$condition)

# Save per-comparison dds
save(dds, file = snakemake@output[["dds"]])

# Get log2 fold changes and adjusted p-values
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "uid") %>%
  select(-c(baseMean, lfcSE, stat, pvalue))

# Add normalised counts from the subsetted dds
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "uid")

res_df <- left_join(res_df, norm_counts, by = "uid") %>%
  separate(uid, into = c("fid", "ftype"), sep = "\\$") %>%
  arrange(padj)


# Annotate features with chromosome, start, end, strand
# -----------------------------------------------
genome <- snakemake@config$genome
dbfolder <- snakemake@config$tesmall$dbfolder
annotation_files <- Sys.glob(file.path(
  dbfolder,
  "genomes",
  genome,
  "annotation",
  "*.bed"
))
if (length(annotation_files) != 8) {
  print("Warning: No annotation files found for feature annotation!")
  print("Saving results without annotation.")
  write_csv(res_df, snakemake@output[["csv"]])
  final_df <- res_df
} else {
  annotation_list <- list()
  for (i in seq_along(annotation_files)) {
    temp_annotation <- read_delim(
      annotation_files[i],
      col_names = FALSE,
      show_col_types = FALSE
    )
    names(temp_annotation) <- c("chr", "start", "end", "fid", "score", "strand")
    annotation_list[[i]] <- temp_annotation %>%
      select(fid, chr, start, end, strand)
  }
  annotation_df <- bind_rows(annotation_list) %>%
    distinct(fid, .keep_all = TRUE)

  final_df <- left_join(res_df, annotation_df, by = "fid")
  write_csv(final_df, snakemake@output[["csv"]])
}

# Create volcano plot
# -----------------------------------------------
p <- ggplot(
  final_df,
  aes(x = log2FoldChange, y = -log10(padj), fill = ftype)
) +
  geom_point(
    aes(color = ftype),
    alpha = 0.6,
    size = 3,
    shape = 21,
    colour = "black"
  ) +
  theme_cowplot(16) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_fill_brewer(palette = "Set1")

ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)
