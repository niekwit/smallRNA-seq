# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(DESeq2)


# Prepare data for DESeq2 analysis
# -----------------------------------------------

# Get reference and test condition names
comparison <- snakemake@wildcards[["comparison"]]
reference_condition <- str_split(comparison, "_vs_")[[1]][2]
test_condition <- str_split(comparison, "_vs_")[[1]][1]

# Get input files
tesmall_files <- snakemake@input[["txt"]]

# Create named vector with sample names and conditions
samples <- vector()
for (file in tesmall_files) {
  sample <- basename(dirname(file))
  condition <- snakemake@config$samples[[sample]]
  samples[sample] <- condition
}

# Reorder vector to have reference condition first
reference_condition <- snakemake@config$reference_condition
ref_samples <- samples[names(samples)[samples == reference_condition]]
samples <- c(
  ref_samples,
  samples[names(samples)[samples != reference_condition]]
)

# Create vector of factors using conditions
condition_factors <- factor(samples, levels = unique(samples))

# Open first file to establish data frame structure
df <- read_delim(tesmall_files[1])

# Loop through remaining files and combine into single data frame
for (file in tesmall_files[2:length(tesmall_files)]) {
  temp_df <- read_delim(file)
  df <- full_join(df, temp_df, by = c("fid", "ftype"))
}

# Replace NAs with zeros
df[is.na(df)] <- 0

# Create unique identifier for each feature
df$uid <- paste0(df$fid, "$", df$ftype)

# Collapse anti_TE and sense_TE ftype into single TE ftype
df$ftype <- str_replace_all(
  df$ftype,
  c("anti_TE" = "TE", "sense_TE" = "TE")
)

# Create matrix of counts with fid as rownames
# no ftype column (extract from df later)
count_matrix <- unique(df) %>%
  select(-ftype, -fid) %>%
  column_to_rownames(var = "uid") %>%
  as.matrix()

# Reorder columns of count matrix to match samples vector
count_matrix <- count_matrix[, names(samples)]


# Run DESeq2 analysis
# -----------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(
    condition = condition_factors
  ),
  design = ~condition
)
dds <- DESeq(dds)
res <- results(dds)

# Save DESeq2 object for downstream analysis
save(dds, file = snakemake@output[["rds"]])

# Remove conditions from dds that are not in the current comparison
# (in case multiple comparisons are being run)
keep_samples <- rownames(colData(dds))[
  colData(dds)$condition %in% c(reference_condition, test_condition)
]
dds <- dds[, keep_samples]
dds$condition <- droplevels(dds$condition)

# Get log2 fold changes and adjusted p-values
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "uid") %>%
  select(-c(baseMean, lfcSE, stat, pvalue))

# Add normalised counts
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column(var = "uid")

res_df <- left_join(res_df, norm_counts, by = "uid") %>%
  # Separate uid into fid and ftype
  separate(uid, into = c("fid", "ftype"), sep = "\\$") %>%
  arrange(padj)


# Annotate features with chromosome, start, end, strand
# -----------------------------------------------

# First check if annotation files exists
# This should be there as TEsmall creates/uses it during its run
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

  # For plotting purposes, set final_df to res_df
  final_df <- res_df
} else {
  annotation_list <- list()

  for (i in seq_along(annotation_files)) {
    temp_annotation <- read_delim(
      annotation_files[i],
      col_names = FALSE
    )
    names(temp_annotation) <- c(
      "chr",
      "start",
      "end",
      "fid",
      "score",
      "strand"
    )
    temp_annotation <- temp_annotation %>%
      select(fid, chr, start, end, strand)
    annotation_list[[i]] <- temp_annotation
  }
  # Combine all annotations into single data frame
  # and remove duplicates (some elements may be present in multiple files)
  annotation_df <- bind_rows(annotation_list) %>%
    distinct(fid, .keep_all = TRUE)

  # Merge results with annotation data
  final_df <- left_join(res_df, annotation_df, by = "fid")

  # Save results with annotation
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

# Save plot to file
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)
