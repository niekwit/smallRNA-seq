log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(DESeq2)

tesmall_files <- snakemake@input[["txt"]]
reference_condition <- snakemake@config$reference_condition

# Build sample -> condition vector, reference condition first
samples <- setNames(
  sapply(tesmall_files, function(f) {
    snakemake@config$samples[[basename(dirname(f))]]
  }),
  sapply(tesmall_files, function(f) basename(dirname(f)))
)
ref_samples <- samples[samples == reference_condition]
samples <- c(ref_samples, samples[samples != reference_condition])

condition_factors <- factor(samples, levels = unique(samples))

# Build count matrix from all files
df <- read_delim(tesmall_files[1], show_col_types = FALSE)
for (file in tesmall_files[2:length(tesmall_files)]) {
  temp_df <- read_delim(file, show_col_types = FALSE)
  df <- full_join(df, temp_df, by = c("fid", "ftype"))
}
df[is.na(df)] <- 0
df$uid <- paste0(df$fid, "$", df$ftype)

count_matrix <- unique(df) %>%
  select(-ftype, -fid) %>%
  column_to_rownames(var = "uid") %>%
  as.matrix()
count_matrix <- count_matrix[, names(samples)]

# Run DESeq2 on all conditions
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = data.frame(condition = condition_factors),
  design = ~condition
)
dds <- DESeq(dds)

save(dds, file = snakemake@output[["rds"]])
