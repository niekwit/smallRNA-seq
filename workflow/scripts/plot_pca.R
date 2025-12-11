# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(DESeq2)

# Load dds
load(snakemake@input[["rds"]])

# Transform counts for PCA
vsd <- vst(dds, blind = FALSE)

# Set colours for plotting
vsd_df <- as.data.frame(colData(vsd))

if (length(unique(vsd_df$condition)) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- RColorBrewer::brewer.pal(
    n = length(unique(vsd_df$condition)),
    name = "Set3"
  )
}

# Create PCA plot
p <- DESeq2::plotPCA(vsd, intgroup = "condition") +
  coord_fixed(ratio = NULL) +
  labs(
    title = NULL,
    colour = "Genotype"
  ) +
  theme_cowplot(12) +
  scale_color_manual(
    values = c(
      "Morc2a_WT" = "#cccccc",
      "Morc2a_KO" = "#dd3b3b"
    )
  )

# Save PCA plot to file
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 5,
  height = 2.5
)
