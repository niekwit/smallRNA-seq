# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(cowplot)
library(DESeq2)
library(ggrepel)

# Load dds
load(snakemake@input[["rds"]])

# Extracting transformed values
if (nrow(dds) < 1000) {
  print("Fewer than 1000 features. Using VST directly.")
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  # Try-catch to handle the 'nsub' error if data is sparse
  vsd <- tryCatch(
    {
      vst(dds, blind = FALSE)
    },
    error = function(e) {
      print(
        "vst failed due to low row counts, falling back to varianceStabilizingTransformation"
      )
      varianceStabilizingTransformation(dds, blind = FALSE)
    }
  )
}

# Set colours for plotting
vsd_df <- as.data.frame(colData(vsd))

if (length(unique(vsd_df$condition)) == 2) {
  colours <- c("#cccccc", "#dd3b3b")
} else {
  colours <- c(
    "#cccccc",
    RColorBrewer::brewer.pal(n = length(vsd_df$condition) - 1, name = "Set1")
  )
}

# Create PCA plot
p <- DESeq2::plotPCA(vsd, intgroup = "condition") +
  coord_cartesian() +
  labs(title = NULL) +
  theme_cowplot(12) +
  theme(legend.position = "none") +
  scale_color_manual(
    values = colours
  ) +
  geom_text_repel(aes(label = name))

# Save PCA plot to file
ggsave(
  filename = snakemake@output[["pdf"]],
  plot = p,
  width = 5,
  height = 2.5
)
