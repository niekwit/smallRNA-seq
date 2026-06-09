# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Rsamtools)
library(dplyr)
library(ggplot2)
library(cowplot)

# Inputs from Snakemake
bam_input <- snakemake@input[["bam"]]
te <- snakemake@wildcards[["te"]]
out_pdf <- snakemake@output[["pdf"]]

# Build named list: condition -> character vector of BAM paths
all_samples <- snakemake@config$samples # named list: sample -> condition
conditions <- unique(unlist(all_samples))

bam_files <- lapply(conditions, function(cond) {
  samples_in_cond <- names(all_samples)[all_samples == cond]
  bam_input[
    sapply(samples_in_cond, function(s) grepl(s, bam_input)) |>
      rowSums() >
      0
  ]
})
names(bam_files) <- conditions

# Helper: read relevant BAM fields for a single file ---
read_bam <- function(bam_path) {
  param <- ScanBamParam(what = c("qname", "flag", "rname", "pos", "qwidth"))
  bam <- scanBam(bam_path, param = param)[[1]]
  data.frame(
    qname = bam$qname,
    flag = bam$flag,
    rname = as.character(bam$rname),
    pos = bam$pos,
    qwidth = bam$qwidth,
    stringsAsFactors = FALSE
  )
}

parse_count <- function(qnames) {
  as.integer(sub(".*-(\\d+)$", "\\1", qnames))
}

# Get TE sequence length from BAM header
header <- scanBamHeader(bam_input[[1]])
seq_len <- header[[1]]$targets[[te]]
if (is.null(seq_len)) stop("No BAM header entry found for TE: ", te)

# Load and merge replicates; keep only FLAG 0 (fwd) and FLAG 16 (rev)
# Filter to the target TE immediately to reduce memory usage
print("Loading BAM files...")
reads <- lapply(names(bam_files), function(cond) {
  df <- do.call(rbind, lapply(bam_files[[cond]], read_bam))
  df$condition <- cond
  df
})
reads <- do.call(rbind, reads)
reads <- reads[reads$flag %in% c(0L, 16L) & reads$rname == te, ]

reads$count <- parse_count(reads$qname)
reads <- reads[!is.na(reads$count) & !is.na(reads$pos), ]
reads$strand <- ifelse(reads$flag == 0L, "forward", "reverse")

# For reverse-strand reads, BAM pos is the leftmost (3') coordinate.
# The 5' end is pos + read_length - 1.
reads$five_p <- ifelse(
  reads$flag == 16L,
  reads$pos + reads$qwidth - 1,
  reads$pos
)

agg <- reads %>%
  group_by(condition, strand, five_p) %>%
  summarise(raw = sum(count), .groups = "drop") %>%
  group_by(condition, strand) %>%
  mutate(norm = raw / sum(raw) * 100) %>%
  ungroup() %>%
  mutate(
    norm = ifelse(strand == "reverse", -norm, norm),
    condition = factor(condition, levels = names(bam_files))
  )

p <- ggplot(
  agg,
  aes(x = five_p, xend = five_p, y = 0, yend = norm, colour = strand)
) +
  geom_segment(linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey60") +
  facet_wrap(~condition, ncol = 1) +
  scale_x_continuous(
    limits = c(0, seq_len),
    expand = c(0, 0),
    breaks = seq(0, seq_len, by = 1000),
  ) +
  scale_colour_manual(values = c(forward = "#2196F3", reverse = "#E53935")) +
  labs(
    title = te,
    x = "Position (bp)",
    y = "Normalised read count (x100)",
    colour = "Strand"
  ) +
  theme_cowplot(12) +
  theme(legend.position = "none")

ggsave(out_pdf, p, width = 10, height = 5)
print("Done.")
