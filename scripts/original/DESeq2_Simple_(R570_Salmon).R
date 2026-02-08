
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

proj_dir <- "/QRISdata/Q9062"
quant_dir <- file.path(proj_dir, "07_salmon_R570", "quants")
meta_file <- file.path(proj_dir, "08_deseq2", "sample_metadata.csv")
tx2gene_file <- file.path(proj_dir, "06_ref", "tx2gene_R570_1os2g.tsv")
out_dir <- file.path(proj_dir, "08_deseq2_R570")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load metadata and tx2gene
meta <- read_csv(meta_file, show_col_types = FALSE)

tx2gene <- read_tsv(
  tx2gene_file,
  col_names = c("TXNAME", "GENEID"),
  show_col_types = FALSE
)

files <- file.path(quant_dir, meta$sample, "quant.sf")
names(files) <- meta$sample

missing <- names(files)[!file.exists(files)]
if (length(missing) > 0) {
  stop(paste("Missing quant.sf for:", paste(missing, collapse = ", ")))
}

# Run one contrast (RKN vs C) within genotype + time
run_one <- function(meta_sub, label) {
meta_sub <- meta_sub %>% arrange(treatment, sample)
f <- files[meta_sub$sample]

txi <- tximport(
  f,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = meta_sub,
  design = ~ treatment
)

dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("treatment", "RKN", "C"))
res_df <- as.data.frame(res)
res_df$id <- rownames(res_df)
res_df <- res_df %>% relocate(id)

norm <- counts(dds, normalized = TRUE)
norm_df <- as.data.frame(norm)
norm_df$id <- rownames(norm_df)
norm_df <- norm_df %>% relocate(id)

lab_dir <- file.path(out_dir, label)
dir.create(lab_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(res_df, file.path(lab_dir, paste0(label, "_DESeq2.DE.results.csv")))
write_csv(norm_df, file.path(lab_dir, paste0(label, "_DESeq2.norm_counts.csv")))

sig <- 0.001
res_sig <- res_df %>% filter(!is.na(padj), padj < sig)

write_csv(
  res_sig,
  file.path(lab_dir,
            paste0(label, "*DESeq2.DE.results_padjust*", sig, ".csv"))
)

write_csv(
  tibble::tibble(id = res_sig$id),
  file.path(lab_dir,
            paste0(label, "*tags_DESeq2.DE.results_padjust*", sig, ".csv"))
)

write_csv(
  tibble::tibble(id = res_sig$id[res_sig$log2FoldChange > 0]),
  file.path(lab_dir,
            paste0(label, "*tags_DESeq2.DE.results_up.reg_padjust*", sig, ".csv"))
)

write_csv(
  tibble::tibble(id = res_sig$id[res_sig$log2FoldChange < 0]),
  file.path(lab_dir,
            paste0(label, "*tags_DESeq2.DE.results_down.reg_padjust*", sig, ".csv"))
)

png(file.path(lab_dir, paste0(label, "_DESeq2_MA_plot.png")))
plotMA(res, ylim = c(-4, 4))
dev.off()

vsd <- vst(dds, blind = FALSE)

png(file.path(lab_dir, paste0(label, "_DESeq2_PCA_plot.png")))
print(plotPCA(vsd, intgroup = "treatment"))
dev.off()

sel <- order(rowMeans(norm), decreasing = TRUE)
sel <- sel[seq_len(min(30, length(sel)))]

png(file.path(lab_dir, paste0(label, "_DESeq2_HeatMap.png")))
pheatmap(assay(vsd)[sel, ], show_rownames = FALSE)
dev.off()
}
# Run all genotype x time comparisons
groups <- meta %>% distinct(genotype, time) %>% arrange(genotype, time)

for (i in seq_len(nrow(groups))) {
  g <- as.character(groups$genotype[i])
  t <- as.character(groups$time[i])
  
  meta_sub <- meta %>% filter(genotype == g, time == t)
  
  if (length(unique(meta_sub$treatment)) < 2) {
    message("Skip ", g, "_", t, " (only one treatment)")
    next
  }
  
  label <- paste0(g, "_", t, "_RKN_vs_C")
  message("Running: ", label)
  
  run_one(meta_sub, label)
}



