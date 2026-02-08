suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
})

read_cfg <- function(path) {
  if (!file.exists(path)) {
    stop("Missing config: ", path)
  }

  lines <- readLines(path, warn = FALSE)

  lines <- gsub("\\r", "", lines)
  lines <- trimws(lines)

  lines <- lines[lines != ""]
  lines <- lines[!grepl("^#", lines)]

  keys <- character()
  vals <- character()

  for (ln in lines) {
    kv <- strsplit(ln, "=", fixed = TRUE)[[1]]
    if (length(kv) < 2) {
      next
    }

    k <- trimws(kv[1])
    v <- trimws(paste(kv[-1], collapse = "="))

    v <- gsub('^"|"$', '', v)

    keys <- c(keys, k)
    vals <- c(vals, v)
  }

  as.list(setNames(vals, keys))
}

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- "config/config.env"

if (length(args) >= 1) {
  cfg_path <- args[1]
}

cfg <- read_cfg(cfg_path)

need <- c(
  "SALMON_DIR",
  "DESEQ_DIR",
  "TX2GENE_TSV",
  "CONTROL_LABEL",
  "CASE_LABEL",
  "DESEQ_PADJ"
)

miss <- need[!need %in% names(cfg)]
if (length(miss) > 0) {
  stop("Missing in config.env: ", paste(miss, collapse = ", "))
}

quant_dir <- file.path(cfg$SALMON_DIR, "quants")
out_dir <- cfg$DESEQ_DIR
tx2gene_file <- cfg$TX2GENE_TSV

meta1 <- file.path(out_dir, "sample_metadata.csv")
meta2 <- "metadata.csv"
meta_file <- if (file.exists(meta1)) meta1 else meta2

if (!file.exists(meta_file)) {
  stop("Missing metadata: ", meta1, " or ", meta2)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- read_csv(meta_file, show_col_types = FALSE)

if (!("sample_id" %in% names(meta)) && ("sample" %in% names(meta))) {
  meta <- meta %>% rename(sample_id = sample)
}

if (!("time_point" %in% names(meta)) && ("time" %in% names(meta))) {
  meta <- meta %>% rename(time_point = time)
}

if (!("replicate" %in% names(meta)) && ("rep" %in% names(meta))) {
  meta <- meta %>% rename(replicate = rep)
}

need_cols <- c("sample_id", "genotype", "time_point", "treatment")
miss <- need_cols[!need_cols %in% names(meta)]
if (length(miss) > 0) {
  stop("Missing columns in metadata: ", paste(miss, collapse = ", "))
}

meta <- meta %>%
  mutate(
    treatment = factor(
      treatment,
      levels = c(cfg$CONTROL_LABEL, cfg$CASE_LABEL)
    )
  )

tx2gene <- read_tsv(
  tx2gene_file,
  col_names = c("TXNAME", "GENEID"),
  show_col_types = FALSE
)

files <- file.path(quant_dir, meta$sample_id, "quant.sf")
names(files) <- meta$sample_id

missing <- names(files)[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing quant.sf for: ", paste(missing, collapse = ", "))
}

padj_cut <- as.numeric(cfg$DESEQ_PADJ)

run_one <- function(meta_sub, label) {
  meta_sub <- meta_sub %>% arrange(treatment, sample_id)
  f <- files[meta_sub$sample_id]

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

  res <- results(
    dds,
    contrast = c("treatment", cfg$CASE_LABEL, cfg$CONTROL_LABEL)
  )

  res_df <- as.data.frame(res)
  res_df$id <- rownames(res_df)
  res_df <- res_df %>% relocate(id)

  norm <- counts(dds, normalized = TRUE)
  norm_df <- as.data.frame(norm)
  norm_df$id <- rownames(norm_df)
  norm_df <- norm_df %>% relocate(id)

  lab_dir <- file.path(out_dir, label)
  dir.create(lab_dir, showWarnings = FALSE, recursive = TRUE)

  write_csv(res_df, file.path(lab_dir, paste0(label, "_de.csv")))
  write_csv(norm_df, file.path(lab_dir, paste0(label, "_norm.csv")))

  res_sig <- res_df %>% filter(!is.na(padj), padj < padj_cut)
  write_csv(
    res_sig,
    file.path(lab_dir, paste0(label, "_de_padj_", padj_cut, ".csv"))
  )

  write_csv(
    tibble::tibble(id = res_sig$id),
    file.path(lab_dir, paste0(label, "_tags_padj_", padj_cut, ".csv"))
  )

  write_csv(
    tibble::tibble(id = res_sig$id[res_sig$log2FoldChange > 0]),
    file.path(
      lab_dir,
      paste0(label, "_tags_up_padj_", padj_cut, ".csv")
    )
  )

  write_csv(
    tibble::tibble(id = res_sig$id[res_sig$log2FoldChange < 0]),
    file.path(
      lab_dir,
      paste0(label, "_tags_down_padj_", padj_cut, ".csv")
    )
  )

  png(file.path(lab_dir, paste0(label, "_ma.png")))
  plotMA(res, ylim = c(-4, 4))
  dev.off()

  vsd <- vst(dds, blind = FALSE)

  png(file.path(lab_dir, paste0(label, "_pca.png")))
  print(plotPCA(vsd, intgroup = "treatment"))
  dev.off()

  sel <- order(rowMeans(norm), decreasing = TRUE)
  sel <- sel[seq_len(min(30, length(sel)))]

  png(file.path(lab_dir, paste0(label, "_heatmap.png")))
  pheatmap(assay(vsd)[sel, ], show_rownames = FALSE)
  dev.off()

  tibble::tibble(
    comparison = label,
    n_samples = nrow(meta_sub),
    n_sig = nrow(res_sig)
  )
}

groups <- meta %>%
  distinct(genotype, time_point) %>%
  arrange(genotype, time_point)

summary_list <- list()

for (i in seq_len(nrow(groups))) {
  g <- as.character(groups$genotype[i])
  t <- as.character(groups$time_point[i])

  meta_sub <- meta %>% filter(genotype == g, time_point == t)

  if (length(unique(meta_sub$treatment)) < 2) {
    message("Skip ", g, "_", t, " (only one treatment)")
    next
  }

  label <- paste0(g, "_", t, "_", cfg$CASE_LABEL, "_vs_", cfg$CONTROL_LABEL)
  message("Running: ", label)

  summary_list[[length(summary_list) + 1]] <- run_one(meta_sub, label)
}

if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  write_csv(summary_df, file.path(out_dir, "deseq2_summary.csv"))
}

sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
# Placeholder wrapper.
# Use your original DESeq2 script under scripts/original/ if preferred.
