suppressWarnings(suppressMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Missing R package: ggplot2")
  }
}))

# ----------------------------
# USER PATHS (ALL EXPLICIT)
# ----------------------------
BASE_DIR <- "/QRISdata/Q9062"

ANN_PATH <- file.path(BASE_DIR, "06_ref", "gene_annotation_R570P14.tsv")
KO2PW_PATH <- file.path(BASE_DIR, "06_ref", "kegg", "ko2pathway.tsv")
PWLIST_PATH <- file.path(BASE_DIR, "06_ref", "kegg", "pathway_list.tsv")

DESEQ_ROOT <- file.path(BASE_DIR, "08_deseq2_R570")
OUT_ROOT <- file.path(BASE_DIR, "10_enrich")

CONTRASTS <- c(
  "Q208_12w_RKN_vs_C",
  "Q208_21d_RKN_vs_C",
  "Q208_7d_RKN_vs_C",
  "SES208_12w_RKN_vs_C"
)

SIG_PADJ <- 0.001
Q_CUT <- 0.05
MIN_GS <- 10
MAX_GS <- 5000
MIN_K <- 2

# ----------------------------
# HELPERS
# ----------------------------
split_ids <- function(x) {
  x <- gsub("[;|]", ",", x)
  x <- gsub("[[:space:]]+", ",", x)
  v <- unlist(strsplit(x, ",", fixed = TRUE))
  v <- trimws(v)
  v <- v[nzchar(v)]
  unique(v)
}

read_annotation_gene2ko <- function(path) {
  ann <- read.delim(path, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, quote = "")
  if (!("gene_id" %in% names(ann))) stop("ANN missing gene_id column")
  if (!("KO" %in% names(ann))) stop("ANN missing KO column")
  
  ann$KO <- trimws(ann$KO)
  ann$KO[ann$KO == ""] <- NA
  
  idx <- which(!is.na(ann$KO))
  gene <- character(0)
  ko <- character(0)
  
  for (i in idx) {
    ks <- split_ids(ann$KO[i])
    if (!length(ks)) next
    gene <- c(gene, rep(ann$gene_id[i], length(ks)))
    ko <- c(ko, ks)
  }
  
  data.frame(gene_id = gene, ko = ko, stringsAsFactors = FALSE)
}

read_ko2pathway <- function(path) {
  x <- read.delim(path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("ko2pathway.tsv must have >=2 columns")
  
  ko <- sub("^ko:", "", x[[1]])
  pw <- sub("^path:", "", x[[2]])
  
  data.frame(pathway = pw, ko = ko, stringsAsFactors = FALSE)
}


read_pathway_names <- function(path) {
  x <- read.delim(path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, quote = "")
  colnames(x) <- c("pathway", "name")
  x$pathway <- sub("^path:", "", x$pathway)
  x
}

read_deseq_results <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("id", "padj", "log2FoldChange")
  miss <- setdiff(need, names(x))
  if (length(miss)) stop("DESeq2 results missing: ", paste(miss, collapse = ", "))
  x
}

ora_hyper <- function(deg, uni, term2gene, term2name,
                      q_cut, min_gs, max_gs, min_k) {
  
  deg <- unique(deg)
  uni <- unique(uni)
  
  term2gene <- term2gene[term2gene$gene_id %in% uni, ]
  all_terms <- unique(term2gene$term_id)
  
  M <- length(uni)
  N <- length(deg)
  
  out <- list()
  
  for (t in all_terms) {
    genes_t <- unique(term2gene$gene_id[term2gene$term_id == t])
    m <- length(genes_t)
    if (m < min_gs || m > max_gs) next
    
    hit <- intersect(deg, genes_t)
    k <- length(hit)
    if (k < min_k) next
    
    p <- phyper(q = k - 1, m = m, n = M - m, k = N, lower.tail = FALSE)
    
    out[[t]] <- data.frame(
      term_id = t,
      term_name = ifelse(t %in% names(term2name), term2name[[t]], t),
      k = k,
      m = m,
      N = N,
      M = M,
      pvalue = p,
      padj = NA_real_,
      gene_ratio = k / N,
      bg_ratio = m / M,
      genes = paste(hit, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  if (!length(out)) return(NULL)
  
  res <- do.call(rbind, out)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$padj, -res$k, res$pvalue), ]
  res <- res[res$padj <= q_cut, ]
  if (!nrow(res)) return(NULL)
  res
}

plot_dot <- function(df, out_png, top_n, title) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  df <- df[seq_len(min(top_n, nrow(df))), ]
  df$term_name <- factor(df$term_name, levels = rev(df$term_name))
  df$mlog10 <- -log10(df$padj)
  
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = gene_ratio, y = term_name,
                                    size = k, color = mlog10)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = title, x = "Gene ratio (k/N)", y = NULL) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 160)
}

plot_bar <- function(df, out_png, top_n, title) {
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  df <- df[seq_len(min(top_n, nrow(df))), ]
  df$term_name <- factor(df$term_name, levels = rev(df$term_name))
  
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = term_name, y = -log10(padj))) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "-log10(FDR)") +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 160)
}

# ----------------------------
# MAIN
# ----------------------------
cat("Using paths:\n")
cat("  ANN:", ANN_PATH, "\n")
cat("  KO2PW:", KO2PW_PATH, "\n")
cat("  PWLIST:", PWLIST_PATH, "\n")
cat("  DESEQ_ROOT:", DESEQ_ROOT, "\n")
cat("  OUT_ROOT:", OUT_ROOT, "\n\n")

if (!file.exists(ANN_PATH)) stop("Missing: ", ANN_PATH)
if (!file.exists(KO2PW_PATH)) stop("Missing: ", KO2PW_PATH)
if (!file.exists(PWLIST_PATH)) stop("Missing: ", PWLIST_PATH)

gene2ko <- read_annotation_gene2ko(ANN_PATH)
# Filter for plant reasonable explanation
ko2pw <- read_ko2pathway(KO2PW_PATH)
ko2pw <- ko2pw[grepl("^map", ko2pw$pathway), ]
ko2pw <- ko2pw[!grepl("^map05", ko2pw$pathway), ]

pw_names <- read_pathway_names(PWLIST_PATH)
term2name <- setNames(pw_names$name, pw_names$pathway)

gene2pw <- merge(gene2ko, ko2pw, by = "ko")
gene2pw <- gene2pw[, c("gene_id", "pathway")]
colnames(gene2pw) <- c("gene_id", "term_id")

cat("Sanity check:\n")
cat("  gene2ko pairs:", nrow(gene2ko), "\n")
cat("  ko2pw pairs:", nrow(ko2pw), "\n")
cat("  gene2pw pairs:", nrow(gene2pw), "\n")
cat("  unique pathways:", length(unique(gene2pw$term_id)), "\n\n")

for (ct in CONTRASTS) {
  cat("==", ct, "==\n")
  
  res_path <- file.path(DESEQ_ROOT, ct, paste0(ct, "_DESeq2.DE.results.csv"))
  if (!file.exists(res_path)) {
    stop("Missing DESeq2 file: ", res_path)
  }
  
  out_dir <- file.path(OUT_ROOT, ct)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  res <- read_deseq_results(res_path)
  uni <- unique(res$id)
  
  deg_all <- res$id[!is.na(res$padj) & res$padj < SIG_PADJ]
  deg_up <- res$id[!is.na(res$padj) & res$padj < SIG_PADJ &
                     res$log2FoldChange > 0]
  deg_dn <- res$id[!is.na(res$padj) & res$padj < SIG_PADJ &
                     res$log2FoldChange < 0]
  
  run_one <- function(deg, tag) {
    df <- ora_hyper(
      deg = deg,
      uni = uni,
      term2gene = gene2pw,
      term2name = term2name,
      q_cut = Q_CUT,
      min_gs = MIN_GS,
      max_gs = MAX_GS,
      min_k = MIN_K
    )
    
    csv_out <- file.path(out_dir, paste0("KEGGpw_", tag, ".csv"))
    if (is.null(df) || !nrow(df)) {
      cat("  KEGGpw_", tag, ": 0 terms\n", sep = "")
      return(invisible(NULL))
    }
    
    write.csv(df, csv_out, row.names = FALSE)
    
    plot_dot(df,
             file.path(out_dir, paste0("KEGGpw_", tag, ".dotplot.png")),
             top_n = 20,
             title = paste(ct, "KEGG pathway", tag))
    
    plot_bar(df,
             file.path(out_dir, paste0("KEGGpw_", tag, ".barplot.png")),
             top_n = 20,
             title = paste(ct, "KEGG pathway", tag))
    
    cat("  KEGGpw_", tag, ": ", nrow(df), " terms\n", sep = "")
  }
  
  run_one(deg_all, "all")
  run_one(deg_up, "up")
  run_one(deg_dn, "down")
}

cat("\nDone.\n")
