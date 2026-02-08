params <- list(
  base_dir = "/QRISdata/Q9062/08_deseq2_R570",
  ann_tsv = "/QRISdata/Q9062/06_ref/gene_annotation_R570P14.tsv",
  out_dir = "/QRISdata/Q9062/10_enrich",
  alpha = 0.001,
  min_gs_size = 5,
  max_gs_size = 5000,
  q_cut = 0.05
)

# BiocManager::install("clusterProfiler")
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
})

read_ann <- function(path) {
  ann <- read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
    mutate(
      GO = ifelse(is.na(GO), "", GO),
      KO = ifelse(is.na(KO), "", KO)
    )
  ann
}

make_term2gene_go <- function(ann) {
  go <- ann %>%
    select(gene_id, GO) %>%
    filter(GO != "") %>%
    separate_rows(GO, sep = "\s+") %>%
    filter(str_detect(GO, "^GO:\d+")) %>%
    distinct(GO, gene_id) %>%
    rename(term = GO, gene = gene_id)
  go
}

make_term2gene_ko <- function(ann) {
  ko <- ann %>%
    select(gene_id, KO) %>%
    filter(KO != "") %>%
    separate_rows(KO, sep = "\s+") %>%
    filter(str_detect(KO, "^K\d+")) %>%
    distinct(KO, gene_id) %>%
    rename(term = KO, gene = gene_id)
  ko
}

run_enricher <- function(genes, universe, term2gene, q_cut,
                         min_gs, max_gs) {
  genes <- unique(genes)
  universe <- unique(universe)
  
  if (length(genes) < min_gs) {
    return(NULL)
  }
  
  eg <- enricher(
    gene = genes,
    universe = universe,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    qvalueCutoff = q_cut,
    minGSSize = min_gs,
    maxGSSize = max_gs
  )
  eg
}

save_enrich <- function(eg, out_prefix) {
  if (is.null(eg)) {
    return(invisible(NULL))
  }
  df <- as.data.frame(eg)
  write.csv(df, paste0(out_prefix, ".csv"), row.names = FALSE)
  
  p1 <- barplot(eg, showCategory = 20)
  ggsave(paste0(out_prefix, ".barplot.png"), p1,
         width = 10, height = 6, dpi = 200)
  
  p2 <- dotplot(eg, showCategory = 20)
  ggsave(paste0(out_prefix, ".dotplot.png"), p2,
         width = 10, height = 6, dpi = 200)
  
  invisible(df)
}
