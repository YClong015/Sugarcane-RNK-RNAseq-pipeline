suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(UpSetR)
})

# ---------------------------
# Paths (Hardcoded to your current directories)
# ---------------------------
BASE_DEG    <- "/QRISdata/Q9062/08_deseq2_R570"
BASE_ENR    <- "/QRISdata/Q9062/10_enrich"
META_CSV    <- "/QRISdata/Q9062/08_deseq2_R570/sample_metadata.csv"
GENESET_RGA <- "/QRISdata/Q9062/10_enrich/_gene_sets/rga_genes.txt"
GENESET_PR1 <- "/QRISdata/Q9062/10_enrich/_gene_sets/pr1_genes.txt"

OUTDIR <- "/QRISdata/Q9062/11_reports/figs"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Your 4 contrasts (based on your directory structure)
CONTRASTS <- c(
  "Q208_7d_RKN_vs_C",
  "Q208_12w_RKN_vs_C",
  "Q208_21d_RKN_vs_C",
  "SES208_12w_RKN_vs_C"
)

# ---------------------------
# Helper functions
# ---------------------------
read_deg_tables <- function(contrast) {
  ddir <- file.path(BASE_DEG, contrast)
  # Look for the DESeq2 results and normalized counts files
  res_file <- Sys.glob(file.path(ddir, "*DESeq2.DE.results.csv"))
  cnt_file <- Sys.glob(file.path(ddir, "*DESeq2.norm_counts.csv"))

  if (length(res_file) != 1) stop(paste("Cannot find unique results.csv:", contrast))
  if (length(cnt_file) != 1) stop(paste("Cannot find unique norm_counts.csv:", contrast))

  res <- read_csv(res_file, show_col_types = FALSE)
  cnt <- read_csv(cnt_file, show_col_types = FALSE)

  # Validate column names
  if (!"id" %in% names(res)) stop("results missing id")
  if (!"log2FoldChange" %in% names(res)) stop("results missing log2FoldChange")
  if (!"padj" %in% names(res)) stop("results missing padj")

  list(res = res, cnt = cnt, res_file = res_file, cnt_file = cnt_file)
}

get_deg_sets <- function(res, padj_cut = 0.001, lfc_cut = 1) {
  # Add direction column (up/down/ns)
  res2 <- res %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(dir = case_when(
      padj < padj_cut & log2FoldChange >=  lfc_cut ~ "up",
      padj < padj_cut & log2FoldChange <= -lfc_cut ~ "down",
      TRUE ~ "ns"
    ))

  list(
    all  = res2 %>% filter(padj < padj_cut, abs(log2FoldChange) >= lfc_cut) %>% pull(id),
    up   = res2 %>% filter(dir == "up") %>% pull(id),
    down = res2 %>% filter(dir == "down") %>% pull(id),
    res2 = res2
  )
}

save_png <- function(path, w=8, h=6, dpi=300) {
  ggsave(filename = path, width = w, height = h, dpi = dpi)
}

plot_volcano <- function(res2, title, out_png) {
  p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.35, size = 0.9) +
    labs(title = title, x = "log2FoldChange", y = "-log10(adjusted p)") +
    theme_bw()
  ggsave(out_png, p, width = 7.5, height = 6, dpi = 300)
}

plot_deg_bar <- function(deg_sets, title, out_png) {
  df <- data.frame(
    group = c("up", "down"),
    n = c(length(deg_sets$up), length(deg_sets$down))
  )
  p <- ggplot(df, aes(x = group, y = n)) +
    geom_col() +
    labs(title = title, x = "", y = "DEG count (padj<0.001, |LFC|>=1)") +
    theme_bw()
  ggsave(out_png, p, width = 5.5, height = 4.5, dpi = 300)
}

plot_sample_corr_heatmap <- function(cnt, title, out_png) {
  mat <- as.data.frame(cnt)
  rownames(mat) <- mat$id
  mat$id <- NULL
  mat <- as.matrix(mat)

  # log transform to stabilize variance
  mat <- log2(mat + 1)

  # sample-sample correlation (Pearson)
  cmat <- cor(mat, method = "pearson")

  png(out_png, width = 2000, height = 1700, res = 250)
  pheatmap(cmat, main = title, border_color = NA)
  dev.off()
}

plot_top_gene_heatmap <- function(cnt, res, top_n = 50, title, out_png) {
  mat <- as.data.frame(cnt)
  rownames(mat) <- mat$id
  mat$id <- NULL
  mat <- as.matrix(mat)

  # Select top significant genes by padj
  res2 <- res %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(top_n)
  ids <- res2$id
  ids <- ids[ids %in% rownames(mat)]
  if (length(ids) < 2) return(invisible(NULL))

  sub <- mat[ids, , drop = FALSE]
  sub <- log2(sub + 1)

  png(out_png, width = 2200, height = 2200, res = 250)
  pheatmap(sub, scale = "row", main = title, border_color = NA)
  dev.off()
}

read_enrich <- function(contrast, fname) {
  f <- file.path(BASE_ENR, contrast, fname)
  if (!file.exists(f)) return(NULL)
  
  df <- read_csv(f, show_col_types = FALSE)
  
  # === 1. Fix padj (standardize column name) ===
  # ClusterProfiler usually outputs 'p.adjust', but we need 'padj'
  if ("p.adjust" %in% names(df)) {
    df <- df %>% rename(padj = p.adjust)
  } else if ("qvalue" %in% names(df) && !"padj" %in% names(df)) {
    df <- df %>% rename(padj = qvalue)
  }
  
  # === 2. Fix term_name (standardize column name) ===
  # ClusterProfiler outputs 'Description', we need 'term_name'
  if ("Description" %in% names(df)) {
    df <- df %>% rename(term_name = Description)
  }
  
  # === 3. Fix k (Count) ===
  # ClusterProfiler outputs 'Count', we need 'k' for dot sizes
  if ("Count" %in% names(df)) {
    df <- df %>% rename(k = Count)
  }
  
  # === 4. Fix gene_ratio (Rename + Numeric Conversion) ===
  
  # Case A: Rename 'GeneRatio' to 'gene_ratio'
  if ("GeneRatio" %in% names(df)) {
    df <- df %>% rename(gene_ratio = GeneRatio)
  }
  
  # Case B: Convert fraction string "20/100" -> numeric 0.2
  # This is essential; otherwise, ggplot treats ratios as text and messes up the axis.
  if ("gene_ratio" %in% names(df) && is.character(df$gene_ratio)) {
    
    # Helper function to parse "a/b" strings
    parse_ratio <- function(x) {
      parts <- as.numeric(unlist(strsplit(as.character(x), "/")))
      if(length(parts) == 2 && parts[2] != 0) return(parts[1] / parts[2])
      return(NA) # Return NA if parsing fails
    }
    
    # Apply conversion to the column
    df$gene_ratio <- sapply(df$gene_ratio, parse_ratio)
  }
  
  if ("genes" %in% names(df)) {
    df <- df %>% rename(geneID = genes)
  } else if ("Genes" %in% names(df)) {
    df <- df %>% rename(geneID = Genes)
  } else if ("gene_id" %in% names(df)) {
    df <- df %>% rename(geneID = gene_id)
  } else if ("core_enrichment" %in% names(df)) {
    df <- df %>% rename(geneID = core_enrichment)
  }

  
  return(df)
}

plot_enrich_dot <- function(df, title, out_png, top_n = 20) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  # Select top terms
  df2 <- df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(top_n) %>%
    mutate(term_name = factor(term_name, levels = rev(unique(term_name))))

  p <- ggplot(df2, aes(x = gene_ratio, y = term_name, size = k, color = padj)) +
    geom_point() +
    labs(title = title, x = "Gene ratio", y = "") +
    theme_bw()
  ggsave(out_png, p, width = 9.5, height = 6.5, dpi = 300)
}

load_geneset <- function(path) {
  if (!file.exists(path)) return(character(0))
  x <- readLines(path, warn = FALSE)
  x <- x[nchar(x) > 0]
  unique(x)
}

plot_geneset_heatmap <- function(cnt, genes, title, out_png, max_genes = 80) {
  if (length(genes) == 0) return(invisible(NULL))

  mat <- as.data.frame(cnt)
  rownames(mat) <- mat$id
  mat$id <- NULL
  mat <- as.matrix(mat)

  # Intersect geneset with expression matrix
  genes_in <- genes[genes %in% rownames(mat)]
  if (length(genes_in) < 2) return(invisible(NULL))

  # Limit number of genes to avoid overcrowded heatmaps
  if (length(genes_in) > max_genes) genes_in <- genes_in[1:max_genes]

  sub <- mat[genes_in, , drop = FALSE]
  sub <- log2(sub + 1)

  png(out_png, width = 2400, height = 2400, res = 250)
  pheatmap(sub, scale = "row", main = title, border_color = NA)
  dev.off()
}

# ---------------------------
# Figure Summary:
# Fig2: Overview (Top 50 heatmap) + DEG counts
# Fig3: Sample correlation heatmap
# Fig4: Volcano plot
# Fig5: DEG up/down bar plot
# Fig6: UpSet plot across contrasts (intersecting DEG sets)
# Fig7: GO enrichment (dotplot; all/up/down)
# Fig8: KEGG pathway enrichment (dotplot; all/up/down)
# Fig10: RGA genes heatmap
# Fig13: PR1 genes heatmap
# Fig11/12: Templates for future use (specific pathway/gene plots)
# ---------------------------

all_deg_sets <- list()

for (ct in CONTRASTS) {
  message("Processing: ", ct)
  x <- read_deg_tables(ct)
  deg <- get_deg_sets(x$res, padj_cut = 0.001, lfc_cut = 1)
  all_deg_sets[[ct]] <- deg$all

  # Fig3: Sample correlation
  plot_sample_corr_heatmap(
    x$cnt,
    paste0("Fig3 Sample correlation: ", ct),
    file.path(OUTDIR, paste0(ct, "_Fig3_sample_corr.png"))
  )

  # Fig4: Volcano plot
  plot_volcano(
    deg$res2,
    paste0("Fig4 Volcano: ", ct),
    file.path(OUTDIR, paste0(ct, "_Fig4_volcano.png"))
  )

  # Fig5: DEG bar plot
  plot_deg_bar(
    deg,
    paste0("Fig5 DEG counts: ", ct),
    file.path(OUTDIR, paste0(ct, "_Fig5_deg_bar.png"))
  )

  # Fig2: Top genes heatmap (Overview of expression changes)
  plot_top_gene_heatmap(
    x$cnt, x$res, top_n = 50,
    title = paste0("Fig2 Top50 DE genes (by padj): ", ct),
    out_png = file.path(OUTDIR, paste0(ct, "_Fig2_top50_heatmap.png"))
  )

  # Fig7: GO enrichment
  go_all  <- read_enrich(ct, "GO_all.csv")
  go_up   <- read_enrich(ct, "GO_up.csv")
  go_down <- read_enrich(ct, "GO_down.csv")

  plot_enrich_dot(go_all,  paste0("Fig7 GO (all): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig7_GO_all_dot.png")))
  plot_enrich_dot(go_up,   paste0("Fig7 GO (up): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig7_GO_up_dot.png")))
  plot_enrich_dot(go_down, paste0("Fig7 GO (down): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig7_GO_down_dot.png")))

  # Fig8: KEGG pathway enrichment
  kp_all  <- read_enrich(ct, "KEGGpw_all.csv")
  kp_up   <- read_enrich(ct, "KEGGpw_up.csv")
  kp_down <- read_enrich(ct, "KEGGpw_down.csv")

  plot_enrich_dot(kp_all,  paste0("Fig8 KEGGpw (all): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig8_KEGGpw_all_dot.png")))
  plot_enrich_dot(kp_up,   paste0("Fig8 KEGGpw (up): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig8_KEGGpw_up_dot.png")))
  plot_enrich_dot(kp_down, paste0("Fig8 KEGGpw (down): ", ct),
                  file.path(OUTDIR, paste0(ct, "_Fig8_KEGGpw_down_dot.png")))

  # Fig10: RGA heatmap
  rga <- load_geneset(GENESET_RGA)
  plot_geneset_heatmap(
    x$cnt, rga,
    paste0("Fig10 RGA genes heatmap: ", ct),
    file.path(OUTDIR, paste0(ct, "_Fig10_RGA_heatmap.png")),
    max_genes = 80
  )

  # Fig13: PR1 heatmap
  pr1 <- load_geneset(GENESET_PR1)
  plot_geneset_heatmap(
    x$cnt, pr1,
    paste0("Fig13 PR1 genes heatmap: ", ct),
    file.path(OUTDIR, paste0(ct, "_Fig13_PR1_heatmap.png")),
    max_genes = 80
  )
}

# Fig6: UpSet plot (Intersection of DEG sets across contrasts)
upset_in <- fromList(all_deg_sets)
png(file.path(OUTDIR, "Fig6_DEG_UpSet_across_contrasts.png"),
    width = 2400, height = 1600, res = 250)
upset(upset_in, nsets = length(CONTRASTS), order.by = "freq")
dev.off()

# Note on Fig11/12 Templates:
# 1) Fig11: Select a KEGG pathway (e.g., Starch and sucrose metabolism map00500),
#    extract its genes, and plot heatmap or line plot from norm_counts.
# 2) Fig12: Select a set of "key genes to explain" and plot bar/line charts.
# =================================================================
# NEW SECTION: Generate Gene-Concept Network Plot (Cnetplot)
# =================================================================
library(ggraph)
library(tidygraph)
library(tidyr)
library(stringr)
library(GO.db)

# --- Pro Version: English Annotations & Custom Colors ---
plot_cnet_pro <- function(enrich_df, title, out_png, top_n = 5, gene_color = "#8DA0CB") {
  
  # 1. Data Check
  if (is.null(enrich_df) || nrow(enrich_df) == 0) return(invisible(NULL))
  
  # 2. Smart Label Selection
  if ("Description" %in% names(enrich_df)) {
    enrich_df <- enrich_df %>% mutate(plot_name = Description)
  } else if ("term_name" %in% names(enrich_df)) {
    enrich_df <- enrich_df %>% mutate(plot_name = term_name)
  } else {
    enrich_df <- enrich_df %>% mutate(plot_name = term_id) 
  }
  
  # 3. Translate GO IDs
  if (nrow(enrich_df) > 0 && str_starts(enrich_df$plot_name[1], "GO:")) {
    message("Found GO IDs. Translating...")
    go_ids <- enrich_df$plot_name
    english_terms <- tryCatch({
      AnnotationDbi::Term(go_ids)
    }, error = function(e) { return(NULL) })
    
    if (!is.null(english_terms)) {
      enrich_df$plot_name <- ifelse(is.na(english_terms), go_ids, english_terms)
    }
  }
  
  # 4. Text Formatting
  enrich_df$plot_name <- str_remove(enrich_df$plot_name, "^GO:\\d+ ") 
  enrich_df$plot_name <- str_to_sentence(enrich_df$plot_name)
  
  # 5. Filter Top N
  top_terms <- enrich_df %>%
    arrange(padj) %>%
    head(top_n) %>%
    dplyr::select(plot_name, geneID) 
  
  if (nrow(top_terms) == 0) return(invisible(NULL))
  
  # 6. Expand Data
  edges <- top_terms %>%
    mutate(geneID = str_split(geneID, "/")) %>% 
    unnest(geneID) %>%
    rename(from = plot_name, to = geneID)
  
  # 7. Create Graph
  graph <- as_tbl_graph(edges) %>%
    mutate(Degree = centrality_degree(mode = 'all'))
  
  # 8. Plotting
  p <- ggraph(graph, layout = "stress") + 
    
    # Edges
    geom_edge_link(alpha = 0.4, color = "grey70") + 
    
    # Pathways (Yellow/Beige)
    geom_node_point(aes(filter = name %in% top_terms$plot_name, size = Degree), 
                    color = "#FDBF6F", alpha = 1) + 
    
    # Genes (Red/Blue)
    geom_node_point(aes(filter = !(name %in% top_terms$plot_name)), 
                    color = gene_color, size = 2.5, alpha = 0.9) + 
    
    # Labels
    geom_node_text(aes(filter = name %in% top_terms$plot_name, 
                       label = str_wrap(name, width = 20)), 
                   repel = TRUE, size = 3.5, fontface = "bold", 
                   bg.color = "white", bg.r = 0.15) +
    

    theme_bw() +
    xlab("x") +
    ylab("y") + 
    
    labs(title = title, 
         caption = paste0("Yellow: Pathways | ", 
                          ifelse(gene_color=="#E41A1C", "Red: Up-regulated", "Blue: Down-regulated"))) +
    
    theme(
      plot.title = element_text(hjust = 0.5, face="bold", size = 14),
      plot.caption = element_text(size = 10, color = "grey30"),
      panel.grid.major = element_line(color = "grey90"), 
      panel.grid.minor = element_blank(),               
      axis.text = element_text(color = "grey50")        
    )
  
  print(p)
  ggsave(out_png, p, width = 11, height = 9, dpi = 300)
}

# --- Execution Loops ---

message("Generating Network Plots with Coordinates (GO)...")
for (ct in CONTRASTS) {
  go_up   <- read_enrich(ct, "GO_up.csv")
  go_down <- read_enrich(ct, "GO_down.csv")
  
  if (!is.null(go_up)) {
    plot_cnet_pro(go_up, 
                  title = paste0("Upregulated Gene Network: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig9_Network_UP_English.png")), 
                  top_n = 5, gene_color = "#E41A1C") 
  }
  
  if (!is.null(go_down)) {
    plot_cnet_pro(go_down, 
                  title = paste0("Downregulated Gene Network: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig9_Network_DOWN_English.png")), 
                  top_n = 5, gene_color = "#377EB8") 
  }
}

message("Generating Network Plots with Coordinates (KEGG)...")
for (ct in CONTRASTS) {
  kp_up   <- read_enrich(ct, "KEGGpw_up.csv")
  kp_down <- read_enrich(ct, "KEGGpw_down.csv")
  
  if (!is.null(kp_up)) {
    plot_cnet_pro(kp_up, 
                  title = paste0("KEGG Network UP: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig11_KEGG_Network_UP.png")), 
                  top_n = 5, gene_color = "#E41A1C")
  }
  
  if (!is.null(kp_down)) {
    plot_cnet_pro(kp_down, 
                  title = paste0("KEGG Network DOWN: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig11_KEGG_Network_DOWN.png")), 
                  top_n = 5, gene_color = "#377EB8")
  }
}

message("Done! Networks now have axes and grids.")
# --- Loop to generate KEGG Network Plots (Red/Blue + English Annotations) ---
message("Generating English Annotation Network Plots (KEGG)...")

for (ct in CONTRASTS) {
  kp_up   <- read_enrich(ct, "KEGGpw_up.csv")
  kp_down <- read_enrich(ct, "KEGGpw_down.csv")
  
  # 1. Upregulated -> Red
  if (!is.null(kp_up)) {
    plot_cnet_pro(kp_up, 
                  title = paste0("KEGG Network UP: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig11_KEGG_Network_UP.png")), 
                  top_n = 5,
                  gene_color = "#E41A1C")
  }
  
  # 2. Downregulated -> Blue
  if (!is.null(kp_down)) {
    plot_cnet_pro(kp_down, 
                  title = paste0("KEGG Network DOWN: ", ct), 
                  out_png = file.path(OUTDIR, paste0(ct, "_Fig11_KEGG_Network_DOWN.png")), 
                  top_n = 5,
                  gene_color = "#377EB8")
  }
}

# =================================================================
# NEW SECTION: Generate Professional Volcano Plot
# =================================================================
library(ggrepel) 
library(dplyr)

plot_volcano_pro <- function(res, title, out_file) {
  
  # 1. Prepare Data
  df <- as.data.frame(res)
  df$geneID <- rownames(df) 
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  
  # 2. Define Thresholds
  padj_cut <- 0.001
  lfc_cut <- 1
  
  # 3. Grouping
  df$Group <- "NS"
  df$Group[df$padj < padj_cut & df$log2FoldChange > lfc_cut] <- "Up"
  df$Group[df$padj < padj_cut & df$log2FoldChange < -lfc_cut] <- "Down"
  df$Group <- factor(df$Group, levels = c("Up", "Down", "NS"))
  
  # 4. Select Top 10 genes for labeling
  top_genes <- df %>%
    filter(Group != "NS") %>%
    arrange(padj) %>%
    head(10)
  
  # 5. Plotting
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Group)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey80")) + 
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_text_repel(data = top_genes, aes(label = geneID),
                    size = 3, box.padding = 0.5, max.overlaps = Inf,
                    color = "black", show.legend = FALSE) + 
    labs(title = title, x = "log2 Fold Change", y = "-log10 (adj. P-value)",
         subtitle = paste0("p-adj < ", padj_cut, ", |log2FC| > ", lfc_cut)) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top")
  
  print(p)
  ggsave(out_file, p, width = 8, height = 7, dpi = 300)
}

message("Generating Professional Volcano Plots...")
for (ct in CONTRASTS) {
  x <- read_deg_tables(ct) 
  if (!is.null(x$res)) {
    plot_volcano_pro(x$res, paste0("Volcano: ", ct), file.path(OUTDIR, paste0(ct, "_Fig4_volcano_pro.png")))
  }
}

# =================================================================
# NEW SECTION: Generate Combined Summary Plot (Figure 5 Style)
# =================================================================
summary_df <- data.frame()

# Iterate through contrasts to collect data
for (ct in CONTRASTS) {
  x <- read_deg_tables(ct)
  deg <- get_deg_sets(x$res, padj_cut = 0.001, lfc_cut = 1)
  summary_df <- rbind(summary_df, data.frame(Contrast = ct, Direction = "Upregulated", Count = length(deg$up)))
  summary_df <- rbind(summary_df, data.frame(Contrast = ct, Direction = "Downregulated", Count = length(deg$down)))
}

# Set factor levels for correct ordering
summary_df$Contrast <- factor(summary_df$Contrast, levels = CONTRASTS)
summary_df$Direction <- factor(summary_df$Direction, levels = c("Upregulated", "Downregulated"))

# Create Bar Plot
p_combined <- ggplot(summary_df, aes(x = Contrast, y = Count, fill = Direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) + 
  labs(title = "Differentially Expressed Genes by Contrast",
       subtitle = "padj < 0.001 and |log2FC| >= 1", x = "Comparison", y = "Number of DEGs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), 
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top", panel.grid.major.x = element_blank())

print(p_combined)
out_file <- file.path(OUTDIR, "Fig5_Combined_DEG_Barplot.png")
ggsave(out_file, p_combined, width = 10, height = 6, dpi = 300)

message("Combined barplot saved to: ", out_file)
message("All plotting tasks completed successfully!")

