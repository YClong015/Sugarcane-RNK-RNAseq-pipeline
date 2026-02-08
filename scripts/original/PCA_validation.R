library(tximport)
library(DESeq2)
library(ggplot2)

meta <- read.csv(
  "/QRISdata/Q9062/08_deseq2_R570/sample_metadata.csv",
  stringsAsFactors = TRUE
)

files <- file.path(
  "/QRISdata/Q9062/07_salmon_R570/quants",
  meta$sample,
  "quant.sf"
)
names(files) <- meta$sample

tx2gene <- read.delim(
  "/QRISdata/Q9062/06_ref/tx2gene_R570_1os2g.tsv",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(tx2gene) <- c("TXNAME", "GENEID")

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

dds0 <- DESeqDataSetFromTximport(
  txi,
  colData = meta,
  design = ~ genotype + time + treatment
)

dds0 <- dds0[rowSums(counts(dds0)) > 10, ]
vsd0 <- vst(dds0, blind = TRUE)

pca <- plotPCA(
  vsd0,
  intgroup = c("genotype", "time", "treatment"),
  returnData = TRUE
)

percentVar <- round(100 * attr(pca, "percentVar"))
pca$sample <- rownames(pca)
pca$group <- paste(pca$genotype, pca$time, pca$treatment, sep = ":")

pca$group <- factor(
  pca$group,
  levels = c(
    "Q208:7d:C",
    "Q208:7d:RKN",
    "Q208:21d:C",
    "Q208:21d:RKN",
    "Q208:12w:C",
    "Q208:12w:RKN",
    "SES208:12w:C",
    "SES208:12w:RKN"
  )
)

cols <- c(
  "Q208:7d:C"      = "#56B4E9",
  "Q208:7d:RKN"    = "#0072B2",
  "Q208:21d:C"     = "#009E73",
  "Q208:21d:RKN"   = "#004D40",
  "Q208:12w:C"     = "#E69F00",
  "Q208:12w:RKN"   = "#D55E00",
  "SES208:12w:C"   = "#CC79A7",
  "SES208:12w:RKN" = "#9B0066"
)

xr <- range(pca$PC1, na.rm = TRUE)
yr <- range(pca$PC2, na.rm = TRUE)
xpad <- diff(xr) * 0.18
ypad <- diff(yr) * 0.18
xlim <- c(xr[1] - xpad, xr[2] + xpad)
ylim <- c(yr[1] - ypad, yr[2] + ypad)

g <- ggplot(pca, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = group, shape = treatment),
    size = 3.2,
    stroke = 0.3
  ) +
  stat_ellipse(
    aes(fill = group),
    geom = "polygon",
    level = 0.95,
    alpha = 0.28,      
    color = NA
  ) +
  stat_ellipse(
    aes(color = group),
    level = 0.95,
    linewidth = 1.2    
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = "PCA of RNA-seq samples",
    subtitle = "VST-transformed counts with 95% confidence ellipses",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Group",
    fill = "Group",
    shape = "Treatment"
  ) +
  theme_classic(base_size = 14)

print(g)
ggsave(
  "/QRISdata/Q9062/08_deseq2_R570/QC_PCA_publication.pdf",
  g,
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  "/QRISdata/Q9062/08_deseq2_R570/QC_PCA_publication.png",
  g,
  width = 8,
  height = 6,
  dpi = 300
)


