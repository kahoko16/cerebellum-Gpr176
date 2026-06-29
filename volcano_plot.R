# Volcano plot for Gpr176 KO vs WT microarray data
# Data: 8 samples (2x WT-CT6, 2x WT-CT18, 2x KO-CT6, 2x KO-CT18)
# Required packages: readxl, ggplot2, ggrepel

library(readxl)
library(ggplot2)
library(ggrepel)

# ---- Parameters ----
INPUT_FILE   <- "ALL-presence.xlsx"   # path to Excel file
SHEET        <- 1                      # sheet number or name
FC_CUTOFF    <- 1.5                    # fold-change threshold (linear)
PVAL_CUTOFF  <- 0.05                   # adjusted p-value threshold
TOP_N_LABEL  <- 20                     # number of top genes to label
OUTPUT_FILE  <- "volcano_plot.pdf"

# ---- Column indices (1-based) ----
# Signal columns: cols 2, 5, 8, 11 (WT-CT6 x2, WT-CT18 x2) and
#                 cols 14, 17, 20, 23 (KO-CT6 x2, KO-CT18 x2)
# Each sample block: [signal, Detection, Detection p-value]
WT_COLS <- c(2, 5, 8, 11)   # 1-WT-CT6, 2-WT-CT6, 3-WT-CT18, 4-WT-CT18
KO_COLS <- c(14, 17, 20, 23) # 5-KO-CT6, 6-KO-CT6, 7-KO-CT18, 8-KO-CT18

GENE_SYMBOL_COL <- "Gene Symbol"
PROBE_ID_COL    <- "Probe Set ID"

# ---- Load data ----
message("Reading: ", INPUT_FILE)
raw <- read_excel(INPUT_FILE, sheet = SHEET)

probe_ids    <- raw[[PROBE_ID_COL]]
gene_symbols <- raw[[GENE_SYMBOL_COL]]

wt_mat <- as.matrix(raw[, WT_COLS])
ko_mat <- as.matrix(raw[, KO_COLS])

# Ensure numeric
wt_mat <- apply(wt_mat, 2, as.numeric)
ko_mat <- apply(ko_mat, 2, as.numeric)

# ---- Compute log2FC and p-value (Welch t-test per probe) ----
log2fc <- log2(rowMeans(ko_mat, na.rm = TRUE) / rowMeans(wt_mat, na.rm = TRUE))

pvals <- vapply(seq_len(nrow(raw)), function(i) {
  tryCatch(
    t.test(ko_mat[i, ], wt_mat[i, ])$p.value,
    error = function(e) NA_real_
  )
}, numeric(1))

# BH correction
padj <- p.adjust(pvals, method = "BH")

# ---- Build results table ----
res <- data.frame(
  probe_id    = probe_ids,
  gene_symbol = gene_symbols,
  log2FC      = log2fc,
  pval        = pvals,
  padj        = padj,
  stringsAsFactors = FALSE
)

# Remove rows with NA
res <- res[complete.cases(res), ]

# Significance classification
res$sig <- "NS"
res$sig[abs(res$log2FC) >= log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Significant"
res$sig[res$log2FC >=  log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Up in KO"
res$sig[res$log2FC <= -log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Down in KO"

res$sig <- factor(res$sig, levels = c("Up in KO", "Down in KO", "NS"))

# Top genes to label (by adjusted p-value among significant)
sig_rows  <- res[res$sig != "NS", ]
top_genes <- sig_rows[order(sig_rows$padj)[seq_len(min(TOP_N_LABEL, nrow(sig_rows)))], ]
label_col <- ifelse(top_genes$gene_symbol %in% c("---", "", NA), top_genes$probe_id, top_genes$gene_symbol)
top_genes$label <- label_col

# ---- Plot ----
color_map <- c("Up in KO" = "#E41A1C", "Down in KO" = "#377EB8", "NS" = "grey70")

p <- ggplot(res, aes(x = log2FC, y = -log10(padj), color = sig)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-log2(FC_CUTOFF), log2(FC_CUTOFF)), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(x = log2FC, y = -log10(padj), label = label),
    color = "black", size = 2.5, max.overlaps = 30,
    box.padding = 0.3, point.padding = 0.2
  ) +
  scale_color_manual(values = color_map, name = NULL) +
  labs(
    title    = "Volcano Plot: Gpr176 KO vs WT",
    subtitle = sprintf("FC cutoff: %.1fx  |  FDR cutoff: %.2f  |  Significant: %d probes",
                       FC_CUTOFF, PVAL_CUTOFF, sum(res$sig != "NS")),
    x        = expression(log[2]~"(KO / WT)"),
    y        = expression(-log[10]~"(adjusted p-value)")
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position  = "top",
    plot.subtitle    = element_text(size = 9, color = "grey40")
  )

# Save
ggsave(OUTPUT_FILE, plot = p, width = 7, height = 6)
message("Saved: ", OUTPUT_FILE)

# Save results table
write.csv(res[order(res$padj), ], "volcano_results.csv", row.names = FALSE)
message("Saved: volcano_results.csv")

print(p)
