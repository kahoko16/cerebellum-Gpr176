# Volcano plot for Gpr176 KO vs WT microarray data
# Produces separate plots for CT6 and CT18
# Required packages: readxl, ggplot2, ggrepel

library(readxl)
library(ggplot2)
library(ggrepel)

# ---- Parameters ----
INPUT_FILE  <- "ALL-presence.xlsx"      # Excelファイルのパス
SHEET       <- "Allプレゼンスのみ"       # シート名
FC_CUTOFF   <- 1.5                       # fold-change閾値（線形）
PVAL_CUTOFF <- 0.05                      # adjusted p-value閾値
TOP_N_LABEL <- 20                        # ラベル表示する上位遺伝子数

# ---- Signal column indices (1-based) ----
# 各サンプルブロック: [signal, Detection, Detection p-value] の3列構成
# 1-WT-CT6  : col 2
# 2-WT-CT6  : col 5
# 3-WT-CT18 : col 8
# 4-WT-CT18 : col 11
# 5-KO-CT6  : col 14
# 6-KO-CT6  : col 17
# 7-KO-CT18 : col 20
# 8-KO-CT18 : col 23
WT_CT6_COLS  <- c(2, 5)
WT_CT18_COLS <- c(8, 11)
KO_CT6_COLS  <- c(14, 17)
KO_CT18_COLS <- c(20, 23)

GENE_SYMBOL_COL <- "Gene Symbol"
PROBE_ID_COL    <- "Probe Set ID"

# ---- Load data ----
message("Reading: ", INPUT_FILE, "  sheet: ", SHEET)
raw <- read_excel(INPUT_FILE, sheet = SHEET)

probe_ids    <- raw[[PROBE_ID_COL]]
gene_symbols <- raw[[GENE_SYMBOL_COL]]

to_mat <- function(cols) {
  m <- as.matrix(raw[, cols])
  apply(m, 2, as.numeric)
}

wt_ct6  <- to_mat(WT_CT6_COLS)
wt_ct18 <- to_mat(WT_CT18_COLS)
ko_ct6  <- to_mat(KO_CT6_COLS)
ko_ct18 <- to_mat(KO_CT18_COLS)

# ---- Helper: compute stats and build result table ----
compute_results <- function(ko_mat, wt_mat) {
  log2fc <- log2(rowMeans(ko_mat, na.rm = TRUE) / rowMeans(wt_mat, na.rm = TRUE))

  pvals <- vapply(seq_len(nrow(raw)), function(i) {
    tryCatch(
      t.test(ko_mat[i, ], wt_mat[i, ])$p.value,
      error = function(e) NA_real_
    )
  }, numeric(1))

  # n=2のためBH補正は使わず raw p値を使用
  padj <- pvals

  res <- data.frame(
    probe_id    = probe_ids,
    gene_symbol = gene_symbols,
    log2FC      = log2fc,
    pval        = pvals,
    padj        = padj,
    stringsAsFactors = FALSE
  )
  res <- res[complete.cases(res), ]

  res$sig <- "NS"
  res$sig[res$log2FC >=  log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Up in KO"
  res$sig[res$log2FC <= -log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Down in KO"
  res$sig <- factor(res$sig, levels = c("Up in KO", "Down in KO", "NS"))

  res
}

# ---- Helper: draw volcano plot ----
draw_volcano <- function(res, title_label) {
  color_map <- c("Up in KO" = "#E41A1C", "Down in KO" = "#377EB8", "NS" = "grey70")

  sig_rows <- res[res$sig != "NS", ]
  n_top    <- min(TOP_N_LABEL, nrow(sig_rows))
  top_genes <- sig_rows[order(sig_rows$padj)[seq_len(n_top)], ]
  top_genes$label <- ifelse(
    is.na(top_genes$gene_symbol) | top_genes$gene_symbol %in% c("---", ""),
    top_genes$probe_id,
    top_genes$gene_symbol
  )

  ggplot(res, aes(x = log2FC, y = -log10(padj), color = sig)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed",
               color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-log2(FC_CUTOFF), log2(FC_CUTOFF)),
               linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_text_repel(
      data = top_genes,
      aes(x = log2FC, y = -log10(padj), label = label),
      color = "black", size = 2.5, max.overlaps = 30,
      box.padding = 0.3, point.padding = 0.2
    ) +
    scale_color_manual(values = color_map, name = NULL) +
    labs(
      title    = sprintf("Volcano Plot: KO vs WT  [%s]", title_label),
      subtitle = sprintf("FC > %.1fx  |  p < %.2f (raw, n=2)  |  Up: %d  Down: %d",
                         FC_CUTOFF, PVAL_CUTOFF,
                         sum(res$sig == "Up in KO"),
                         sum(res$sig == "Down in KO")),
      x = expression(log[2]~"(KO / WT)"),
      y = expression(-log[10]~"(p-value, raw)")
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "top",
      plot.subtitle   = element_text(size = 9, color = "grey40")
    )
}

# ---- CT6 解析 ----
message("--- CT6 解析 ---")
res_ct6 <- compute_results(ko_ct6, wt_ct6)
p_ct6   <- draw_volcano(res_ct6, "CT6")

ggsave("volcano_CT6.pdf", plot = p_ct6, width = 7, height = 6)
write.csv(res_ct6[order(res_ct6$padj), ], "volcano_results_CT6.csv", row.names = FALSE)
message("Saved: volcano_CT6.pdf / volcano_results_CT6.csv")

# ---- CT18 解析 ----
message("--- CT18 解析 ---")
res_ct18 <- compute_results(ko_ct18, wt_ct18)
p_ct18   <- draw_volcano(res_ct18, "CT18")

ggsave("volcano_CT18.pdf", plot = p_ct18, width = 7, height = 6)
write.csv(res_ct18[order(res_ct18$padj), ], "volcano_results_CT18.csv", row.names = FALSE)
message("Saved: volcano_CT18.pdf / volcano_results_CT18.csv")

# ---- 並べて表示 ----
library(patchwork)
combined <- p_ct6 + p_ct18 + plot_layout(ncol = 2)
ggsave("volcano_CT6_CT18_combined.pdf", plot = combined, width = 14, height = 6)
message("Saved: volcano_CT6_CT18_combined.pdf")

print(combined)

# ---- 網羅的統合リスト（CT6・CT18全プローブをマージ） ----
# CT6とCT18の結果を全プローブで横並びにする
all6  <- res_ct6[,  c("probe_id", "gene_symbol", "log2FC", "pval", "sig")]
all18 <- res_ct18[, c("probe_id", "gene_symbol", "log2FC", "pval", "sig")]
colnames(all6)  <- c("probe_id", "gene_symbol", "log2FC_CT6",  "pval_CT6",  "sig_CT6")
colnames(all18) <- c("probe_id", "gene_symbol", "log2FC_CT18", "pval_CT18", "sig_CT18")

merged <- merge(all6, all18, by = c("probe_id", "gene_symbol"))

# CT6・CT18それぞれでのUp/Down判定
merged$up_CT6   <- merged$sig_CT6  == "Up in KO"
merged$up_CT18  <- merged$sig_CT18 == "Up in KO"
merged$down_CT6 <- merged$sig_CT6  == "Down in KO"
merged$down_CT18<- merged$sig_CT18 == "Down in KO"

# カテゴリ分類
merged$category <- "NS"
merged$category[merged$up_CT6   & merged$up_CT18  ] <- "Up_both"
merged$category[merged$down_CT6 & merged$down_CT18] <- "Down_both"
merged$category[merged$up_CT6   & !merged$up_CT18 & !merged$down_CT18] <- "Up_CT6only"
merged$category[merged$up_CT18  & !merged$up_CT6  & !merged$down_CT6 ] <- "Up_CT18only"
merged$category[merged$down_CT6 & !merged$down_CT18& !merged$up_CT18 ] <- "Down_CT6only"
merged$category[merged$down_CT18& !merged$down_CT6 & !merged$up_CT6  ] <- "Down_CT18only"

# pval_mean（CT6・CT18の幾何平均）でソート
merged$pval_mean <- sqrt(merged$pval_CT6 * merged$pval_CT18)

# 全遺伝子テーブル保存
out_cols <- c("probe_id", "gene_symbol",
              "log2FC_CT6", "pval_CT6", "sig_CT6",
              "log2FC_CT18", "pval_CT18", "sig_CT18",
              "category", "pval_mean")
write.csv(merged[order(merged$pval_mean), out_cols], "all_genes_CT6_CT18.csv", row.names = FALSE)
message(sprintf("Saved: all_genes_CT6_CT18.csv  (%d probes total)", nrow(merged)))

# ---- Up遺伝子リスト（リガンド候補） ----
up_genes <- merged[merged$category %in% c("Up_both", "Up_CT6only", "Up_CT18only"), ]
up_genes <- up_genes[order(up_genes$pval_mean), out_cols]

write.csv(up_genes, "up_genes_KO_comprehensive.csv", row.names = FALSE)
message(sprintf("Saved: up_genes_KO_comprehensive.csv  (%d probes)", nrow(up_genes)))

# サマリー
message("\n===== Up遺伝子 サマリー（リガンド候補） =====")
message(sprintf("  CT6・CT18 両方でUp : %d", sum(merged$category == "Up_both")))
message(sprintf("  CT6のみでUp        : %d", sum(merged$category == "Up_CT6only")))
message(sprintf("  CT18のみでUp       : %d", sum(merged$category == "Up_CT18only")))
message(sprintf("  合計               : %d", nrow(up_genes)))

message("\n--- Up遺伝子 Top30（pval_meanでソート） ---")
print(head(up_genes[, c("gene_symbol", "log2FC_CT6", "log2FC_CT18",
                         "pval_CT6", "pval_CT18", "category")], 30))
