# ============================================================
# volcano_plot_combined.R
# CT6・CT18を統合してWT(n=4) vs KO(n=4)のボルケーノプロットを作成する
#
# 【なぜ統合するか】
#   小脳ではGpr176の発現リズムが弱い → CT6/CT18を分ける生物学的根拠が薄い
#   統合することでサンプル数が n=2 → n=4 になり、検定の信頼性が上がる
#
# 【volcano_plot.R との違い】
#   n=4になるのでBH補正（多重比較補正）を試みる
#   ただし遺伝子数が多い場合はBH補正後も有意遺伝子が出ないことがある
#   → スクリプトが自動判定してフォールバック
# ============================================================

library(readxl)
library(ggplot2)
library(ggrepel)

# ---- パラメータ ----
INPUT_FILE  <- "allpresence.xlsx"
SHEET       <- "Allプレゼンスのみ"
FC_CUTOFF   <- 1.5
PVAL_CUTOFF <- 0.05
TOP_N_LABEL <- 30  # CT別より多め（統合なので候補が増える）

# ---- シグナル値の列番号 ----
# CT6・CT18を区別せず、WT全4サンプル・KO全4サンプルとして扱う
WT_COLS <- c(2, 5, 8, 11)    # 1-WT-CT6, 2-WT-CT6, 3-WT-CT18, 4-WT-CT18
KO_COLS <- c(14, 17, 20, 23) # 5-KO-CT6, 6-KO-CT6, 7-KO-CT18, 8-KO-CT18

GENE_SYMBOL_COL <- "Gene Symbol"
PROBE_ID_COL    <- "Probe Set ID"

# ---- データ読み込み ----
message("Reading: ", INPUT_FILE, "  sheet: ", SHEET)
raw <- read_excel(INPUT_FILE, sheet = SHEET)

probe_ids    <- raw[[PROBE_ID_COL]]
gene_symbols <- raw[[GENE_SYMBOL_COL]]

to_mat <- function(cols) {
  m <- as.matrix(raw[, cols])
  apply(m, 2, as.numeric)
}

wt_mat <- to_mat(WT_COLS)
ko_mat <- to_mat(KO_COLS)

# ---- 統計計算 ----
# n=4 vs n=4 の Welch t検定
log2fc <- log2(rowMeans(ko_mat, na.rm = TRUE) / rowMeans(wt_mat, na.rm = TRUE))

pvals <- vapply(seq_len(nrow(raw)), function(i) {
  tryCatch(
    t.test(ko_mat[i, ], wt_mat[i, ])$p.value,
    error = function(e) NA_real_
  )
}, numeric(1))

# BH補正（Benjamini-Hochberg法）による多重比較補正
# 数万プローブを同時に検定するため、偽陽性を制御する必要がある
# BH補正: FDR（偽発見率）を指定した割合以下に抑える
# n=4でも検出力が低い場合はpadj ≈ 1 になることがある
padj <- p.adjust(pvals, method = "BH")

# ---- 結果テーブル ----
res <- data.frame(
  probe_id    = probe_ids,
  gene_symbol = gene_symbols,
  log2FC      = log2fc,
  pval        = pvals,
  padj        = padj,
  stringsAsFactors = FALSE
)
res <- res[complete.cases(res), ]

# BH補正版とraw p値版の両方で有意遺伝子を分類
res$sig_bh <- "NS"
res$sig_bh[res$log2FC >=  log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Up in KO"
res$sig_bh[res$log2FC <= -log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Down in KO"

res$sig_raw <- "NS"
res$sig_raw[res$log2FC >=  log2(FC_CUTOFF) & res$pval < PVAL_CUTOFF] <- "Up in KO"
res$sig_raw[res$log2FC <= -log2(FC_CUTOFF) & res$pval < PVAL_CUTOFF] <- "Down in KO"

# 有意遺伝子数を確認
message(sprintf("BH補正あり → Up: %d  Down: %d",
                sum(res$sig_bh == "Up in KO"), sum(res$sig_bh == "Down in KO")))
message(sprintf("raw p値    → Up: %d  Down: %d",
                sum(res$sig_raw == "Up in KO"), sum(res$sig_raw == "Down in KO")))

# BH補正で10個以上有意遺伝子があればBH採用、少なければraw p値にフォールバック
use_bh <- sum(res$sig_bh != "NS") >= 10
res$sig <- if (use_bh) res$sig_bh else res$sig_raw

# Y軸の値とラベルをBH/rawで切り替え
y_label <- if (use_bh) {
  expression(-log[10]~"(adjusted p-value, BH)")
} else {
  expression(-log[10]~"(p-value, raw)")
}
res$y          <- if (use_bh) -log10(res$padj) else -log10(res$pval)
subtitle_note  <- if (use_bh) "BH補正適用 (n=4)" else "raw p値使用 (n=4)"

res$sig <- factor(res$sig, levels = c("Up in KO", "Down in KO", "NS"))

# ---- ラベル遺伝子の選定 ----
sig_rows  <- res[res$sig != "NS", ]
n_top     <- min(TOP_N_LABEL, nrow(sig_rows))
top_genes <- sig_rows[order(if (use_bh) sig_rows$padj else sig_rows$pval)[seq_len(n_top)], ]
top_genes$label <- ifelse(
  is.na(top_genes$gene_symbol) | top_genes$gene_symbol %in% c("---", ""),
  top_genes$probe_id,
  top_genes$gene_symbol
)

# ---- プロット ----
color_map <- c("Up in KO" = "#E41A1C", "Down in KO" = "#377EB8", "NS" = "grey70")

p <- ggplot(res, aes(x = log2FC, y = y, color = sig)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-log2(FC_CUTOFF), log2(FC_CUTOFF)),
             linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text_repel(
    data = top_genes,
    aes(x = log2FC, y = y, label = label),
    color = "black", size = 2.5, max.overlaps = 40,
    box.padding = 0.3, point.padding = 0.2
  ) +
  scale_color_manual(values = color_map, name = NULL) +
  labs(
    title    = "Volcano Plot: KO vs WT  [CT6 + CT18 統合, n=4]",
    subtitle = sprintf("FC > %.1fx  |  p < %.2f (%s)  |  Up: %d  Down: %d",
                       FC_CUTOFF, PVAL_CUTOFF, subtitle_note,
                       sum(res$sig == "Up in KO"),
                       sum(res$sig == "Down in KO")),
    x = expression(log[2]~"(KO / WT)"),
    y = y_label
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

ggsave("volcano_combined.pdf", plot = p, width = 7, height = 6)
message("Saved: volcano_combined.pdf")
print(p)

# ---- 結果CSV ----
# sig_bh: BH補正での有意判定  sig_raw: raw p値での有意判定
out <- res[order(if (use_bh) res$padj else res$pval),
           c("probe_id", "gene_symbol", "log2FC", "pval", "padj", "sig_bh", "sig_raw")]
write.csv(out, "volcano_results_combined.csv", row.names = FALSE)
message("Saved: volcano_results_combined.csv")

# Up遺伝子リスト（いずれかの基準で有意なもの）
up <- out[out$sig_bh == "Up in KO" | out$sig_raw == "Up in KO", ]
write.csv(up, "up_genes_combined.csv", row.names = FALSE)
message(sprintf("Saved: up_genes_combined.csv  (%d probes)", nrow(up)))

message("\n--- Up遺伝子 Top20 ---")
print(head(up[, c("gene_symbol", "log2FC", "pval", "padj")], 20))
