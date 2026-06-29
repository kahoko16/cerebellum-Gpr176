# ============================================================
# volcano_plot.R
# マイクロアレイデータからCT6・CT18別ボルケーノプロットを作成する
#
# 【ボルケーノプロットとは】
#   X軸: log2(発現変化倍率) ... 右がKOでUp、左がKOでDown
#   Y軸: -log10(p値)        ... 上ほど統計的に有意
#   左右かつ上にある点 = 発現変化が大きく、かつ有意な遺伝子
#
# 【なぜlog2スケールか】
#   2倍と0.5倍（1/2倍）を±1で対称に表現できるため
#   例: KO/WT = 2倍 → log2(2) = +1
#       KO/WT = 0.5倍 → log2(0.5) = -1
# ============================================================

# ---- パッケージ読み込み ----
# readxl  : Excelファイルを読み込む
# ggplot2 : グラフ作成の定番パッケージ
# ggrepel : ラベルが重ならないよう自動配置してくれる
library(readxl)
library(ggplot2)
library(ggrepel)

# ---- パラメータ（ここを変えると解析条件が変わる） ----
INPUT_FILE  <- "allpresence.xlsx"   # 読み込むExcelファイル名
SHEET       <- "Allプレゼンスのみ"  # 使用するシート名
FC_CUTOFF   <- 1.5    # 発現変化の閾値（線形）: 1.5倍以上の変化を「有意」とみなす
PVAL_CUTOFF <- 0.05   # p値の閾値: 0.05未満を「有意」とみなす（n=2のためraw p値を使用）
TOP_N_LABEL <- 20     # ラベルを表示する上位遺伝子数

# ---- シグナル値の列番号（1始まり）----
# Excelの列構成:
#   1列目: Probe Set ID
#   2列目: 1-WT-CT6のシグナル値  3列目: Detection  4列目: Detection p値
#   5列目: 2-WT-CT6のシグナル値  6列目: Detection  7列目: Detection p値
#   ... 以下同様に3列ずつ
# → シグナル値は 2, 5, 8, 11, 14, 17, 20, 23 列目
WT_CT6_COLS  <- c(2, 5)    # WT-CT6の2サンプル
WT_CT18_COLS <- c(8, 11)   # WT-CT18の2サンプル
KO_CT6_COLS  <- c(14, 17)  # KO-CT6の2サンプル
KO_CT18_COLS <- c(20, 23)  # KO-CT18の2サンプル

GENE_SYMBOL_COL <- "Gene Symbol"  # 遺伝子名の列名
PROBE_ID_COL    <- "Probe Set ID" # プローブIDの列名

# ---- Excelデータの読み込み ----
message("Reading: ", INPUT_FILE, "  sheet: ", SHEET)
raw <- read_excel(INPUT_FILE, sheet = SHEET)

# 遺伝子名とプローブIDを取り出す
probe_ids    <- raw[[PROBE_ID_COL]]
gene_symbols <- raw[[GENE_SYMBOL_COL]]

# 指定した列番号のデータを数値の行列に変換する関数
to_mat <- function(cols) {
  m <- as.matrix(raw[, cols])
  apply(m, 2, as.numeric)  # 文字として読まれた場合に備えてas.numericで変換
}

wt_ct6  <- to_mat(WT_CT6_COLS)
wt_ct18 <- to_mat(WT_CT18_COLS)
ko_ct6  <- to_mat(KO_CT6_COLS)
ko_ct18 <- to_mat(KO_CT18_COLS)

# ---- 統計計算・結果テーブル作成の関数 ----
compute_results <- function(ko_mat, wt_mat) {

  # log2(KO平均 / WT平均) を各プローブで計算
  # rowMeans: 行（=プローブ）ごとに列（=サンプル）の平均を取る
  log2fc <- log2(rowMeans(ko_mat, na.rm = TRUE) / rowMeans(wt_mat, na.rm = TRUE))

  # Welch t検定: 各プローブについてKO群とWT群の平均値が異なるかを検定
  # vapply: 各プローブ（行）に対してループ処理を行う高速版sapply
  # tryCatch: エラーが起きても止まらないようにする（分散が0の場合など）
  pvals <- vapply(seq_len(nrow(raw)), function(i) {
    tryCatch(
      t.test(ko_mat[i, ], wt_mat[i, ])$p.value,
      error = function(e) NA_real_  # エラーならNAを返す
    )
  }, numeric(1))

  # 【注意】n=2では多重比較補正（BH法）を使うとp値がすべて1になる
  # 理由: 検出力が低すぎて、数万プローブを補正すると全部非有意になる
  # → raw p値（補正なし）をそのまま使用する
  padj <- pvals  # n=2のためBH補正なし

  # 結果をデータフレームにまとめる
  res <- data.frame(
    probe_id    = probe_ids,
    gene_symbol = gene_symbols,
    log2FC      = log2fc,
    pval        = pvals,
    padj        = padj,
    stringsAsFactors = FALSE
  )

  # NAが含まれる行を除外（t検定が計算できなかったプローブ）
  res <- res[complete.cases(res), ]

  # 有意性の分類
  # FC閾値を log2スケールに変換: 1.5倍 → log2(1.5) ≈ 0.585
  res$sig <- "NS"  # まず全部 Not Significant
  res$sig[res$log2FC >=  log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Up in KO"
  res$sig[res$log2FC <= -log2(FC_CUTOFF) & res$padj < PVAL_CUTOFF] <- "Down in KO"
  res$sig <- factor(res$sig, levels = c("Up in KO", "Down in KO", "NS"))

  res
}

# ---- ボルケーノプロット描画の関数 ----
draw_volcano <- function(res, title_label) {
  color_map <- c("Up in KO" = "#E41A1C", "Down in KO" = "#377EB8", "NS" = "grey70")

  # p値が低い順に上位N遺伝子を選んでラベル表示
  sig_rows  <- res[res$sig != "NS", ]
  n_top     <- min(TOP_N_LABEL, nrow(sig_rows))
  top_genes <- sig_rows[order(sig_rows$padj)[seq_len(n_top)], ]

  # 遺伝子名が "---"（未同定）の場合はプローブIDで代替
  top_genes$label <- ifelse(
    is.na(top_genes$gene_symbol) | top_genes$gene_symbol %in% c("---", ""),
    top_genes$probe_id,
    top_genes$gene_symbol
  )

  # ggplot2でプロット作成
  ggplot(res, aes(x = log2FC, y = -log10(padj), color = sig)) +
    geom_point(size = 1, alpha = 0.6) +  # 全点を描画（半透明で重なりを見やすく）

    # 水平の破線: p値カットオフ（-log10(0.05) ≈ 1.3）
    geom_hline(yintercept = -log10(PVAL_CUTOFF), linetype = "dashed",
               color = "black", linewidth = 0.5) +

    # 垂直の破線: FCカットオフ（±log2(1.5)）
    geom_vline(xintercept = c(-log2(FC_CUTOFF), log2(FC_CUTOFF)),
               linetype = "dashed", color = "black", linewidth = 0.5) +

    # 重ならないようにラベルを自動配置（ggrepelパッケージ）
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

# ---- CT6解析 ----
message("--- CT6 解析 ---")
res_ct6 <- compute_results(ko_ct6, wt_ct6)
p_ct6   <- draw_volcano(res_ct6, "CT6")

ggsave("volcano_CT6.pdf", plot = p_ct6, width = 7, height = 6)
write.csv(res_ct6[order(res_ct6$padj), ], "volcano_results_CT6.csv", row.names = FALSE)
message("Saved: volcano_CT6.pdf / volcano_results_CT6.csv")

# ---- CT18解析 ----
message("--- CT18 解析 ---")
res_ct18 <- compute_results(ko_ct18, wt_ct18)
p_ct18   <- draw_volcano(res_ct18, "CT18")

ggsave("volcano_CT18.pdf", plot = p_ct18, width = 7, height = 6)
write.csv(res_ct18[order(res_ct18$padj), ], "volcano_results_CT18.csv", row.names = FALSE)
message("Saved: volcano_CT18.pdf / volcano_results_CT18.csv")

# ---- CT6・CT18を横並びで表示（patchworkパッケージ） ----
library(patchwork)
# + 演算子でプロットを横に並べる（patchwork記法）
combined <- p_ct6 + p_ct18 + plot_layout(ncol = 2)
ggsave("volcano_CT6_CT18_combined.pdf", plot = combined, width = 14, height = 6)
message("Saved: volcano_CT6_CT18_combined.pdf")

print(combined)

# ---- CT6・CT18の網羅的統合リスト ----
# 全プローブのCT6・CT18結果を横に並べてカテゴリ分類する

# 列名に "_CT6" / "_CT18" のサフィックスをつけてマージ
all6  <- res_ct6[,  c("probe_id", "gene_symbol", "log2FC", "pval", "sig")]
all18 <- res_ct18[, c("probe_id", "gene_symbol", "log2FC", "pval", "sig")]
colnames(all6)  <- c("probe_id", "gene_symbol", "log2FC_CT6",  "pval_CT6",  "sig_CT6")
colnames(all18) <- c("probe_id", "gene_symbol", "log2FC_CT18", "pval_CT18", "sig_CT18")

merged <- merge(all6, all18, by = c("probe_id", "gene_symbol"))

# 各解析での方向性フラグ（TRUE/FALSE）
merged$up_CT6    <- merged$sig_CT6  == "Up in KO"
merged$up_CT18   <- merged$sig_CT18 == "Up in KO"
merged$down_CT6  <- merged$sig_CT6  == "Down in KO"
merged$down_CT18 <- merged$sig_CT18 == "Down in KO"

# カテゴリ分類
# Up_both: CT6・CT18両方でUp → 最も信頼性が高い
merged$category <- "NS"
merged$category[merged$up_CT6  & merged$up_CT18 ]  <- "Up_both"
merged$category[merged$down_CT6 & merged$down_CT18] <- "Down_both"
merged$category[merged$up_CT6  & !merged$up_CT18  & !merged$down_CT18] <- "Up_CT6only"
merged$category[merged$up_CT18 & !merged$up_CT6   & !merged$down_CT6 ] <- "Up_CT18only"
merged$category[merged$down_CT6 & !merged$down_CT18 & !merged$up_CT18] <- "Down_CT6only"
merged$category[merged$down_CT18 & !merged$down_CT6 & !merged$up_CT6 ] <- "Down_CT18only"

# CT6とCT18のp値の幾何平均（両時間帯で安定して有意なものほど小さい値になる）
merged$pval_mean <- sqrt(merged$pval_CT6 * merged$pval_CT18)

# 全遺伝子テーブル保存
out_cols <- c("probe_id", "gene_symbol",
              "log2FC_CT6", "pval_CT6", "sig_CT6",
              "log2FC_CT18", "pval_CT18", "sig_CT18",
              "category", "pval_mean")
write.csv(merged[order(merged$pval_mean), out_cols], "all_genes_CT6_CT18.csv", row.names = FALSE)
message(sprintf("Saved: all_genes_CT6_CT18.csv  (%d probes total)", nrow(merged)))

# Up遺伝子リスト（リガンド候補の絞り込み用）
up_genes <- merged[merged$category %in% c("Up_both", "Up_CT6only", "Up_CT18only"), ]
up_genes <- up_genes[order(up_genes$pval_mean), out_cols]
write.csv(up_genes, "up_genes_KO_comprehensive.csv", row.names = FALSE)
message(sprintf("Saved: up_genes_KO_comprehensive.csv  (%d probes)", nrow(up_genes)))

# サマリーをコンソールに表示
message("\n===== Up遺伝子 サマリー =====")
message(sprintf("  CT6・CT18 両方でUp : %d", sum(merged$category == "Up_both")))
message(sprintf("  CT6のみでUp        : %d", sum(merged$category == "Up_CT6only")))
message(sprintf("  CT18のみでUp       : %d", sum(merged$category == "Up_CT18only")))

message("\n--- Up遺伝子 Top30（pval_meanでソート） ---")
print(head(up_genes[, c("gene_symbol", "log2FC_CT6", "log2FC_CT18",
                         "pval_CT6", "pval_CT18", "category")], 30))
