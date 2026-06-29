# ============================================================
# summary_consistency.R
# CT6別・CT18別・統合の3解析を横断して「一致性」をまとめる
#
# 【なぜ一致性を見るか】
#   n=2の解析はノイズが多い
#   複数の解析で同じ方向に変動している遺伝子 = 偶然ではなく本物の変化
#   の可能性が高い → 候補遺伝子の優先順位づけに使う
#
# 【入力ファイル】（先に3つのvolcano_plot*.Rを実行しておく）
#   volcano_results_CT6.csv
#   volcano_results_CT18.csv
#   volcano_results_combined.csv
#
# 【出力ファイル】
#   consistency_summary.csv  : 全遺伝子の統合テーブル
#   consistency_plot.pdf     : CT6 vs CT18 log2FC散布図
#   consistency_heatmap.pdf  : 上位遺伝子のlog2FCヒートマップ
# ============================================================

library(ggplot2)
library(ggrepel)

# ---- 各解析結果の読み込み ----
ct6  <- read.csv("volcano_results_CT6.csv",      stringsAsFactors = FALSE)
ct18 <- read.csv("volcano_results_CT18.csv",     stringsAsFactors = FALSE)
comb <- read.csv("volcano_results_combined.csv", stringsAsFactors = FALSE)

# 列名にサフィックスをつけて区別できるようにする
colnames(ct6)[colnames(ct6)   == "padj"]   <- "padj_ct6"
colnames(ct18)[colnames(ct18) == "padj"]   <- "padj_ct18"
colnames(ct6)[colnames(ct6)   == "pval"]   <- "pval_ct6"
colnames(ct18)[colnames(ct18) == "pval"]   <- "pval_ct18"
colnames(comb)[colnames(comb) == "pval"]   <- "pval_comb"
colnames(comb)[colnames(comb) == "padj"]   <- "padj_comb"
colnames(ct6)[colnames(ct6)   == "log2FC"] <- "log2FC_ct6"
colnames(ct18)[colnames(ct18) == "log2FC"] <- "log2FC_ct18"
colnames(comb)[colnames(comb) == "log2FC"] <- "log2FC_comb"
colnames(ct6)[colnames(ct6)   == "sig"]    <- "sig_ct6"
colnames(ct18)[colnames(ct18) == "sig"]    <- "sig_ct18"

# ---- 3つのテーブルをプローブIDと遺伝子名でマージ ----
# merge(..., all = TRUE): 片方にしかない行もNAで残す（outer join）
m <- merge(
  ct6[,  c("probe_id","gene_symbol","log2FC_ct6","pval_ct6","sig_ct6")],
  ct18[, c("probe_id","gene_symbol","log2FC_ct18","pval_ct18","sig_ct18")],
  by = c("probe_id","gene_symbol"), all = TRUE
)
m <- merge(
  m,
  comb[, c("probe_id","gene_symbol","log2FC_comb","pval_comb","sig_raw")],
  by = c("probe_id","gene_symbol"), all = TRUE
)
colnames(m)[colnames(m) == "sig_raw"] <- "sig_comb"

# ---- 各解析でUp/Downかどうかのフラグ（TRUE/FALSE） ----
m$up_ct6   <- !is.na(m$sig_ct6)  & m$sig_ct6  == "Up in KO"
m$up_ct18  <- !is.na(m$sig_ct18) & m$sig_ct18 == "Up in KO"
m$up_comb  <- !is.na(m$sig_comb) & m$sig_comb == "Up in KO"
m$dn_ct6   <- !is.na(m$sig_ct6)  & m$sig_ct6  == "Down in KO"
m$dn_ct18  <- !is.na(m$sig_ct18) & m$sig_ct18 == "Down in KO"
m$dn_comb  <- !is.na(m$sig_comb) & m$sig_comb == "Down in KO"

# 何解析で有意か（0〜3の整数）
# as.integer(TRUE)=1, as.integer(FALSE)=0 を利用して合計する
m$n_up <- as.integer(m$up_ct6) + as.integer(m$up_ct18) + as.integer(m$up_comb)
m$n_dn <- as.integer(m$dn_ct6) + as.integer(m$dn_ct18) + as.integer(m$dn_comb)

# 一致カテゴリの割り当て（n_upやn_dnが大きいほど信頼性が高い）
m$consistency <- "NS"
m$consistency[m$n_up == 3] <- "Up_全解析一致"
m$consistency[m$n_up == 2] <- "Up_2解析一致"
m$consistency[m$n_up == 1] <- "Up_1解析のみ"
m$consistency[m$n_dn == 3] <- "Down_全解析一致"
m$consistency[m$n_dn == 2] <- "Down_2解析一致"
m$consistency[m$n_dn == 1] <- "Down_1解析のみ"

# 3解析のp値の幾何平均（全解析で安定して低いp値ほど小さい値になる）
# NAはp値=1（非有意）として扱う
m$pval_mean <- exp(rowMeans(log(cbind(
  ifelse(is.na(m$pval_ct6),  1, m$pval_ct6),
  ifelse(is.na(m$pval_ct18), 1, m$pval_ct18),
  ifelse(is.na(m$pval_comb), 1, m$pval_comb)
)), na.rm = TRUE))

# ---- 統合テーブルのCSV保存 ----
out_cols <- c("probe_id","gene_symbol",
              "log2FC_ct6","pval_ct6","sig_ct6",
              "log2FC_ct18","pval_ct18","sig_ct18",
              "log2FC_comb","pval_comb","sig_comb",
              "n_up","n_dn","consistency","pval_mean")
m_sorted <- m[order(m$pval_mean), out_cols]

write.csv(m_sorted, "consistency_summary.csv", row.names = FALSE)
message("Saved: consistency_summary.csv")

# サマリーをコンソールに表示
message("\n===== 一致性サマリー =====")
for (cat in c("Up_全解析一致","Up_2解析一致","Down_全解析一致","Down_2解析一致")) {
  genes <- m$gene_symbol[m$consistency == cat & !is.na(m$gene_symbol)]
  genes <- genes[genes != "---"]
  message(sprintf("  %s (%d件): %s",
                  cat, sum(m$consistency == cat),
                  paste(head(genes, 10), collapse = ", ")))
}

# ---- 図1: CT6 vs CT18 log2FC 散布図 ----
# X軸: CT6でのlog2FC、Y軸: CT18でのlog2FC
# 対角線（傾き1）上の点 = CT6とCT18で同じ方向・同じ大きさの変動
# 点の色 = 一致性カテゴリ

highlight <- m[m$consistency %in% c("Up_全解析一致","Up_2解析一致",
                                     "Down_全解析一致","Down_2解析一致"), ]
highlight <- highlight[order(highlight$pval_mean), ]
highlight$label <- ifelse(
  is.na(highlight$gene_symbol) | highlight$gene_symbol %in% c("---",""),
  highlight$probe_id, highlight$gene_symbol
)

color_map <- c(
  "Up_全解析一致"   = "#B22222",  # 濃い赤: 最も信頼性高い
  "Up_2解析一致"    = "#F08080",  # 薄い赤
  "Down_全解析一致" = "#00008B",  # 濃い青
  "Down_2解析一致"  = "#6495ED",  # 薄い青
  "NS"              = "grey80"    # 灰色: 非有意
)
size_map <- c(
  "Up_全解析一致" = 3, "Up_2解析一致" = 2,
  "Down_全解析一致" = 3, "Down_2解析一致" = 2, "NS" = 0.8
)

m$consistency <- factor(m$consistency,
  levels = c("Up_全解析一致","Up_2解析一致","Down_全解析一致","Down_2解析一致","NS"))

p <- ggplot(m, aes(x = log2FC_ct6, y = log2FC_ct18,
                   color = consistency, size = consistency)) +
  geom_point(alpha = 0.5) +
  # X=0, Y=0 の基準線
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  # 対角線: この線上 = CT6とCT18で同じ変動
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.4) +
  geom_text_repel(
    data = head(highlight, 40),
    aes(x = log2FC_ct6, y = log2FC_ct18, label = label),
    color = "black", size = 2.3, max.overlaps = 40,
    box.padding = 0.3, point.padding = 0.2, inherit.aes = FALSE
  ) +
  scale_color_manual(values = color_map, name = "一致性") +
  scale_size_manual(values = size_map, name = "一致性") +
  labs(
    title    = "CT6 vs CT18 log2FC 比較（KO/WT）",
    subtitle = "対角線上 = CT6・CT18で同方向変動  |  色 = 解析間の一致性",
    x = expression(log[2]~"FC (CT6: KO/WT)"),
    y = expression(log[2]~"FC (CT18: KO/WT)")
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave("consistency_plot.pdf", plot = p, width = 8, height = 7)
message("Saved: consistency_plot.pdf")
print(p)

# ---- 図2: 上位遺伝子のlog2FCヒートマップ ----
# 2解析以上で一致した遺伝子について、CT6/CT18/統合の3列でlog2FCを色表示
# 赤 = Up（KOで増加）、青 = Down（KOで減少）

top30 <- head(highlight, 30)
top30 <- top30[!is.na(top30$log2FC_ct6) & !is.na(top30$log2FC_ct18), ]
top30$label <- ifelse(
  is.na(top30$gene_symbol) | top30$gene_symbol %in% c("---",""),
  top30$probe_id, top30$gene_symbol)
# factor化して表示順を固定（reverseにするとY軸の上からpval_mean順になる）
top30$label <- factor(top30$label, levels = rev(top30$label))

# wide形式 → long形式に変換（ggplot2はlong形式が基本）
long <- rbind(
  data.frame(label=top30$label, log2FC=top30$log2FC_ct6,  analysis="CT6"),
  data.frame(label=top30$label, log2FC=top30$log2FC_ct18, analysis="CT18"),
  data.frame(label=top30$label, log2FC=top30$log2FC_comb, analysis="統合(n=4)")
)
long$analysis <- factor(long$analysis, levels = c("CT6","CT18","統合(n=4)"))

p2 <- ggplot(long, aes(x = analysis, y = label, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # 各タイルにlog2FC値をテキスト表示
  geom_text(aes(label = sprintf("%.2f", log2FC)), size = 2.5, color = "black") +
  # scale_fill_gradient2: 0を白、正を赤、負を青で表現
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0, name = "log2FC") +
  labs(
    title    = "2解析以上で一致した遺伝子のlog2FC",
    subtitle = "CT6 / CT18 / 統合(n=4) 比較",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 9),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave("consistency_heatmap.pdf", plot = p2, width = 6,
       height = 2 + nrow(top30) * 0.32)
message("Saved: consistency_heatmap.pdf")
print(p2)
