# GO enrichment analysis for Gpr176 KO vs WT differentially expressed genes
# 入力: consistency_summary.csv
# Required packages: clusterProfiler, org.Mm.eg.db, ggplot2, enrichplot

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)

# ---- Parameters ----
INPUT_FILE    <- "consistency_summary.csv"
PVAL_CUTOFF   <- 0.05   # GO解析のp値カットオフ
QVAL_CUTOFF   <- 0.2    # q値（FDR）カットオフ（サンプル少なめなので緩め）
MIN_GENE_SET  <- 5      # 最小遺伝子セットサイズ
TOP_N_PLOT    <- 20     # プロットする上位カテゴリ数

# ---- データ読み込み ----
message("Reading: ", INPUT_FILE)
dat <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# ---- 遺伝子リスト作成 ----
# 2解析以上で一致した遺伝子を対象
up_genes   <- dat$gene_symbol[dat$n_up >= 2 & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"]
down_genes <- dat$gene_symbol[dat$n_dn >= 2 & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"]
all_sig    <- dat$gene_symbol[
  (dat$n_up >= 1 | dat$n_dn >= 1) & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"
]

# 背景遺伝子（全プローブのうちgene symbolがあるもの）
background <- dat$gene_symbol[!is.na(dat$gene_symbol) & dat$gene_symbol != "---"]
background <- unique(background)

message(sprintf("Up遺伝子（2解析以上一致）: %d", length(up_genes)))
message(sprintf("Down遺伝子（2解析以上一致）: %d", length(down_genes)))
message(sprintf("背景遺伝子: %d", length(background)))

# Gene Symbol → Entrez ID 変換
symbol_to_entrez <- function(symbols) {
  ids <- mapIds(org.Mm.eg.db, keys = unique(symbols),
                column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  ids <- ids[!is.na(ids)]
  ids
}

bg_entrez   <- symbol_to_entrez(background)
up_entrez   <- symbol_to_entrez(up_genes)
down_entrez <- symbol_to_entrez(down_genes)
all_entrez  <- symbol_to_entrez(all_sig)

message(sprintf("EntrezID変換: Up=%d, Down=%d, 背景=%d",
                length(up_entrez), length(down_entrez), length(bg_entrez)))

# ---- GO解析関数 ----
run_go <- function(gene_ids, bg_ids, label, ontology = "BP") {
  if (length(gene_ids) < 3) {
    message(sprintf("  [%s] 遺伝子数が少なすぎます（%d個）", label, length(gene_ids)))
    return(NULL)
  }
  tryCatch({
    enrichGO(
      gene          = gene_ids,
      universe      = bg_ids,
      OrgDb         = org.Mm.eg.db,
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = PVAL_CUTOFF,
      qvalueCutoff  = QVAL_CUTOFF,
      minGSSize     = MIN_GENE_SET,
      readable      = TRUE
    )
  }, error = function(e) { message("  エラー: ", e$message); NULL })
}

# ---- 各グループでGO解析（BP: Biological Process） ----
message("\n--- GO解析（Biological Process） ---")

go_up   <- run_go(up_entrez,   bg_entrez, "Up遺伝子")
go_down <- run_go(down_entrez, bg_entrez, "Down遺伝子")
go_all  <- run_go(all_entrez,  bg_entrez, "全変動遺伝子")

# ---- 結果CSV保存 ----
save_go_csv <- function(go_res, filename) {
  if (!is.null(go_res) && nrow(go_res) > 0) {
    write.csv(as.data.frame(go_res), filename, row.names = FALSE)
    message(sprintf("Saved: %s  (%d terms)", filename, nrow(go_res)))
  } else {
    message(sprintf("  有意なGOタームなし → %s は出力しません", filename))
  }
}

save_go_csv(go_up,   "go_results_up.csv")
save_go_csv(go_down, "go_results_down.csv")
save_go_csv(go_all,  "go_results_all.csv")

# ---- プロット関数 ----
plot_go <- function(go_res, title_label, filename_prefix) {
  if (is.null(go_res) || nrow(go_res) == 0) {
    message(sprintf("  [%s] プロットするタームがありません", title_label))
    return(invisible(NULL))
  }

  # dotplot
  p1 <- dotplot(go_res, showCategory = TOP_N_PLOT, font.size = 9) +
    labs(title = sprintf("GO Biological Process: %s", title_label)) +
    theme(plot.title = element_text(size = 11))

  ggsave(sprintf("%s_dotplot.pdf", filename_prefix),
         plot = p1, width = 8, height = 2 + min(TOP_N_PLOT, nrow(go_res)) * 0.28)
  message(sprintf("Saved: %s_dotplot.pdf", filename_prefix))

  # barplot
  p2 <- barplot(go_res, showCategory = TOP_N_PLOT, font.size = 9) +
    labs(title = sprintf("GO Biological Process: %s", title_label)) +
    theme(plot.title = element_text(size = 11))

  ggsave(sprintf("%s_barplot.pdf", filename_prefix),
         plot = p2, width = 8, height = 2 + min(TOP_N_PLOT, nrow(go_res)) * 0.28)
  message(sprintf("Saved: %s_barplot.pdf", filename_prefix))

  # emap (enrichment map) - タームが5個以上の場合
  if (nrow(go_res) >= 5) {
    go_sim <- pairwise_termsim(go_res)
    p3 <- emapplot(go_sim, showCategory = min(30, nrow(go_res))) +
      labs(title = sprintf("GO Term Network: %s", title_label))
    ggsave(sprintf("%s_network.pdf", filename_prefix),
           plot = p3, width = 8, height = 7)
    message(sprintf("Saved: %s_network.pdf", filename_prefix))
  }
}

plot_go(go_up,   "KOでUp（2解析以上）",   "go_up")
plot_go(go_down, "KOでDown（2解析以上）", "go_down")
plot_go(go_all,  "KO変動遺伝子（全体）",  "go_all")

# ---- MF（Molecular Function）とCC（Cellular Component）も実施 ----
for (ont in c("MF", "CC")) {
  message(sprintf("\n--- GO解析（%s） ---", ont))
  go_all_ont <- run_go(all_entrez, bg_entrez, sprintf("全変動遺伝子_%s", ont), ont)
  save_go_csv(go_all_ont, sprintf("go_results_all_%s.csv", tolower(ont)))
  plot_go(go_all_ont, sprintf("全変動遺伝子 [%s]", ont), sprintf("go_all_%s", tolower(ont)))
}

message("\n===== GO解析完了 =====")
message("出力ファイル一覧:")
message("  go_results_up/down/all.csv    ← 結果テーブル（BP）")
message("  go_up/down/all_dotplot.pdf    ← ドットプロット")
message("  go_up/down/all_barplot.pdf    ← バープロット")
message("  go_up/down/all_network.pdf    ← タームネットワーク図")
message("  go_results_all_mf/cc.csv      ← MF・CC結果テーブル")
