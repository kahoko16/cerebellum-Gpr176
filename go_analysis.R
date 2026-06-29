# ============================================================
# go_analysis.R
# 変動遺伝子のGene Ontology（GO）エンリッチメント解析
#
# 【GO解析とは】
#   変動遺伝子リストが、どんな生物学的カテゴリ（経路・機能・局在）に
#   偏っているかを統計的に検定する
#
# 【GOの3カテゴリ】
#   BP (Biological Process)  : 生物学的プロセス（シグナル伝達、代謝など）
#   MF (Molecular Function)  : 分子機能（受容体活性、酵素活性など）
#   CC (Cellular Component)  : 細胞内局在（核、細胞膜、繊毛など）
#
# 【エンリッチメントとは】
#   「変動遺伝子リスト」の中にあるGOタームの割合が
#   「背景遺伝子全体」の中の割合より統計的に高い（偏っている）かを検定
#   → 偶然以上に多く含まれるGOタームを「エンリッチされている」という
#
# 【入力ファイル】（先にsummary_consistency.Rを実行しておく）
#   consistency_summary.csv
# ============================================================

library(clusterProfiler)  # GO/KEGG解析の定番パッケージ
library(org.Mm.eg.db)     # マウス(Mus musculus)のアノテーションデータベース
library(ggplot2)
library(enrichplot)        # GO解析結果の可視化

# ---- パラメータ ----
INPUT_FILE    <- "consistency_summary.csv"
PVAL_CUTOFF   <- 0.05  # GO解析のp値カットオフ
QVAL_CUTOFF   <- 0.2   # q値（FDR）カットオフ
               # ※サンプル数が少なく遺伝子数も少ないため、通常より緩めに設定
MIN_GENE_SET  <- 5     # GOタームに含まれる最小遺伝子数（少なすぎるタームを除外）
TOP_N_PLOT    <- 20    # プロットする上位タームの数

# ---- データ読み込みと遺伝子リストの作成 ----
message("Reading: ", INPUT_FILE)
dat <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)

# 2解析以上で一致した遺伝子を「信頼性の高い変動遺伝子」として使用
up_genes   <- dat$gene_symbol[dat$n_up >= 2 & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"]
down_genes <- dat$gene_symbol[dat$n_dn >= 2 & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"]

# 1解析以上で有意な遺伝子（網羅的なリスト）
all_sig <- dat$gene_symbol[
  (dat$n_up >= 1 | dat$n_dn >= 1) & !is.na(dat$gene_symbol) & dat$gene_symbol != "---"
]

# 背景遺伝子: 解析に使った全プローブのうちgene symbolが付いているもの
# GOエンリッチメントは「変動遺伝子/背景遺伝子」の比率で検定するため、
# 適切な背景遺伝子を設定することが重要
background <- unique(dat$gene_symbol[!is.na(dat$gene_symbol) & dat$gene_symbol != "---"])

message(sprintf("Up遺伝子（2解析以上一致）: %d", length(up_genes)))
message(sprintf("Down遺伝子（2解析以上一致）: %d", length(down_genes)))
message(sprintf("背景遺伝子: %d", length(background)))

# ---- Gene Symbol → Entrez ID 変換 ----
# clusterProfilerはEntrez IDを使うため変換が必要
# org.Mm.eg.db: マウスの遺伝子アノテーションデータベース
symbol_to_entrez <- function(symbols) {
  ids <- mapIds(org.Mm.eg.db,
                keys     = unique(symbols),
                column   = "ENTREZID",  # 変換先
                keytype  = "SYMBOL",    # 変換元
                multiVals = "first")    # 1対多の場合は最初の1つを使用
  ids <- ids[!is.na(ids)]  # 変換できなかったものを除外
  ids
}

bg_entrez   <- symbol_to_entrez(background)
up_entrez   <- symbol_to_entrez(up_genes)
down_entrez <- symbol_to_entrez(down_genes)
all_entrez  <- symbol_to_entrez(all_sig)

message(sprintf("EntrezID変換: Up=%d, Down=%d, 背景=%d",
                length(up_entrez), length(down_entrez), length(bg_entrez)))

# ---- GO解析を実行する関数 ----
run_go <- function(gene_ids, bg_ids, label, ontology = "BP") {
  if (length(gene_ids) < 3) {
    message(sprintf("  [%s] 遺伝子数が少なすぎます（%d個）", label, length(gene_ids)))
    return(NULL)
  }
  tryCatch({
    enrichGO(
      gene          = gene_ids,       # 対象遺伝子（Entrez ID）
      universe      = bg_ids,         # 背景遺伝子（Entrez ID）
      OrgDb         = org.Mm.eg.db,  # アノテーションDB（マウス）
      ont           = ontology,       # "BP", "MF", "CC" のいずれか
      pAdjustMethod = "BH",           # 多重比較補正（BH法）
      pvalueCutoff  = PVAL_CUTOFF,
      qvalueCutoff  = QVAL_CUTOFF,
      minGSSize     = MIN_GENE_SET,
      readable      = TRUE            # Entrez ID → Gene Symbolに戻して表示
    )
  }, error = function(e) { message("  エラー: ", e$message); NULL })
}

# ---- GO解析（BP: Biological Process）を各グループで実施 ----
message("\n--- GO解析（Biological Process） ---")
go_up   <- run_go(up_entrez,   bg_entrez, "Up遺伝子")
go_down <- run_go(down_entrez, bg_entrez, "Down遺伝子")
go_all  <- run_go(all_entrez,  bg_entrez, "全変動遺伝子")

# ---- 結果CSVを保存する関数 ----
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

# ---- 可視化する関数 ----
plot_go <- function(go_res, title_label, filename_prefix) {
  if (is.null(go_res) || nrow(go_res) == 0) {
    message(sprintf("  [%s] プロットするタームがありません", title_label))
    return(invisible(NULL))
  }

  # dotplot: X軸=遺伝子比率、点の大きさ=遺伝子数、色=p値
  # → エンリッチされたタームの「強さ」と「有意性」を同時に表現できる
  p1 <- dotplot(go_res, showCategory = TOP_N_PLOT, font.size = 9) +
    labs(title = sprintf("GO Biological Process: %s", title_label)) +
    theme(plot.title = element_text(size = 11))
  ggsave(sprintf("%s_dotplot.pdf", filename_prefix),
         plot = p1, width = 8, height = 2 + min(TOP_N_PLOT, nrow(go_res)) * 0.28)
  message(sprintf("Saved: %s_dotplot.pdf", filename_prefix))

  # barplot: X軸=遺伝子数、色=p値 → シンプルで見やすい
  p2 <- barplot(go_res, showCategory = TOP_N_PLOT, font.size = 9) +
    labs(title = sprintf("GO Biological Process: %s", title_label)) +
    theme(plot.title = element_text(size = 11))
  ggsave(sprintf("%s_barplot.pdf", filename_prefix),
         plot = p2, width = 8, height = 2 + min(TOP_N_PLOT, nrow(go_res)) * 0.28)
  message(sprintf("Saved: %s_barplot.pdf", filename_prefix))

  # emapplot（エンリッチメントマップ）: タームを点、類似タームを線で結ぶネットワーク図
  # → 関連するGOタームがクラスターとしてまとまって見える
  if (nrow(go_res) >= 5) {
    go_sim <- pairwise_termsim(go_res)  # タームの類似度を計算
    p3 <- emapplot(go_sim, showCategory = min(30, nrow(go_res))) +
      labs(title = sprintf("GO Term Network: %s", title_label))
    ggsave(sprintf("%s_network.pdf", filename_prefix), plot = p3, width = 8, height = 7)
    message(sprintf("Saved: %s_network.pdf", filename_prefix))
  }
}

plot_go(go_up,   "KOでUp（2解析以上）",   "go_up")
plot_go(go_down, "KOでDown（2解析以上）", "go_down")
plot_go(go_all,  "KO変動遺伝子（全体）",  "go_all")

# ---- MF（Molecular Function）とCC（Cellular Component）も解析 ----
# MF: 遺伝子産物の分子レベルの機能（受容体、酵素、転写因子など）
# CC: 遺伝子産物がどこにあるか（核、細胞膜、一次繊毛など）
for (ont in c("MF", "CC")) {
  message(sprintf("\n--- GO解析（%s） ---", ont))
  go_all_ont <- run_go(all_entrez, bg_entrez, sprintf("全変動遺伝子_%s", ont), ont)
  save_go_csv(go_all_ont, sprintf("go_results_all_%s.csv", tolower(ont)))
  plot_go(go_all_ont, sprintf("全変動遺伝子 [%s]", ont), sprintf("go_all_%s", tolower(ont)))
}

message("\n===== GO解析完了 =====")
message("出力ファイル一覧:")
message("  go_results_up/down/all.csv      ← 結果テーブル（BP）")
message("  go_up/down/all_dotplot.pdf      ← ドットプロット")
message("  go_up/down/all_barplot.pdf      ← バープロット")
message("  go_up/down/all_network.pdf      ← タームネットワーク図")
message("  go_results_all_mf/cc.csv        ← MF・CC結果テーブル")
