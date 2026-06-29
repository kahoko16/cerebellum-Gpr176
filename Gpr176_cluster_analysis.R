# ============================================================
# Gpr176_cluster_analysis.R
#
# 論文: Kozareva et al. 2021, Nature
# "A transcriptomic atlas of mouse cerebellar cortex
#  comprehensively defines cell types"
# GEO: GSE165805
#
# 目的: Gpr176陽性クラスターにGnaz（Gz）とRGS16が発現するかを検証
#
# 解析の流れ:
#   1. GEOからメタデータ（クラスターラベル）と発現行列を取得
#   2. Gpr176の発現でクラスターを定義（陽性 vs 陰性）
#   3. 各クラスター内でのGnaz・RGS16の発現量を可視化
# ============================================================

# ---- パッケージ ----
required_pkgs <- c("GEOquery", "Seurat", "ggplot2", "dplyr", "Matrix",
                   "patchwork", "viridis", "stringr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("GEOquery", "Seurat")) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(pkg, character.only = TRUE)
}

# ---- 設定 ----
GEO_ID       <- "GSE165805"
GENES_OF_INT <- c("Gpr176", "Gnaz", "Rgs16")  # Gz = Gnaz
OUT_DIR      <- "Gpr176_cluster_results"
dir.create(OUT_DIR, showWarnings = FALSE)

# Gpr176「陽性」の閾値（正規化カウント > この値なら陽性）
GPR176_THRESHOLD <- 0  # 0より大きければ陽性（検出可）

# ---- Step 1: GEOからデータ取得 ----
message("\n=== Step 1: GEOデータ取得 ===")
# GEOqueryでシリーズ情報を取得（SupplementaryFilesのURLを確認するため）
gse <- GEOquery::getGEO(GEO_ID, GSEMatrix = FALSE)

# GSE165805はSeuratオブジェクト形式で提供されている場合があるため
# まずGSEMatrixで試みる
gse_matrix <- GEOquery::getGEO(GEO_ID, GSEMatrix = TRUE)
message("GSEMatrix取得完了")

# ---- Step 2: 補足ファイルのダウンロード ----
# GSE165805は複数のSupplementaryファイルを持つ
# count matrix + cell metadata の取得を試みる
message("\n=== Step 2: 補足ファイルのダウンロード ===")

supp_dir <- file.path(OUT_DIR, "supp_files")
dir.create(supp_dir, showWarnings = FALSE)

# GEOquery::getGEOSuppFiles でダウンロード
supp_files <- GEOquery::getGEOSuppFiles(GEO_ID, baseDir = supp_dir, fetch_files = TRUE)
message("補足ファイルリスト:")
print(rownames(supp_files))

# ---- Step 3: ファイル特定とロード ----
message("\n=== Step 3: ファイルのロードと解析 ===")

supp_paths <- rownames(supp_files)

# ファイル種別を特定する関数
find_file <- function(paths, patterns) {
  for (pat in patterns) {
    m <- grep(pat, paths, ignore.case = TRUE, value = TRUE)
    if (length(m) > 0) return(m[1])
  }
  return(NULL)
}

# Seuratオブジェクト(.rds)があればそれを使う
rds_file <- find_file(supp_paths, c("\\.rds$", "\\.RDS$", "\\.Rds$"))

if (!is.null(rds_file)) {
  # .gz圧縮なら解凍
  if (grepl("\\.gz$", rds_file)) {
    rds_ungz <- sub("\\.gz$", "", rds_file)
    if (!file.exists(rds_ungz)) R.utils::gunzip(rds_file, destname = rds_ungz, remove = FALSE)
    rds_file <- rds_ungz
  }
  message("Seuratオブジェクトを読み込み中: ", rds_file)
  seurat_obj <- readRDS(rds_file)

} else {
  # matrix.mtx / barcodes.tsv / features.tsv の組み合わせを探す
  mtx_file      <- find_file(supp_paths, c("matrix.*\\.mtx", "count.*\\.mtx"))
  barcode_file  <- find_file(supp_paths, c("barcode", "cell"))
  feature_file  <- find_file(supp_paths, c("feature", "gene"))
  meta_file     <- find_file(supp_paths, c("meta", "cluster", "annotation", "label"))

  # .gz解凍
  decompress <- function(f) {
    if (!is.null(f) && grepl("\\.gz$", f)) {
      out <- sub("\\.gz$", "", f)
      if (!file.exists(out)) R.utils::gunzip(f, destname = out, remove = FALSE)
      return(out)
    }
    f
  }
  mtx_file     <- decompress(mtx_file)
  barcode_file <- decompress(barcode_file)
  feature_file <- decompress(feature_file)
  meta_file    <- decompress(meta_file)

  if (!is.null(mtx_file) && !is.null(barcode_file) && !is.null(feature_file)) {
    message("10x形式matrixを読み込み中...")
    counts <- Seurat::ReadMtx(mtx = mtx_file,
                              cells = barcode_file,
                              features = feature_file)
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts, project = GEO_ID)

    if (!is.null(meta_file)) {
      meta_df <- read.csv(meta_file, row.names = 1, check.names = FALSE)
      # cellIDが一致するものだけ追加
      common_cells <- intersect(colnames(seurat_obj), rownames(meta_df))
      seurat_obj <- seurat_obj[, common_cells]
      seurat_obj <- Seurat::AddMetaData(seurat_obj, meta_df[common_cells, , drop = FALSE])
    }

  } else {
    # TSV/CSV形式のカウントマトリクスを探す
    count_csv <- find_file(supp_paths, c("count.*\\.csv", "expr.*\\.csv",
                                         "count.*\\.tsv", "expr.*\\.tsv"))
    count_csv <- decompress(count_csv)
    if (!is.null(count_csv)) {
      message("CSV/TSV形式のカウントマトリクスを読み込み中...")
      sep_char <- ifelse(grepl("\\.tsv", count_csv), "\t", ",")
      counts_df <- read.table(count_csv, sep = sep_char, header = TRUE,
                              row.names = 1, check.names = FALSE)
      counts_mat <- as.matrix(counts_df)
      seurat_obj <- Seurat::CreateSeuratObject(counts = counts_mat, project = GEO_ID)
    } else {
      stop("認識できる発現行列ファイルが見つかりませんでした。\n",
           "ダウンロードされたファイル:\n",
           paste(supp_paths, collapse = "\n"))
    }
  }
}

message("Seuratオブジェクト概要:")
print(seurat_obj)

# ---- Step 4: 遺伝子の存在確認 ----
message("\n=== Step 4: 対象遺伝子の確認 ===")

all_genes <- rownames(seurat_obj)

# 大文字小文字を考慮した検索
find_gene <- function(gene, all_genes) {
  # 完全一致
  if (gene %in% all_genes) return(gene)
  # 大文字小文字無視
  m <- all_genes[tolower(all_genes) == tolower(gene)]
  if (length(m) > 0) return(m[1])
  # 部分一致
  m <- grep(paste0("^", gene, "$"), all_genes, value = TRUE, ignore.case = TRUE)
  if (length(m) > 0) return(m[1])
  return(NA)
}

gene_map <- setNames(sapply(GENES_OF_INT, find_gene, all_genes = all_genes), GENES_OF_INT)
message("遺伝子名マッピング:")
print(gene_map)

missing <- gene_map[is.na(gene_map)]
if (length(missing) > 0) {
  message("⚠ 以下の遺伝子が発現行列に見つかりません: ", paste(names(missing), collapse = ", "))
}

found_genes <- gene_map[!is.na(gene_map)]
if (length(found_genes) == 0) stop("解析対象遺伝子が1つも見つかりません")

# ---- Step 5: 正規化 ----
message("\n=== Step 5: データ正規化 ===")
seurat_obj <- Seurat::NormalizeData(seurat_obj,
                                    normalization.method = "LogNormalize",
                                    scale.factor = 10000,
                                    verbose = FALSE)
message("LogNormalize完了（log(CP10K + 1)）")

# ---- Step 6: Gpr176陽性細胞の同定 ----
message("\n=== Step 6: Gpr176陽性クラスターの同定 ===")

gpr176_gene <- found_genes["Gpr176"]
if (is.na(gpr176_gene)) {
  stop("Gpr176が発現行列に見つからないため解析を中断します")
}

# 正規化済み発現量を取得
norm_data  <- Seurat::GetAssayData(seurat_obj, slot = "data")
gpr176_exp <- norm_data[gpr176_gene, ]

# Gpr176陽性フラグ
seurat_obj$Gpr176_positive <- gpr176_exp > GPR176_THRESHOLD
seurat_obj$Gpr176_expr     <- as.numeric(gpr176_exp)

n_pos <- sum(seurat_obj$Gpr176_positive)
n_tot <- ncol(seurat_obj)
message(sprintf("Gpr176陽性細胞: %d / %d (%.1f%%)", n_pos, n_tot, 100 * n_pos / n_tot))

# クラスターアノテーションが存在すればクラスター別に集計
cluster_cols <- c("cluster", "ClusterID", "cell_type", "CellType", "seurat_clusters",
                  "leiden", "louvain", "Cluster", "cluster_label", "subtype")
cluster_col  <- cluster_cols[cluster_cols %in% colnames(seurat_obj@meta.data)][1]

if (!is.na(cluster_col)) {
  message("クラスター列を使用: ", cluster_col)
  seurat_obj$cluster_id <- seurat_obj@meta.data[[cluster_col]]

  # クラスター別Gpr176陽性率
  cluster_summary <- seurat_obj@meta.data %>%
    dplyr::group_by(cluster_id) %>%
    dplyr::summarise(
      n_cells       = dplyr::n(),
      n_gpr176_pos  = sum(Gpr176_positive),
      pct_gpr176    = 100 * mean(Gpr176_positive),
      mean_gpr176   = mean(Gpr176_expr),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(pct_gpr176))

  write.csv(cluster_summary,
            file.path(OUT_DIR, "cluster_Gpr176_summary.csv"),
            row.names = FALSE)
  message("クラスター別Gpr176陽性率 (上位10):")
  print(head(cluster_summary, 10))

  # Gpr176陽性率 > 10% のクラスターを「Gpr176陽性クラスター」とする
  gpr176_pos_clusters <- cluster_summary %>%
    dplyr::filter(pct_gpr176 > 10) %>%
    dplyr::pull(cluster_id)

  if (length(gpr176_pos_clusters) == 0) {
    # 閾値を下げて上位5クラスターを採用
    gpr176_pos_clusters <- cluster_summary$cluster_id[1:min(5, nrow(cluster_summary))]
    message("⚠ 陽性率10%超クラスターなし。上位5クラスターを使用: ",
            paste(gpr176_pos_clusters, collapse = ", "))
  } else {
    message("Gpr176陽性クラスター (陽性率>10%): ",
            paste(gpr176_pos_clusters, collapse = ", "))
  }

  seurat_obj$Gpr176_cluster <- ifelse(
    seurat_obj$cluster_id %in% gpr176_pos_clusters, "Gpr176_positive", "Gpr176_negative"
  )

} else {
  message("クラスターアノテーションが見つかりません。細胞単位で陽性/陰性を使用します")
  seurat_obj$cluster_id    <- NA
  seurat_obj$Gpr176_cluster <- ifelse(seurat_obj$Gpr176_positive, "Gpr176_positive", "Gpr176_negative")
}

# ---- Step 7: GnazとRGS16の発現量を抽出 ----
message("\n=== Step 7: GnazおよびRGS16の発現量解析 ===")

target_genes <- found_genes[names(found_genes) %in% c("Gnaz", "Rgs16")]

for (gname in names(target_genes)) {
  actual_name <- target_genes[gname]
  seurat_obj@meta.data[[paste0(gname, "_expr")]] <- as.numeric(norm_data[actual_name, ])
}

# 統計サマリー（Gpr176陽性 vs 陰性）
result_rows <- list()
for (gname in names(target_genes)) {
  expr_col <- paste0(gname, "_expr")
  pos_vals <- seurat_obj@meta.data[seurat_obj$Gpr176_cluster == "Gpr176_positive", expr_col]
  neg_vals <- seurat_obj@meta.data[seurat_obj$Gpr176_cluster == "Gpr176_negative", expr_col]

  wt <- wilcox.test(pos_vals, neg_vals, alternative = "greater")
  pct_pos <- 100 * mean(pos_vals > 0)
  pct_neg <- 100 * mean(neg_vals > 0)

  result_rows[[gname]] <- data.frame(
    gene         = gname,
    actual_name  = target_genes[gname],
    mean_Gpr176pos = round(mean(pos_vals), 4),
    mean_Gpr176neg = round(mean(neg_vals), 4),
    pct_expressed_pos = round(pct_pos, 2),
    pct_expressed_neg = round(pct_neg, 2),
    log2FC       = round(log2((mean(pos_vals) + 1e-6) / (mean(neg_vals) + 1e-6)), 4),
    wilcox_pval  = signif(wt$p.value, 4),
    stringsAsFactors = FALSE
  )

  message(sprintf("\n[%s (%s)] Gpr176+ クラスター vs Gpr176- クラスター", gname, actual_name))
  message(sprintf("  平均発現量: %.4f vs %.4f", mean(pos_vals), mean(neg_vals)))
  message(sprintf("  発現細胞率: %.1f%% vs %.1f%%", pct_pos, pct_neg))
  message(sprintf("  log2FC: %.4f  Wilcoxon p値: %s", log2((mean(pos_vals)+1e-6)/(mean(neg_vals)+1e-6)), signif(wt$p.value,4)))
}

result_df <- do.call(rbind, result_rows)
write.csv(result_df,
          file.path(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"),
          row.names = FALSE)
message("\n結果をCSVに保存: ", file.path(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"))
print(result_df)

# ---- Step 8: 可視化 ----
message("\n=== Step 8: 可視化 ===")

# 8-A: バイオリンプロット（Gpr176陽性 vs 陰性クラスターでの各遺伝子発現）
plot_list <- list()

all_target <- c("Gpr176", names(target_genes))
for (gname in all_target) {
  actual <- if (gname == "Gpr176") gpr176_gene else target_genes[gname]
  if (is.na(actual)) next

  expr_vals <- as.numeric(norm_data[actual, ])
  plot_df <- data.frame(
    expr    = expr_vals,
    group   = seurat_obj$Gpr176_cluster,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c("Gpr176_positive" = "#E64B35",
                                          "Gpr176_negative" = "#4DBBD5")) +
    ggplot2::labs(title = sprintf("%s (%s)", gname, actual),
                  x = NULL, y = "log(CP10K + 1)",
                  subtitle = sprintf("Gpr176+ vs Gpr176- clusters")) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold"))

  plot_list[[gname]] <- p
}

combined_violin <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.pdf"),
                combined_violin, width = 4 * length(plot_list), height = 5)
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.png"),
                combined_violin, width = 4 * length(plot_list), height = 5, dpi = 150)
message("バイオリンプロット保存完了")

# 8-B: クラスター別ドットプロット（クラスターアノテーションがある場合）
if (!is.na(cluster_col) && length(gpr176_pos_clusters) > 0) {

  dot_genes <- c(gpr176_gene, unname(target_genes))
  dot_genes <- dot_genes[!is.na(dot_genes)]

  # Gpr176陽性クラスター + 発現量上位5陰性クラスターに絞る
  top_neg <- cluster_summary %>%
    dplyr::filter(!cluster_id %in% gpr176_pos_clusters) %>%
    dplyr::slice_head(n = 5) %>%
    dplyr::pull(cluster_id)

  show_clusters <- c(as.character(gpr176_pos_clusters), as.character(top_neg))
  cells_to_show <- colnames(seurat_obj)[seurat_obj$cluster_id %in% show_clusters]

  if (length(cells_to_show) > 1) {
    sub_obj <- seurat_obj[, cells_to_show]
    Idents(sub_obj) <- sub_obj$cluster_id

    dp <- Seurat::DotPlot(sub_obj, features = dot_genes) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Gpr176 / Gnaz / Rgs16 発現（上位クラスター）",
                    x = "遺伝子", y = "クラスター") +
      ggplot2::theme_classic(base_size = 11) +
      viridis::scale_color_viridis(option = "plasma")

    ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.pdf"),
                    dp, width = 10, height = 5)
    ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.png"),
                    dp, width = 10, height = 5, dpi = 150)
    message("ドットプロット保存完了")
  }
}

# 8-C: Gpr176陽性クラスターにおけるGnaz・RGS16の発現率の棒グラフ
if (nrow(result_df) > 0) {
  bar_df <- result_df %>%
    dplyr::select(gene, pct_expressed_pos, pct_expressed_neg) %>%
    tidyr::pivot_longer(cols = c(pct_expressed_pos, pct_expressed_neg),
                        names_to = "group", values_to = "pct") %>%
    dplyr::mutate(group = dplyr::recode(group,
                                        "pct_expressed_pos" = "Gpr176+ cluster",
                                        "pct_expressed_neg" = "Gpr176- cluster"))

  bar_p <- ggplot2::ggplot(bar_df, ggplot2::aes(x = gene, y = pct, fill = group)) +
    ggplot2::geom_col(position = "dodge", width = 0.6) +
    ggplot2::scale_fill_manual(values = c("Gpr176+ cluster" = "#E64B35",
                                          "Gpr176- cluster" = "#4DBBD5")) +
    ggplot2::labs(title = "Gpr176陽性クラスターにおける発現細胞率",
                  x = "遺伝子", y = "発現細胞率 (%)", fill = NULL) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(legend.position = "top")

  ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.pdf"),
                  bar_p, width = 6, height = 5)
  ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.png"),
                  bar_p, width = 6, height = 5, dpi = 150)
  message("棒グラフ保存完了")
}

# ---- 最終サマリー ----
message("\n========================================")
message("         解析完了サマリー")
message("========================================")
message(sprintf("データセット : %s (Kozareva et al. 2021)", GEO_ID))
message(sprintf("総細胞数     : %d", ncol(seurat_obj)))
message(sprintf("Gpr176陽性   : %d 細胞 (%.1f%%)", n_pos, 100 * n_pos / n_tot))
message("")
message("Gpr176陽性クラスターでの発現:")
for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  sig <- if (r$wilcox_pval < 0.001) "***" else if (r$wilcox_pval < 0.01) "**" else if (r$wilcox_pval < 0.05) "*" else "ns"
  message(sprintf("  %s: 発現率 %.1f%% (vs %.1f%%), log2FC=%.2f, p=%s %s",
                  r$gene, r$pct_expressed_pos, r$pct_expressed_neg,
                  r$log2FC, r$wilcox_pval, sig))
}
message("")
message(sprintf("出力フォルダ : %s/", OUT_DIR))
message("  - cluster_Gpr176_summary.csv")
message("  - Gnaz_RGS16_in_Gpr176clusters.csv")
message("  - violin_Gnaz_RGS16_Gpr176clusters.pdf/png")
message("  - dotplot_Gpr176_Gnaz_RGS16.pdf/png")
message("  - barplot_expression_pct.pdf/png")
message("========================================")
