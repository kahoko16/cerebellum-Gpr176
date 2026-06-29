# ============================================================
# Gpr176_cluster_analysis.R
#
# 論文: Kozareva et al. 2021, Nature
# "A transcriptomic atlas of mouse cerebellar cortex
#  comprehensively defines cell types"
# GEO: GSE165371
#
# 目的: プルキンエ細胞クラスターにGpr176が発現するか、
#       またGpr176陽性クラスターにGnaz（Gz）とRGS16が
#       共発現するかを検証
#
# 事前準備:
#   GSE165371_cb_adult_mouse.tar.gz をこのスクリプトと
#   同じフォルダに置いてから実行してください。
# ============================================================

# ---- パッケージ ----
required_cran <- c("ggplot2", "dplyr", "patchwork", "viridis",
                   "tidyr", "Matrix", "stringr")
required_bioc <- c("GEOquery")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}
# Seurat は任意（あれば使う）
has_seurat <- requireNamespace("Seurat", quietly = TRUE)
if (has_seurat) library(Seurat)

# ---- 設定 ----
TAR_FILE         <- "GSE165371_cb_adult_mouse.tar.gz"  # ダウンロードしたファイル
GENES_OF_INT     <- c("Gpr176", "Gnaz", "Rgs16")
OUT_DIR          <- "Gpr176_cluster_results"
GPR176_THRESHOLD <- 0      # log-norm > 0 で陽性
POS_PCT_THR      <- 10     # クラスター内Gpr176陽性率(%)がこれ以上を陽性クラスターとする

# 解凍済みファイルを直接指定する場合はここに記入（NULLなら自動検出）
MANUAL_MTX      <- NULL   # 例: "cb_adult_mouse.mtx.gz"
MANUAL_BARCODES <- NULL   # 例: "cb_adult_mouse_barcodes"
MANUAL_GENES    <- NULL   # 例: "cb_adult_mouse_genes"
MANUAL_META     <- NULL   # 例: "cb_adult_mouse_metadata.txt"（なければNULL）

dir.create(OUT_DIR, showWarnings = FALSE)

# ---- ユーティリティ ----
find_gene <- function(gene, all_genes) {
  if (gene %in% all_genes) return(gene)
  m <- all_genes[tolower(all_genes) == tolower(gene)]
  if (length(m) > 0) return(m[1])
  NA_character_
}

# ================================================================
# Step 1: ファイルの探索（解凍済み優先、なければtar.gzを解凍）
# ================================================================
message("\n=== Step 1: ファイルの探索 ===")

# 検索対象ディレクトリ（カレント + extracted/）
search_dirs <- c(".", file.path(OUT_DIR, "extracted"))

# カレントディレクトリ直下に解凍済みファイルがあるか確認
pre_extracted <- list.files(".", recursive = FALSE, full.names = TRUE)
pre_extracted <- c(pre_extracted,
                   list.files(file.path(OUT_DIR, "extracted"),
                              recursive = TRUE, full.names = TRUE))

has_mtx <- any(grepl("\\.mtx(\\.gz)?$", pre_extracted, ignore.case = TRUE))

if (has_mtx) {
  message("解凍済みファイルを検出 → 解凍をスキップ")
  all_files <- pre_extracted
} else if (file.exists(TAR_FILE)) {
  extract_dir <- file.path(OUT_DIR, "extracted")
  dir.create(extract_dir, showWarnings = FALSE)
  message("解凍中: ", TAR_FILE, " → ", extract_dir)
  untar(TAR_FILE, exdir = extract_dir)
  all_files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
} else {
  stop("データファイルが見つかりません。\n",
       "・tar.gzが未解凍の場合: ", TAR_FILE, " をこのスクリプトと同じフォルダに置いてください。\n",
       "・解凍済みの場合: mtx, barcodes, genes ファイルが同じフォルダにあるか確認してください。")
}

message("検出ファイル:")
for (f in head(all_files, 20)) message("  ", basename(f))
if (length(all_files) > 20) message("  ...他 ", length(all_files) - 20, " ファイル")

# ================================================================
# Step 2: データのロード
# ================================================================
message("\n=== Step 2: データのロード ===")

find_file <- function(paths, patterns) {
  for (pat in patterns)
    for (p in paths)
      if (grepl(pat, basename(p), ignore.case = TRUE)) return(p)
  NULL
}

seurat_obj <- NULL

# 手動指定があればそれを優先
if (!is.null(MANUAL_MTX)) {
  mtx_file  <- MANUAL_MTX
  bar_file  <- MANUAL_BARCODES
  feat_file <- MANUAL_GENES
} else {
  # (A) Seurat RDS
  rds_file <- find_file(all_files, c("\\.rds$"))
  if (!is.null(rds_file) && has_seurat) {
    message("Seuratオブジェクトを読み込み: ", basename(rds_file))
    seurat_obj <- readRDS(rds_file)
    message("読み込み完了: ", ncol(seurat_obj), " 細胞 × ", nrow(seurat_obj), " 遺伝子")
  }

  # (B) mtx + barcodes + genes（拡張子なしも対応）
  mtx_file  <- find_file(all_files, c("\\.mtx(\\.gz)?$"))
  bar_file  <- find_file(all_files, c("barcodes(\\.tsv|\\.txt)?(\\.gz)?$",
                                      "barcodes(\\.txt)?$"))
  feat_file <- find_file(all_files, c("genes(\\.tsv|\\.txt)?(\\.gz)?$",
                                      "features(\\.tsv|\\.txt)?(\\.gz)?$",
                                      "genes(\\.txt)?$", "features(\\.txt)?$"))
}

if (is.null(seurat_obj) && !is.null(mtx_file) && !is.null(bar_file) && !is.null(feat_file)) {
  message("10x形式を読み込み中（対象遺伝子のみ抽出・メモリ節約モード）...")
  message("  mtx      : ", basename(mtx_file))
  message("  barcodes : ", basename(bar_file))
  message("  genes    : ", basename(feat_file))

  # 遺伝子リストを読み込む（軽い）
  features   <- read.table(normalizePath(feat_file),
                           header = FALSE, sep = "\t",
                           stringsAsFactors = FALSE)
  gene_names <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  barcodes   <- read.table(normalizePath(bar_file),
                           header = FALSE, stringsAsFactors = FALSE)[[1]]
  message("  総遺伝子数: ", length(gene_names))
  message("  総細胞数  : ", length(barcodes))

  # 対象遺伝子の行インデックスを特定
  target_idx <- which(tolower(gene_names) %in% tolower(GENES_OF_INT))
  found_names <- gene_names[target_idx]
  message("  抽出遺伝子: ", paste(found_names, collapse = ", "),
          " (行インデックス: ", paste(target_idx, collapse = ", "), ")")

  if (length(target_idx) == 0)
    stop("対象遺伝子（", paste(GENES_OF_INT, collapse = ", "),
         "）がgenesファイルに見つかりません。")

  # mtxファイルをスキャンして対象遺伝子の行のみ取り込む
  message("  mtxファイルをスキャン中（数分かかります）...")

  # mtxは .gz の場合 gzcon で開く
  open_mtx <- function(f) {
    if (grepl("\\.gz$", f)) gzcon(file(normalizePath(f), "rb")) else file(normalizePath(f), "r")
  }
  con <- open_mtx(mtx_file)
  on.exit(close(con), add = TRUE)

  # ヘッダー行をスキップ（%で始まる行）
  repeat {
    line <- readLines(con, n = 1)
    if (!startsWith(line, "%")) break
  }
  # 次の行: nrow ncol nnz
  dims    <- as.integer(strsplit(trimws(line), "\\s+")[[1]])
  n_genes <- dims[1]; n_cells <- dims[2]; nnz <- dims[3]
  message(sprintf("  行列次元: %d 遺伝子 × %d 細胞, 非ゼロ要素: %d",
                  n_genes, n_cells, nnz))

  # 対象遺伝子の非ゼロ要素を収集
  rows_list <- vector("list", length(target_idx))
  names(rows_list) <- as.character(target_idx)
  for (i in seq_along(target_idx))
    rows_list[[i]] <- list(col = integer(0), val = numeric(0))

  chunk <- 1e6L
  read_so_far <- 0L
  repeat {
    lines <- readLines(con, n = chunk)
    if (length(lines) == 0) break
    read_so_far <- read_so_far + length(lines)
    if (read_so_far %% 5e6 == 0)
      message(sprintf("    %.0f%% スキャン済み...", 100 * read_so_far / nnz))

    # 数値に変換
    spl <- strsplit(lines, " ", fixed = TRUE)
    mat_chunk <- matrix(as.numeric(unlist(spl, use.names = FALSE)),
                        ncol = 3, byrow = TRUE)
    for (ti in seq_along(target_idx)) {
      keep <- mat_chunk[, 1] == target_idx[ti]
      if (any(keep)) {
        rows_list[[ti]]$col <- c(rows_list[[ti]]$col, as.integer(mat_chunk[keep, 2]))
        rows_list[[ti]]$val <- c(rows_list[[ti]]$val, mat_chunk[keep, 3])
      }
    }
  }
  message("  スキャン完了")

  # スパース行列を構築
  count_list <- lapply(seq_along(target_idx), function(ti) {
    Matrix::sparseMatrix(
      i = rep(1L, length(rows_list[[ti]]$col)),
      j = rows_list[[ti]]$col,
      x = rows_list[[ti]]$val,
      dims = c(1L, length(barcodes)),
      dimnames = list(found_names[ti], barcodes)
    )
  })
  counts <- do.call(rbind, count_list)
  message("  行列構築完了: ", nrow(counts), " 遺伝子 × ", ncol(counts), " 細胞")

  seurat_obj <- list(
    counts = counts,
    meta   = data.frame(row.names = barcodes,
                        cell      = barcodes)
  )
  message("読み込み完了")
}

# (C) h5
if (is.null(seurat_obj)) {
  h5_file <- find_file(all_files, c("\\.h5$"))
  if (!is.null(h5_file) && has_seurat) {
    message("H5形式を読み込み: ", basename(h5_file))
    counts <- Seurat::Read10X_h5(h5_file)
    seurat_obj <- Seurat::CreateSeuratObject(counts = counts,
                                             project = "GSE165371")
  }
}

if (is.null(seurat_obj))
  stop("認識できるデータファイルが見つかりません。\n",
       "検出ファイル:\n", paste(basename(all_files), collapse = "\n"))

# ================================================================
# Step 3: メタデータ（細胞タイプアノテーション）の確認
# ================================================================
message("\n=== Step 3: メタデータの確認 ===")

if (has_seurat && inherits(seurat_obj, "Seurat")) {
  meta <- seurat_obj@meta.data
} else {
  meta <- seurat_obj$meta
}

message("メタデータ列: ", paste(colnames(meta), collapse = ", "))

# 細胞タイプ列を探す
ct_candidates <- c("cell_type", "CellType", "celltype", "cluster", "ClusterID",
                   "seurat_clusters", "leiden", "louvain", "annotation",
                   "SubType", "subtype", "cell_type_label", "cluster_label",
                   "orig.ident", "Cluster")
ct_col <- ct_candidates[ct_candidates %in% colnames(meta)][1]

if (!is.na(ct_col)) {
  message("細胞タイプ列: ", ct_col)
  ct_table <- sort(table(meta[[ct_col]]), decreasing = TRUE)
  message("細胞タイプ一覧:")
  print(ct_table)
} else {
  message("⚠ 細胞タイプ列が見つかりません。利用可能な列:")
  print(head(meta, 3))
}

# メタデータファイルが別途ある場合は読み込む
meta_files <- all_files[grepl("meta|annot|cluster|barcode.*label|cell.*type",
                               basename(all_files), ignore.case = TRUE) &
                         grepl("\\.csv$|\\.tsv$|\\.txt$", all_files)]
if (length(meta_files) > 0 && is.na(ct_col)) {
  message("メタデータファイルを追加読み込み: ", basename(meta_files[1]))
  sep <- if (grepl("\\.tsv$|\\.txt$", meta_files[1])) "\t" else ","
  extra_meta <- read.table(meta_files[1], sep = sep, header = TRUE,
                           row.names = 1, check.names = FALSE)
  message("追加メタデータ列: ", paste(colnames(extra_meta), collapse = ", "))
  if (has_seurat && inherits(seurat_obj, "Seurat")) {
    common <- intersect(colnames(seurat_obj), rownames(extra_meta))
    seurat_obj <- seurat_obj[, common]
    seurat_obj <- Seurat::AddMetaData(seurat_obj,
                                      extra_meta[common, , drop = FALSE])
    meta <- seurat_obj@meta.data
    ct_col <- ct_candidates[ct_candidates %in% colnames(meta)][1]
  }
}

# ================================================================
# Step 4: 正規化
# ================================================================
message("\n=== Step 4: 正規化 ===")

if (has_seurat && inherits(seurat_obj, "Seurat")) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj,
                                      normalization.method = "LogNormalize",
                                      scale.factor = 1e4,
                                      verbose = FALSE)
  norm_mat <- Seurat::GetAssayData(seurat_obj, slot = "data")
  message("LogNormalize完了（log(CP10K + 1)）")
} else {
  # メモリ節約モード: 3遺伝子のみ読み込んでいるため
  # ライブラリサイズが不明 → 生UMIカウントをそのまま使用
  # （発現の有無＝UMI > 0 で判定するため正規化なしでも正確）
  norm_mat <- seurat_obj$counts
  message("生UMIカウントを使用（3遺伝子のみ読み込みのためライブラリサイズ正規化は省略）")
  message("※ 発現の有無（UMI > 0）および平均UMI数で比較します")
}

# ================================================================
# Step 5: 対象遺伝子の確認
# ================================================================
message("\n=== Step 5: 対象遺伝子の確認 ===")

all_genes <- rownames(norm_mat)
gene_map  <- setNames(sapply(GENES_OF_INT, find_gene, all_genes = all_genes),
                      GENES_OF_INT)

for (nm in names(gene_map)) {
  status <- if (!is.na(gene_map[nm])) paste0("✓ (", gene_map[nm], ")") else "✗ 未検出"
  message("  ", nm, " → ", status)
}

if (is.na(gene_map["Gpr176"]))
  stop("Gpr176が発現行列に見つかりません。")

# ================================================================
# Step 6: プルキンエ細胞の特定 と Gpr176陽性クラスターの同定
# ================================================================
message("\n=== Step 6: プルキンエ細胞 & Gpr176陽性クラスターの同定 ===")

gpr176_key <- gene_map["Gpr176"]
gpr176_exp <- as.numeric(norm_mat[gpr176_key, ])

if (has_seurat && inherits(seurat_obj, "Seurat")) {
  meta <- seurat_obj@meta.data
} else {
  meta <- seurat_obj$meta
}

meta$Gpr176_expr     <- gpr176_exp
meta$Gpr176_positive <- gpr176_exp > GPR176_THRESHOLD
meta$cell            <- rownames(meta)

n_pos <- sum(meta$Gpr176_positive)
n_tot <- nrow(meta)
message(sprintf("Gpr176陽性細胞: %d / %d (%.1f%%)", n_pos, n_tot, 100 * n_pos / n_tot))

# プルキンエ細胞フラグ（アノテーション列がある場合）
purkinje_keywords <- c("purkinje", "Purkinje", "PC", "PurkCell", "purk")
if (!is.na(ct_col)) {
  meta$is_purkinje <- grepl(paste(purkinje_keywords, collapse = "|"),
                             meta[[ct_col]], ignore.case = TRUE)
  n_purk <- sum(meta$is_purkinje)
  message(sprintf("プルキンエ細胞数: %d (%.1f%%)", n_purk, 100 * n_purk / n_tot))

  # クラスター別Gpr176陽性率
  cluster_summary <- meta %>%
    dplyr::group_by(cluster = .data[[ct_col]]) %>%
    dplyr::summarise(
      n_cells       = dplyr::n(),
      n_gpr176_pos  = sum(Gpr176_positive),
      pct_gpr176    = 100 * mean(Gpr176_positive),
      mean_gpr176   = mean(Gpr176_expr),
      is_purkinje   = any(grepl(paste(purkinje_keywords, collapse = "|"),
                                dplyr::cur_group()[[1]], ignore.case = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(pct_gpr176))

  write.csv(cluster_summary,
            file.path(OUT_DIR, "cluster_Gpr176_summary.csv"), row.names = FALSE)
  message("\nクラスター別Gpr176陽性率 (上位15):")
  print(as.data.frame(head(cluster_summary, 15)))

  pos_clusters <- cluster_summary$cluster[cluster_summary$pct_gpr176 > POS_PCT_THR]
  if (length(pos_clusters) == 0) {
    pos_clusters <- cluster_summary$cluster[1]
    message("⚠ 閾値超えクラスターなし → 最高陽性率クラスターを使用: ", pos_clusters)
  } else {
    message("Gpr176陽性クラスター (>", POS_PCT_THR, "%): ",
            paste(pos_clusters, collapse = ", "))
  }
  meta$Gpr176_cluster <- ifelse(meta[[ct_col]] %in% pos_clusters,
                                 "Gpr176_positive", "Gpr176_negative")
} else {
  message("細胞タイプアノテーションなし → 細胞単位で陽性/陰性を使用")
  cluster_summary <- NULL
  pos_clusters    <- character(0)
  meta$is_purkinje    <- FALSE
  meta$Gpr176_cluster <- ifelse(meta$Gpr176_positive,
                                 "Gpr176_positive", "Gpr176_negative")
}

# ================================================================
# Step 7: GnazとRGS16の発現量解析
# ================================================================
message("\n=== Step 7: Gnaz・RGS16の発現量解析 ===")

target_genes <- gene_map[names(gene_map) %in% c("Gnaz", "Rgs16")]
target_genes <- target_genes[!is.na(target_genes)]

result_rows <- list()
for (gname in names(target_genes)) {
  actual   <- target_genes[gname]
  expr     <- as.numeric(norm_mat[actual, ])
  pos_vals <- expr[meta$Gpr176_cluster == "Gpr176_positive"]
  neg_vals <- expr[meta$Gpr176_cluster == "Gpr176_negative"]

  wt      <- wilcox.test(pos_vals, neg_vals, alternative = "greater")
  pct_pos <- 100 * mean(pos_vals > 0)
  pct_neg <- 100 * mean(neg_vals > 0)
  log2fc  <- log2((mean(pos_vals) + 1e-6) / (mean(neg_vals) + 1e-6))
  sig     <- dplyr::case_when(wt$p.value < 0.001 ~ "***",
                               wt$p.value < 0.01  ~ "**",
                               wt$p.value < 0.05  ~ "*",
                               TRUE               ~ "ns")

  result_rows[[gname]] <- data.frame(
    gene = gname, actual_name = actual,
    mean_Gpr176pos = round(mean(pos_vals), 4),
    mean_Gpr176neg = round(mean(neg_vals), 4),
    pct_expressed_pos = round(pct_pos, 2),
    pct_expressed_neg = round(pct_neg, 2),
    log2FC = round(log2fc, 4),
    wilcox_pval = signif(wt$p.value, 4),
    significance = sig,
    stringsAsFactors = FALSE
  )

  message(sprintf("\n[%s (%s)]", gname, actual))
  message(sprintf("  Gpr176+: 平均=%.4f, 発現率=%.1f%%", mean(pos_vals), pct_pos))
  message(sprintf("  Gpr176-: 平均=%.4f, 発現率=%.1f%%", mean(neg_vals), pct_neg))
  message(sprintf("  log2FC=%.3f, Wilcoxon p=%s %s", log2fc, signif(wt$p.value, 4), sig))
}

result_df <- do.call(rbind, result_rows)
write.csv(result_df,
          file.path(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"), row.names = FALSE)

# ================================================================
# Step 8: 可視化
# ================================================================
message("\n=== Step 8: 可視化 ===")

colors_group <- c("Gpr176_positive" = "#E64B35", "Gpr176_negative" = "#4DBBD5")

# --- 8-A: バイオリンプロット（Gpr176陽性 vs 陰性） ---
plot_genes <- c("Gpr176", names(target_genes))
actual_map <- c("Gpr176" = gpr176_key, target_genes)

vln_list <- lapply(plot_genes, function(gname) {
  actual <- actual_map[gname]
  if (is.na(actual)) return(NULL)
  expr_vec <- as.numeric(norm_mat[actual, ])
  pdata    <- data.frame(expr = expr_vec,
                         group = meta$Gpr176_cluster,
                         stringsAsFactors = FALSE)
  ggplot2::ggplot(pdata, ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.75) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.3,
                          fill = "white", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = colors_group) +
    ggplot2::scale_x_discrete(labels = c("Gpr176_positive" = "Gpr176+",
                                         "Gpr176_negative" = "Gpr176-")) +
    ggplot2::labs(title = sprintf("%s (%s)", gname, actual),
                  x = NULL, y = "log(CP10K + 1)") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", size = 11))
})
vln_list <- Filter(Negate(is.null), vln_list)

combined_vln <- patchwork::wrap_plots(vln_list, ncol = length(vln_list))
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.pdf"),
                combined_vln, width = 4 * length(vln_list), height = 5)
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.png"),
                combined_vln, width = 4 * length(vln_list), height = 5, dpi = 150)
message("バイオリンプロット保存完了")

# --- 8-B: 発現率棒グラフ ---
bar_df <- result_df %>%
  dplyr::select(gene, pct_expressed_pos, pct_expressed_neg) %>%
  tidyr::pivot_longer(cols = c(pct_expressed_pos, pct_expressed_neg),
                      names_to = "group", values_to = "pct") %>%
  dplyr::mutate(group = dplyr::recode(group,
    "pct_expressed_pos" = "Gpr176+ cluster",
    "pct_expressed_neg" = "Gpr176- cluster"))

bar_p <- ggplot2::ggplot(bar_df, ggplot2::aes(x = gene, y = pct, fill = group)) +
  ggplot2::geom_col(position = "dodge", width = 0.6,
                    color = "black", linewidth = 0.3) +
  ggplot2::geom_text(data = result_df,
    ggplot2::aes(x = gene,
                 y = pmax(pct_expressed_pos, pct_expressed_neg) + 3,
                 label = significance),
    inherit.aes = FALSE, size = 5, fontface = "bold") +
  ggplot2::scale_fill_manual(
    values = c("Gpr176+ cluster" = "#E64B35",
               "Gpr176- cluster" = "#4DBBD5")) +
  ggplot2::labs(
    title    = "Gpr176陽性クラスターにおけるGnaz・RGS16の発現細胞率",
    subtitle = "Kozareva et al. 2021 (GSE165371)",
    x = "遺伝子", y = "発現細胞率 (%)", fill = NULL) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = "top")

ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.pdf"),
                bar_p, width = 6, height = 5)
ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.png"),
                bar_p, width = 6, height = 5, dpi = 150)
message("棒グラフ保存完了")

# --- 8-C: クラスター別ドットプロット ---
if (!is.null(cluster_summary) && !is.na(ct_col)) {
  show_clusters <- cluster_summary$cluster[seq_len(min(20, nrow(cluster_summary)))]
  dot_data <- lapply(show_clusters, function(cl) {
    cl_mask <- meta[[ct_col]] == cl
    lapply(plot_genes, function(gname) {
      actual <- actual_map[gname]
      if (is.na(actual)) return(NULL)
      vals <- as.numeric(norm_mat[actual, cl_mask])
      data.frame(cluster = as.character(cl), gene = gname,
                 mean_expr = mean(vals), pct = 100 * mean(vals > 0))
    })
  })
  dot_df <- do.call(rbind, do.call(c, dot_data))
  dot_df  <- dot_df[!is.null(dot_df), ]

  # クラスターをGpr176陽性率順に並べる
  cl_order <- cluster_summary$cluster[seq_len(min(20, nrow(cluster_summary)))]
  dot_df$cluster <- factor(dot_df$cluster, levels = rev(as.character(cl_order)))

  dot_p <- ggplot2::ggplot(dot_df,
             ggplot2::aes(x = gene, y = cluster, size = pct, color = mean_expr)) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis(option = "plasma",
                                 name = "平均発現量\nlog(CP10K+1)") +
    ggplot2::scale_size_continuous(range = c(0.5, 8),
                                   name = "発現細胞率 (%)") +
    ggplot2::labs(
      title    = "Gpr176 / Gnaz / Rgs16 発現（クラスター別）",
      subtitle = "Kozareva et al. 2021 (GSE165371) — Gpr176陽性率降順",
      x = "遺伝子", y = "クラスター") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))

  h <- max(4, length(unique(dot_df$cluster)) * 0.45 + 2)
  ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.pdf"),
                  dot_p, width = 7, height = h)
  ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.png"),
                  dot_p, width = 7, height = h, dpi = 150)
  message("ドットプロット保存完了")
}

# --- 8-D: プルキンエ細胞に絞ったバイオリンプロット ---
if (any(meta$is_purkinje)) {
  purk_meta <- meta[meta$is_purkinje, ]
  purk_vln_list <- lapply(plot_genes, function(gname) {
    actual <- actual_map[gname]
    if (is.na(actual)) return(NULL)
    expr_vec <- as.numeric(norm_mat[actual, purk_meta$cell])
    pdata    <- data.frame(expr = expr_vec,
                           cluster = purk_meta[[ct_col]],
                           stringsAsFactors = FALSE)
    ggplot2::ggplot(pdata, ggplot2::aes(x = cluster, y = expr, fill = cluster)) +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.75) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.3,
                            fill = "white", alpha = 0.8) +
      ggplot2::labs(title = gname, x = NULL, y = "log(CP10K + 1)") +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(legend.position = "none",
                     axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                     plot.title = ggplot2::element_text(face = "bold.italic"))
  })
  purk_vln_list <- Filter(Negate(is.null), purk_vln_list)

  if (length(purk_vln_list) > 0) {
    purk_combined <- patchwork::wrap_plots(purk_vln_list,
                                           ncol = length(purk_vln_list))
    ggplot2::ggsave(file.path(OUT_DIR, "violin_purkinje_Gpr176_Gnaz_RGS16.pdf"),
                    purk_combined, width = 4 * length(purk_vln_list), height = 5)
    ggplot2::ggsave(file.path(OUT_DIR, "violin_purkinje_Gpr176_Gnaz_RGS16.png"),
                    purk_combined, width = 4 * length(purk_vln_list), height = 5,
                    dpi = 150)
    message("プルキンエ細胞バイオリンプロット保存完了")
  }
}

# ================================================================
# 最終サマリー
# ================================================================
message("\n", strrep("=", 54))
message("           解析完了サマリー")
message(strrep("=", 54))
message("データセット : GSE165371 (Kozareva et al. 2021)")
message(sprintf("総細胞数     : %d", n_tot))
message(sprintf("Gpr176陽性   : %d 細胞 (%.1f%%)", n_pos, 100 * n_pos / n_tot))
if (!is.na(ct_col))
  message(sprintf("陽性クラスター: %s",
                  paste(head(pos_clusters, 5), collapse = ", ")))
message("")
message("Gpr176陽性クラスターでの発現:")
for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  message(sprintf("  %-6s: 発現率 %5.1f%% (vs %5.1f%%)  log2FC=%+.2f  p=%s %s",
                  r$gene, r$pct_expressed_pos, r$pct_expressed_neg,
                  r$log2FC, r$wilcox_pval, r$significance))
}
message("")
message("出力フォルダ: ", OUT_DIR, "/")
message("  cluster_Gpr176_summary.csv")
message("  Gnaz_RGS16_in_Gpr176clusters.csv")
message("  violin_Gnaz_RGS16_Gpr176clusters.pdf/png")
message("  dotplot_Gpr176_Gnaz_RGS16.pdf/png")
message("  barplot_expression_pct.pdf/png")
message("  violin_purkinje_Gpr176_Gnaz_RGS16.pdf/png  (プルキンエ細胞のみ)")
message(strrep("=", 54))
