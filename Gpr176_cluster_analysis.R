# ============================================================
# Gpr176_cluster_analysis.R
#
# 論文: Kozareva et al. 2021, Nature
# "A transcriptomic atlas of mouse cerebellar cortex
#  comprehensively defines cell types"
# GEO: GSE165805
#
# 目的: Gpr176陽性クラスターにGnaz（Gz）とRGS16が発現するかを検証
# ============================================================

# ---- パッケージ ----
required_cran <- c("ggplot2", "dplyr", "patchwork", "viridis",
                   "tidyr", "Matrix", "readxl", "stringr")
required_bioc <- c("GEOquery")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  library(pkg, character.only = TRUE)
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  library(pkg, character.only = TRUE)
}

# ---- 設定 ----
GEO_ID           <- "GSE165805"
GENES_OF_INT     <- c("Gpr176", "Gnaz", "Rgs16")
OUT_DIR          <- "Gpr176_cluster_results"
GPR176_THRESHOLD <- 0      # log-norm > 0 で陽性
POS_PCT_THR      <- 10     # クラスター内Gpr176陽性率(%)がこれ以上を陽性クラスターとする

dir.create(OUT_DIR, showWarnings = FALSE)

# ---- ユーティリティ ----
find_gene <- function(gene, all_genes) {
  if (gene %in% all_genes) return(gene)
  m <- all_genes[tolower(all_genes) == tolower(gene)]
  if (length(m) > 0) return(m[1])
  NA_character_
}

# ================================================================
# Step 1: GEOから補足ファイルをダウンロード
# ================================================================
message("\n=== Step 1: GEOデータ取得 (", GEO_ID, ") ===")

supp_dir <- file.path(OUT_DIR, "supp_files")
dir.create(supp_dir, showWarnings = FALSE)

supp_files <- GEOquery::getGEOSuppFiles(GEO_ID, baseDir = supp_dir, fetch_files = TRUE)
supp_paths <- rownames(supp_files)
message("ダウンロードされたファイル:")
for (p in supp_paths) message("  ", basename(p))

# ================================================================
# Step 2: Excelファイルの読み込みと結合
# ================================================================
message("\n=== Step 2: Excelファイルの読み込み ===")

xlsx_files <- supp_paths[grepl("\\.xlsx$", supp_paths, ignore.case = TRUE)]

if (length(xlsx_files) == 0) {
  stop("Excelファイルが見つかりません。\nダウンロードされたファイル:\n",
       paste(supp_paths, collapse = "\n"))
}

message(length(xlsx_files), "個のExcelファイルを検出:")
for (f in xlsx_files) message("  ", basename(f))

# 各Excelファイルを読み込んで結合
read_geo_excel <- function(path) {
  sheets <- readxl::excel_sheets(path)
  message("  シート: ", paste(sheets, collapse = ", "))

  df <- as.data.frame(
    readxl::read_excel(path, sheet = sheets[1], col_names = TRUE),
    stringsAsFactors = FALSE
  )
  message(sprintf("  読込直後: %d 行 × %d 列", nrow(df), ncol(df)))
  message("  列1-3の型: ", paste(sapply(df[, seq_len(min(3, ncol(df)))], class), collapse = ", "))

  # 文字列列を特定
  is_char <- vapply(df, function(x) is.character(x) || is.factor(x), logical(1L))
  char_cols <- which(is_char)
  num_cols  <- which(!is_char)

  # 文字列列が1列だけ → 遺伝子名列と判断
  if (length(char_cols) == 1) {
    gene_names <- make.unique(as.character(df[[char_cols]]))
    df <- df[, num_cols, drop = FALSE]
  # 文字列列が複数 → 1列目だけ遺伝子名として使う
  } else if (length(char_cols) >= 1) {
    gene_names <- make.unique(as.character(df[[char_cols[1]]]))
    df <- df[, num_cols, drop = FALSE]
  } else {
    # 文字列列なし: 行名をそのまま使う
    gene_names <- rownames(df)
  }

  # 数値行列に変換（dimnames を明示的に設定）
  mat <- matrix(
    as.double(unlist(df, use.names = FALSE)),
    nrow     = nrow(df),
    ncol     = ncol(df),
    dimnames = list(gene_names, colnames(df))
  )
  message(sprintf("  行列サイズ: %d 遺伝子 × %d 細胞", nrow(mat), ncol(mat)))
  mat
}

mat_list <- lapply(xlsx_files, function(f) {
  message("読み込み中: ", basename(f))
  read_geo_excel(f)
})
names(mat_list) <- basename(xlsx_files)

# 行（遺伝子）名を確認
gene_sets <- lapply(mat_list, rownames)
common_genes <- Reduce(intersect, gene_sets)
message("\n共通遺伝子数: ", length(common_genes),
        " (全ファイル共通)")

# ファイルを列方向に結合（細胞を横に並べる）
if (length(mat_list) == 1) {
  count_mat <- mat_list[[1]]
} else {
  mats_aligned <- lapply(mat_list, function(m) m[common_genes, , drop = FALSE])
  count_mat <- do.call(cbind, mats_aligned)
  # 細胞名の重複を避けるためファイル名プレフィックスを付与
  prefixes <- sub("_GEO_processed_data_", "_",
                  sub("\\.xlsx$", "", basename(xlsx_files)))
  new_colnames <- unlist(mapply(function(m, pfx) paste0(pfx, "_", colnames(m)),
                                mats_aligned, prefixes, SIMPLIFY = FALSE))
  colnames(count_mat) <- new_colnames
}

message("発現行列サイズ: ", nrow(count_mat), " 遺伝子 × ", ncol(count_mat), " 細胞")

# ---- データ種別の推定（生カウント vs 正規化済み） ----
sample_vals <- count_mat[count_mat > 0]
is_raw_count <- all(sample_vals == floor(sample_vals)) && max(sample_vals) > 100
message("データ種別: ", if (is_raw_count) "生カウント（正規化を実施）" else "正規化済み（正規化をスキップ）")

# ================================================================
# Step 3: 遺伝子の確認
# ================================================================
message("\n=== Step 3: 対象遺伝子の確認 ===")

all_genes <- rownames(count_mat)
gene_map <- setNames(
  sapply(GENES_OF_INT, find_gene, all_genes = all_genes),
  GENES_OF_INT
)

for (nm in names(gene_map)) {
  status <- if (!is.na(gene_map[nm])) paste0("✓ (", gene_map[nm], ")") else "✗ 未検出"
  message("  ", nm, " -> ", status)
}

if (is.na(gene_map["Gpr176"])) {
  stop("Gpr176が発現行列に見つかりません。解析を中断します。")
}

# ================================================================
# Step 4: 正規化
# ================================================================
message("\n=== Step 4: 正規化 ===")

if (is_raw_count) {
  # log CP10K
  col_sums   <- colSums(count_mat)
  norm_mat   <- sweep(count_mat, 2, col_sums / 1e4, "/")
  norm_mat   <- log1p(norm_mat)
  message("LogNormalize完了（log(CP10K + 1)）")
} else {
  norm_mat <- count_mat
  message("正規化済みデータをそのまま使用")
}

# ================================================================
# Step 5: Gpr176陽性クラスターの同定
# ================================================================
message("\n=== Step 5: Gpr176陽性クラスターの同定 ===")

gpr176_key <- gene_map["Gpr176"]
gpr176_exp <- as.numeric(norm_mat[gpr176_key, ])
names(gpr176_exp) <- colnames(norm_mat)

gpr176_positive <- gpr176_exp > GPR176_THRESHOLD
n_pos <- sum(gpr176_positive)
n_tot <- ncol(norm_mat)
message(sprintf("Gpr176陽性細胞: %d / %d (%.1f%%)", n_pos, n_tot, 100 * n_pos / n_tot))

# ---- クラスター情報の推定 ----
# GSE165805のファイル名には細胞タイプが含まれる (例: _Cd, _InP)
# 細胞名のプレフィックスからクラスターを推定
cell_names <- colnames(norm_mat)

# ファイル由来のクラスターラベルを付与
cluster_from_file <- sub("_[^_]+$", "",   # 最後の_以降（細胞番号）を除く
                         sub("^GSE165805_GEO_processed_data_", "", cell_names))

# ファイルが1つの場合は列名（細胞名）からクラスターを推定
if (length(xlsx_files) > 1) {
  cluster_labels <- cluster_from_file
  message("クラスター（ファイル由来）: ", paste(unique(cluster_labels), collapse = ", "))
} else {
  # 単一ファイル: 細胞名にクラスター情報が含まれていることが多い
  # 例: "AAACCTGAGAAACGCC-1_Granule" -> "Granule"
  cluster_try <- sub(".*_", "", cell_names)
  if (length(unique(cluster_try)) > 1 && length(unique(cluster_try)) < n_tot * 0.5) {
    cluster_labels <- cluster_try
    message("クラスター（細胞名末尾）: ", paste(head(unique(cluster_labels), 10), collapse = ", "))
  } else {
    # クラスター情報なし: Gpr176の発現量でビン分け
    cluster_labels <- rep("all_cells", n_tot)
    message("クラスター情報なし: 細胞単位でGpr176陽性/陰性を使用")
  }
}

# ---- クラスター別Gpr176陽性率 ----
meta_df <- data.frame(
  cell            = cell_names,
  cluster         = cluster_labels,
  Gpr176_expr     = gpr176_exp,
  Gpr176_positive = gpr176_positive,
  stringsAsFactors = FALSE
)

cluster_summary <- meta_df %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_cells      = dplyr::n(),
    n_gpr176_pos = sum(Gpr176_positive),
    pct_gpr176   = 100 * mean(Gpr176_positive),
    mean_gpr176  = mean(Gpr176_expr),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(pct_gpr176))

write.csv(cluster_summary,
          file.path(OUT_DIR, "cluster_Gpr176_summary.csv"),
          row.names = FALSE)
message("\nクラスター別Gpr176陽性率:")
print(as.data.frame(cluster_summary))

# Gpr176陽性クラスターの定義
pos_clusters <- cluster_summary$cluster[cluster_summary$pct_gpr176 > POS_PCT_THR]

if (length(pos_clusters) == 0) {
  # 閾値を下げて最高陽性率クラスターを採用
  pos_clusters <- cluster_summary$cluster[1]
  message(sprintf("⚠ 陽性率%d%%超クラスターなし → 最高陽性率クラスターを使用: %s",
                  POS_PCT_THR, pos_clusters))
} else {
  message("Gpr176陽性クラスター (>", POS_PCT_THR, "%): ",
          paste(pos_clusters, collapse = ", "))
}

meta_df$Gpr176_cluster <- ifelse(meta_df$cluster %in% pos_clusters,
                                  "Gpr176_positive", "Gpr176_negative")

# ================================================================
# Step 6: GnazとRGS16の発現量解析
# ================================================================
message("\n=== Step 6: Gnaz・RGS16の発現量解析 ===")

target_genes <- gene_map[names(gene_map) %in% c("Gnaz", "Rgs16")]
target_genes <- target_genes[!is.na(target_genes)]

result_rows <- list()
for (gname in names(target_genes)) {
  actual  <- target_genes[gname]
  expr    <- as.numeric(norm_mat[actual, ])

  pos_vals <- expr[meta_df$Gpr176_cluster == "Gpr176_positive"]
  neg_vals <- expr[meta_df$Gpr176_cluster == "Gpr176_negative"]

  wt       <- wilcox.test(pos_vals, neg_vals, alternative = "greater")
  pct_pos  <- 100 * mean(pos_vals > 0)
  pct_neg  <- 100 * mean(neg_vals > 0)
  log2fc   <- log2((mean(pos_vals) + 1e-6) / (mean(neg_vals) + 1e-6))
  sig      <- dplyr::case_when(wt$p.value < 0.001 ~ "***",
                               wt$p.value < 0.01  ~ "**",
                               wt$p.value < 0.05  ~ "*",
                               TRUE               ~ "ns")

  result_rows[[gname]] <- data.frame(
    gene              = gname,
    actual_name       = actual,
    mean_Gpr176pos    = round(mean(pos_vals), 4),
    mean_Gpr176neg    = round(mean(neg_vals), 4),
    pct_expressed_pos = round(pct_pos, 2),
    pct_expressed_neg = round(pct_neg, 2),
    log2FC            = round(log2fc, 4),
    wilcox_pval       = signif(wt$p.value, 4),
    significance      = sig,
    stringsAsFactors  = FALSE
  )

  message(sprintf("\n[%s (%s)]", gname, actual))
  message(sprintf("  Gpr176+: 平均=%.4f, 発現率=%.1f%%", mean(pos_vals), pct_pos))
  message(sprintf("  Gpr176-: 平均=%.4f, 発現率=%.1f%%", mean(neg_vals), pct_neg))
  message(sprintf("  log2FC=%.3f, Wilcoxon p=%s %s", log2fc, signif(wt$p.value, 4), sig))
}

result_df <- do.call(rbind, result_rows)
write.csv(result_df,
          file.path(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"),
          row.names = FALSE)
message("\n結果CSV: ", file.path(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"))

# ================================================================
# Step 7: 可視化
# ================================================================
message("\n=== Step 7: 可視化 ===")

all_plot_genes <- c("Gpr176", names(target_genes))
colors_group <- c("Gpr176_positive" = "#E64B35", "Gpr176_negative" = "#4DBBD5")

# --- 7-A: バイオリンプロット ---
plot_list <- list()
for (gname in all_plot_genes) {
  actual <- if (gname == "Gpr176") gpr176_key else target_genes[gname]
  if (is.na(actual)) next

  expr_vec <- as.numeric(norm_mat[actual, ])
  pdata <- data.frame(
    expr  = expr_vec,
    group = meta_df$Gpr176_cluster,
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(pdata, ggplot2::aes(x = group, y = expr, fill = group)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.75) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.3,
                          fill = "white", alpha = 0.8) +
    ggplot2::scale_fill_manual(values = colors_group) +
    ggplot2::scale_x_discrete(labels = c("Gpr176_positive" = "Gpr176+",
                                         "Gpr176_negative" = "Gpr176-")) +
    ggplot2::labs(title = sprintf("%s\n(%s)", gname, actual),
                  x = NULL, y = "log(CP10K + 1)") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", size = 11))

  plot_list[[gname]] <- p
}

combined_violin <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.pdf"),
                combined_violin, width = 4 * length(plot_list), height = 5)
ggplot2::ggsave(file.path(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.png"),
                combined_violin, width = 4 * length(plot_list), height = 5, dpi = 150)
message("バイオリンプロット保存完了")

# --- 7-B: 発現率棒グラフ（p値注釈付き） ---
bar_df <- result_df %>%
  dplyr::select(gene, pct_expressed_pos, pct_expressed_neg) %>%
  tidyr::pivot_longer(cols = c(pct_expressed_pos, pct_expressed_neg),
                      names_to = "group", values_to = "pct") %>%
  dplyr::mutate(group = dplyr::recode(group,
                                      "pct_expressed_pos" = "Gpr176+ cluster",
                                      "pct_expressed_neg" = "Gpr176- cluster"))

bar_p <- ggplot2::ggplot(bar_df, ggplot2::aes(x = gene, y = pct, fill = group)) +
  ggplot2::geom_col(position = "dodge", width = 0.6, color = "black", linewidth = 0.3) +
  ggplot2::geom_text(data = result_df,
                     ggplot2::aes(x = gene, y = pmax(pct_expressed_pos, pct_expressed_neg) + 3,
                                  label = significance),
                     inherit.aes = FALSE, size = 5, fontface = "bold") +
  ggplot2::scale_fill_manual(values = c("Gpr176+ cluster" = "#E64B35",
                                        "Gpr176- cluster" = "#4DBBD5")) +
  ggplot2::labs(title = "Gpr176陽性クラスターにおけるGnaz・RGS16の発現細胞率",
                subtitle = "Kozareva et al. 2021 (GSE165805)",
                x = "遺伝子", y = "発現細胞率 (%)", fill = NULL) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(legend.position = "top")

ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.pdf"),
                bar_p, width = 6, height = 5)
ggplot2::ggsave(file.path(OUT_DIR, "barplot_expression_pct.png"),
                bar_p, width = 6, height = 5, dpi = 150)
message("棒グラフ保存完了")

# --- 7-C: クラスター別ドットプロット ---
dot_genes_actual <- c(gpr176_key, unname(target_genes))
dot_genes_actual <- dot_genes_actual[!is.na(dot_genes_actual)]
dot_gene_labels  <- c("Gpr176", names(target_genes))[
  !is.na(c(gpr176_key, unname(target_genes)))]

dot_data <- lapply(unique(meta_df$cluster), function(cl) {
  cl_cells <- meta_df$cell[meta_df$cluster == cl]
  lapply(seq_along(dot_genes_actual), function(i) {
    vals <- as.numeric(norm_mat[dot_genes_actual[i], cl_cells])
    data.frame(cluster = cl,
               gene    = dot_gene_labels[i],
               mean    = mean(vals),
               pct     = 100 * mean(vals > 0))
  })
})
dot_df <- do.call(rbind, do.call(c, dot_data))

dot_p <- ggplot2::ggplot(dot_df,
         ggplot2::aes(x = gene, y = cluster,
                      size = pct, color = mean)) +
  ggplot2::geom_point() +
  viridis::scale_color_viridis(option = "plasma", name = "平均発現量\nlog(CP10K+1)") +
  ggplot2::scale_size_continuous(range = c(1, 8), name = "発現細胞率 (%)") +
  ggplot2::labs(title = "Gpr176 / Gnaz / Rgs16 発現（クラスター別）",
                subtitle = "Kozareva et al. 2021 (GSE165805)",
                x = "遺伝子", y = "クラスター") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic"))

ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.pdf"),
                dot_p, width = max(5, length(dot_genes_actual) * 2),
                height = max(4, length(unique(meta_df$cluster)) * 0.5 + 2))
ggplot2::ggsave(file.path(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.png"),
                dot_p, width = max(5, length(dot_genes_actual) * 2),
                height = max(4, length(unique(meta_df$cluster)) * 0.5 + 2),
                dpi = 150)
message("ドットプロット保存完了")

# ================================================================
# 最終サマリー
# ================================================================
message("\n", strrep("=", 52))
message("           解析完了サマリー")
message(strrep("=", 52))
message(sprintf("データセット : %s (Kozareva et al. 2021)", GEO_ID))
message(sprintf("総細胞数     : %d", n_tot))
message(sprintf("Gpr176陽性   : %d 細胞 (%.1f%%)", n_pos, 100 * n_pos / n_tot))
message(sprintf("陽性クラスター: %s", paste(pos_clusters, collapse = ", ")))
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
message(strrep("=", 52))
