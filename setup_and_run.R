# ============================================================
# setup_and_run.R
# 必要パッケージのインストール確認 → 解析実行
#
# 事前準備:
#   GSE165371_cb_adult_mouse.tar.gz を
#   このファイルと同じフォルダに置いてから実行してください。
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran_pkgs <- c("ggplot2", "dplyr", "patchwork", "viridis",
               "tidyr", "Matrix", "stringr")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

# Seurat（あれば高速・多機能になるが必須ではない）
if (!requireNamespace("Seurat", quietly = TRUE)) {
  message("Seuratをインストール中（時間がかかる場合があります）...")
  install.packages("Seurat", repos = "https://cloud.r-project.org")
}

cat("パッケージ準備完了\n")
source("Gpr176_cluster_analysis.R")
