# ============================================================
# setup_and_run.R
# 必要パッケージのインストール確認 → 解析実行
# ローカル環境で最初にこのファイルを実行してください
# ============================================================

# ---- BiocManagerのインストール ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# ---- 必要パッケージのインストール ----
cran_pkgs  <- c("ggplot2", "dplyr", "patchwork", "viridis", "tidyr", "Matrix", "R.utils")
bioc_pkgs  <- c("GEOquery", "Seurat")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

cat("パッケージインストール完了\n")

# ---- 解析本体を実行 ----
source("Gpr176_cluster_analysis.R")
