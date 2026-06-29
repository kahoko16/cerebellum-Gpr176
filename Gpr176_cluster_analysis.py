"""
Gpr176_cluster_analysis.py

論文: Kozareva et al. 2021, Nature
  "A transcriptomic atlas of mouse cerebellar cortex comprehensively defines cell types"
GEO: GSE165805

目的: Gpr176陽性クラスターにGnaz（Gz）とRGS16が発現するかを検証

解析の流れ:
  1. GEOから補足ファイル（カウントマトリクス + 細胞メタデータ）を取得
  2. Gpr176の発現に基づいてクラスターを陽性/陰性に分類
  3. 各群内でのGnaz・RGS16の発現量を統計検定および可視化
"""

import os, gzip, shutil, warnings, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.io import mmread
import GEOparse
import scanpy as sc

warnings.filterwarnings("ignore")

# ---- 設定 ----
GEO_ID      = "GSE165805"
GENES       = ["Gpr176", "Gnaz", "Rgs16"]
OUT_DIR     = "Gpr176_cluster_results"
GPR176_THR  = 0          # log-norm > 0 で陽性
POS_PCT_THR = 10.0       # クラスター内Gpr176陽性率(%)がこれ以上を「陽性クラスター」とする

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.verbosity = 1

# ---- ユーティリティ ----
def decompress(path):
    """gzipファイルを解凍して解凍後パスを返す"""
    if path and path.endswith(".gz"):
        out = path[:-3]
        if not os.path.exists(out):
            with gzip.open(path, "rb") as f_in, open(out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return out
    return path

def find_file(paths, patterns):
    """パスリストからpatternに一致する最初のファイルを返す"""
    for pat in patterns:
        for p in paths:
            if re.search(pat, os.path.basename(p), re.IGNORECASE):
                return p
    return None

def find_gene(gene, all_genes):
    """大文字小文字を無視して遺伝子名を検索"""
    # 完全一致
    if gene in all_genes:
        return gene
    # 大文字小文字無視
    lower_map = {g.lower(): g for g in all_genes}
    return lower_map.get(gene.lower())


# ================================================================
# Step 1: GEOから補足ファイルをダウンロード
# ================================================================
print(f"\n=== Step 1: GEOデータ取得 ({GEO_ID}) ===")
supp_dir = os.path.join(OUT_DIR, "supp_files")
os.makedirs(supp_dir, exist_ok=True)

gse = GEOparse.get_GEO(geo=GEO_ID, destdir=supp_dir, silent=True)
print("GSEタイトル:", gse.metadata.get("title", ["不明"])[0])

# 補足ファイルのパスを収集
supp_paths = []
for gsm_name, gsm in gse.gsms.items():
    for url in gsm.metadata.get("supplementary_file", []):
        fname = os.path.basename(url.rstrip("/"))
        local = os.path.join(supp_dir, fname)
        if not os.path.exists(local):
            print(f"  ダウンロード中: {fname}")
            try:
                import urllib.request
                urllib.request.urlretrieve(url, local)
            except Exception as e:
                print(f"  ⚠ ダウンロード失敗: {e}")
        if os.path.exists(local):
            supp_paths.append(local)

# シリーズレベルの補足ファイルも確認
for url in gse.metadata.get("supplementary_file", []):
    fname = os.path.basename(url.rstrip("/"))
    local = os.path.join(supp_dir, fname)
    if not os.path.exists(local):
        print(f"  シリーズ補足ダウンロード中: {fname}")
        try:
            import urllib.request
            urllib.request.urlretrieve(url, local)
        except Exception as e:
            print(f"  ⚠ ダウンロード失敗: {e}")
    if os.path.exists(local):
        supp_paths.append(local)

# 既にsupp_dirにあるファイルも追加
for f in os.listdir(supp_dir):
    fp = os.path.join(supp_dir, f)
    if os.path.isfile(fp) and fp not in supp_paths:
        supp_paths.append(fp)

supp_paths = sorted(set(supp_paths))
print(f"補足ファイル数: {len(supp_paths)}")
for p in supp_paths:
    print(f"  {os.path.basename(p)}")

# ================================================================
# Step 2: ファイル種別の特定とAnnDataオブジェクト構築
# ================================================================
print("\n=== Step 2: データのロード ===")

adata = None

# (A) h5ad形式
h5ad_file = find_file(supp_paths, [r"\.h5ad(\.gz)?$"])
if h5ad_file:
    h5ad_file = decompress(h5ad_file)
    print(f"h5adファイルを読み込み: {h5ad_file}")
    adata = sc.read_h5ad(h5ad_file)

# (B) loom形式
if adata is None:
    loom_file = find_file(supp_paths, [r"\.loom(\.gz)?$"])
    if loom_file:
        loom_file = decompress(loom_file)
        print(f"loomファイルを読み込み: {loom_file}")
        adata = sc.read_loom(loom_file)

# (C) 10x Market Exchange (mtx + barcodes + features)
if adata is None:
    mtx_file = find_file(supp_paths, [r"matrix.*\.mtx(\.gz)?$", r"count.*\.mtx(\.gz)?$"])
    bar_file  = find_file(supp_paths, [r"barcode.*\.(tsv|txt)(\.gz)?$", r"cell.*\.(tsv|txt)(\.gz)?$"])
    feat_file = find_file(supp_paths, [r"feature.*\.(tsv|txt)(\.gz)?$", r"gene.*\.(tsv|txt)(\.gz)?$"])
    if mtx_file and bar_file and feat_file:
        mtx_file  = decompress(mtx_file)
        bar_file  = decompress(bar_file)
        feat_file = decompress(feat_file)
        print(f"10x mtx形式を読み込み中...")
        mat = mmread(mtx_file).T.tocsr()
        barcodes = pd.read_csv(bar_file, header=None, sep="\t")[0].tolist()
        features = pd.read_csv(feat_file, header=None, sep="\t")
        gene_names = features.iloc[:, 1] if features.shape[1] > 1 else features.iloc[:, 0]
        import anndata
        adata = anndata.AnnData(X=mat, obs=pd.DataFrame(index=barcodes),
                                var=pd.DataFrame(index=gene_names.tolist()))
        # メタデータ（ある場合）
        meta_file = find_file(supp_paths, [r"meta.*\.(csv|tsv)(\.gz)?$",
                                            r"cluster.*\.(csv|tsv)(\.gz)?$",
                                            r"annot.*\.(csv|tsv)(\.gz)?$"])
        if meta_file:
            meta_file = decompress(meta_file)
            sep = "\t" if meta_file.endswith(".tsv") else ","
            meta_df = pd.read_csv(meta_file, sep=sep, index_col=0)
            common = adata.obs_names.intersection(meta_df.index)
            if len(common) > 0:
                adata = adata[common]
                adata.obs = meta_df.loc[common]

# (D) CSV/TSV形式のカウントマトリクス
if adata is None:
    count_file = find_file(supp_paths, [r"count.*\.(csv|tsv)(\.gz)?$",
                                         r"expr.*\.(csv|tsv)(\.gz)?$",
                                         r"UMI.*\.(csv|tsv)(\.gz)?$"])
    if count_file:
        count_file = decompress(count_file)
        sep = "\t" if count_file.endswith(".tsv") else ","
        print(f"CSV/TSVカウントマトリクスを読み込み中: {os.path.basename(count_file)}")
        count_df = pd.read_csv(count_file, sep=sep, index_col=0)
        # 行=遺伝子, 列=細胞 の想定
        if count_df.shape[0] > count_df.shape[1]:
            # 転置する
            count_df = count_df.T
        import anndata
        adata = anndata.AnnData(X=count_df.values,
                                obs=pd.DataFrame(index=count_df.index),
                                var=pd.DataFrame(index=count_df.columns))

if adata is None:
    raise RuntimeError(
        "認識できる発現データファイルが見つかりませんでした。\n"
        f"ダウンロードされたファイル:\n" + "\n".join(supp_paths)
    )

print(f"AnnDataオブジェクト: {adata.n_obs} 細胞 × {adata.n_vars} 遺伝子")
print(f"メタデータ列: {list(adata.obs.columns)[:10]}")

# ================================================================
# Step 3: 遺伝子の確認
# ================================================================
print("\n=== Step 3: 対象遺伝子の確認 ===")
all_genes = adata.var_names.tolist()
gene_map = {g: find_gene(g, all_genes) for g in GENES}
print("遺伝子名マッピング:")
for k, v in gene_map.items():
    status = "✓" if v else "✗ (未検出)"
    print(f"  {k} → {v or '---'} {status}")

missing = [g for g, v in gene_map.items() if v is None]
if missing:
    print(f"⚠ 未検出遺伝子: {missing}")

if gene_map.get("Gpr176") is None:
    raise RuntimeError("Gpr176が発現行列に見つかりません。解析を中断します。")

# ================================================================
# Step 4: 正規化
# ================================================================
print("\n=== Step 4: 正規化 (log CP10K) ===")
# 生カウントを保存
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("正規化完了")

# ================================================================
# Step 5: Gpr176陽性クラスターの同定
# ================================================================
print("\n=== Step 5: Gpr176陽性クラスターの同定 ===")
gpr176_key = gene_map["Gpr176"]

if hasattr(adata.X, "toarray"):
    gpr176_exp = adata[:, gpr176_key].X.toarray().flatten()
else:
    gpr176_exp = np.array(adata[:, gpr176_key].X).flatten()

adata.obs["Gpr176_expr"]     = gpr176_exp
adata.obs["Gpr176_positive"] = gpr176_exp > GPR176_THR

n_pos = adata.obs["Gpr176_positive"].sum()
n_tot = adata.n_obs
print(f"Gpr176陽性細胞: {n_pos} / {n_tot} ({100*n_pos/n_tot:.1f}%)")

# クラスターアノテーション列を探す
cluster_candidates = ["cluster", "ClusterID", "cell_type", "CellType", "seurat_clusters",
                      "leiden", "louvain", "Cluster", "cluster_label", "subtype",
                      "celltype", "annotation", "SubType"]
cluster_col = next((c for c in cluster_candidates if c in adata.obs.columns), None)

if cluster_col:
    print(f"クラスター列: {cluster_col}")
    cluster_summary = (
        adata.obs.groupby(cluster_col)
        .agg(
            n_cells      =("Gpr176_positive", "count"),
            n_gpr176_pos =("Gpr176_positive", "sum"),
            pct_gpr176   =("Gpr176_positive", lambda x: 100 * x.mean()),
            mean_gpr176  =("Gpr176_expr",     "mean"),
        )
        .reset_index()
        .sort_values("pct_gpr176", ascending=False)
    )
    cluster_summary.to_csv(os.path.join(OUT_DIR, "cluster_Gpr176_summary.csv"), index=False)
    print("クラスター別Gpr176陽性率 (上位10):")
    print(cluster_summary.head(10).to_string(index=False))

    pos_clusters = cluster_summary.loc[
        cluster_summary["pct_gpr176"] > POS_PCT_THR, cluster_col
    ].tolist()

    if not pos_clusters:
        pos_clusters = cluster_summary[cluster_col].iloc[:5].tolist()
        print(f"⚠ 陽性率{POS_PCT_THR}%超クラスターなし。上位5クラスターを使用: {pos_clusters}")
    else:
        print(f"Gpr176陽性クラスター (>{POS_PCT_THR}%): {pos_clusters}")

    adata.obs["Gpr176_cluster"] = np.where(
        adata.obs[cluster_col].isin(pos_clusters), "Gpr176_positive", "Gpr176_negative"
    )
else:
    print("クラスターアノテーション列が見つかりません。細胞単位で陽性/陰性を使用します。")
    cluster_summary = None
    pos_clusters = []
    adata.obs["Gpr176_cluster"] = np.where(
        adata.obs["Gpr176_positive"], "Gpr176_positive", "Gpr176_negative"
    )

# ================================================================
# Step 6: GnazとRGS16の発現量解析
# ================================================================
print("\n=== Step 6: Gnaz・RGS16の発現量解析 ===")

results = []
target_genes = {g: gene_map[g] for g in ["Gnaz", "Rgs16"] if gene_map.get(g)}

pos_mask = adata.obs["Gpr176_cluster"] == "Gpr176_positive"
neg_mask = ~pos_mask

for gname, actual in target_genes.items():
    if hasattr(adata.X, "toarray"):
        expr = adata[:, actual].X.toarray().flatten()
    else:
        expr = np.array(adata[:, actual].X).flatten()

    adata.obs[f"{gname}_expr"] = expr

    pos_vals = expr[pos_mask]
    neg_vals = expr[neg_mask]

    stat, pval = stats.mannwhitneyu(pos_vals, neg_vals, alternative="greater")
    mean_pos = pos_vals.mean()
    mean_neg = neg_vals.mean()
    pct_pos  = 100 * (pos_vals > 0).mean()
    pct_neg  = 100 * (neg_vals > 0).mean()
    log2fc   = np.log2((mean_pos + 1e-6) / (mean_neg + 1e-6))

    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

    results.append({
        "gene":              gname,
        "actual_name":       actual,
        "mean_Gpr176pos":    round(mean_pos, 4),
        "mean_Gpr176neg":    round(mean_neg, 4),
        "pct_expressed_pos": round(pct_pos, 2),
        "pct_expressed_neg": round(pct_neg, 2),
        "log2FC":            round(log2fc, 4),
        "mannwhitney_pval":  float(f"{pval:.4e}"),
        "significance":      sig,
    })

    print(f"\n[{gname} ({actual})]")
    print(f"  Gpr176+クラスター: 平均={mean_pos:.4f}, 発現率={pct_pos:.1f}%")
    print(f"  Gpr176-クラスター: 平均={mean_neg:.4f}, 発現率={pct_neg:.1f}%")
    print(f"  log2FC={log2fc:.3f}, Mann-Whitney p={pval:.3e} {sig}")

result_df = pd.DataFrame(results)
result_df.to_csv(os.path.join(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"), index=False)
print(f"\n結果CSV保存: {OUT_DIR}/Gnaz_RGS16_in_Gpr176clusters.csv")

# ================================================================
# Step 7: 可視化
# ================================================================
print("\n=== Step 7: 可視化 ===")

plot_genes = ["Gpr176"] + list(target_genes.keys())
actual_map = {"Gpr176": gpr176_key, **target_genes}

# --- 7-A: バイオリンプロット ---
n_genes = len(plot_genes)
fig, axes = plt.subplots(1, n_genes, figsize=(4 * n_genes, 5))
if n_genes == 1:
    axes = [axes]

colors = {"Gpr176_positive": "#E64B35", "Gpr176_negative": "#4DBBD5"}

for ax, gname in zip(axes, plot_genes):
    actual = actual_map.get(gname)
    if actual is None:
        ax.set_visible(False)
        continue

    expr_col = f"{gname}_expr"
    if expr_col not in adata.obs.columns:
        if hasattr(adata.X, "toarray"):
            adata.obs[expr_col] = adata[:, actual].X.toarray().flatten()
        else:
            adata.obs[expr_col] = np.array(adata[:, actual].X).flatten()

    for group, color in colors.items():
        vals = adata.obs.loc[adata.obs["Gpr176_cluster"] == group, expr_col]
        parts = ax.violinplot([vals], positions=[list(colors.keys()).index(group)],
                              showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        parts["cmedians"].set_color("black")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Gpr176+\ncluster", "Gpr176-\ncluster"], fontsize=10)
    ax.set_ylabel("log(CP10K + 1)", fontsize=10)
    ax.set_title(f"{gname}\n({actual})", fontsize=11, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.suptitle("Kozareva et al. 2021 (GSE165805)\nGpr176陽性 vs 陰性クラスターでの遺伝子発現",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(OUT_DIR, "violin_Gnaz_RGS16_Gpr176clusters.png"), dpi=150, bbox_inches="tight")
plt.close()
print("バイオリンプロット保存完了")

# --- 7-B: 発現率棒グラフ ---
if len(result_df) > 0:
    bar_df = result_df[["gene", "pct_expressed_pos", "pct_expressed_neg"]].copy()
    bar_df = bar_df.melt(id_vars="gene", var_name="group", value_name="pct")
    bar_df["group"] = bar_df["group"].map({
        "pct_expressed_pos": "Gpr176+ cluster",
        "pct_expressed_neg": "Gpr176- cluster",
    })

    fig, ax = plt.subplots(figsize=(6, 5))
    x = np.arange(len(result_df))
    w = 0.35
    for i, (grp, color) in enumerate([("Gpr176+ cluster", "#E64B35"), ("Gpr176- cluster", "#4DBBD5")]):
        vals = bar_df.loc[bar_df["group"] == grp, "pct"].tolist()
        bars = ax.bar(x + i * w, vals, w, label=grp, color=color, alpha=0.85, edgecolor="black", linewidth=0.5)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    f"{v:.1f}%", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x + w / 2)
    ax.set_xticklabels(result_df["gene"].tolist(), fontsize=12)
    ax.set_ylabel("発現細胞率 (%)", fontsize=11)
    ax.set_title("Gpr176陽性クラスターにおけるGnaz・RGS16発現率\n(Kozareva et al. 2021, GSE165805)", fontsize=11)
    ax.legend(fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # p値の注釈
    for i, row in enumerate(results):
        sig = row["significance"]
        ymax = max(bar_df.loc[bar_df["gene"] == row["gene"], "pct"]) + 4
        ax.text(x[i] + w / 2, ymax, sig, ha="center", fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "barplot_expression_pct.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "barplot_expression_pct.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("棒グラフ保存完了")

# --- 7-C: ドットプロット（クラスター × 遺伝子）---
if cluster_col and cluster_summary is not None:
    show_clusters = (
        cluster_summary.head(10)[cluster_col].tolist()
    )
    dot_genes_actual = [actual_map[g] for g in plot_genes if actual_map.get(g)]
    obs_sub = adata[adata.obs[cluster_col].isin(show_clusters)].obs.copy()

    # 平均発現量と発現率を計算
    dot_data = []
    for cl in show_clusters:
        cl_mask = obs_sub[cluster_col] == cl
        for gname in plot_genes:
            actual = actual_map.get(gname)
            if actual is None:
                continue
            col = f"{gname}_expr"
            if col not in obs_sub.columns:
                continue
            vals = obs_sub.loc[cl_mask, col]
            dot_data.append({
                "cluster": str(cl),
                "gene":    gname,
                "mean":    vals.mean(),
                "pct":     100 * (vals > 0).mean(),
            })

    dot_df = pd.DataFrame(dot_data)
    if len(dot_df) > 0:
        pivot_mean = dot_df.pivot(index="cluster", columns="gene", values="mean")
        pivot_pct  = dot_df.pivot(index="cluster", columns="gene", values="pct")

        fig, ax = plt.subplots(figsize=(3 * len(plot_genes), max(6, len(show_clusters) * 0.5)))
        for j, gene in enumerate(plot_genes):
            if gene not in pivot_mean.columns:
                continue
            for i, cl in enumerate(pivot_mean.index):
                size  = pivot_pct.loc[cl, gene] * 3   # 発現率をドットサイズに
                color_val = pivot_mean.loc[cl, gene]
                ax.scatter(j, i, s=size, c=[[color_val]], cmap="Reds",
                           vmin=0, vmax=pivot_mean.max().max(), edgecolors="gray", linewidths=0.3)

        ax.set_xticks(range(len(plot_genes)))
        ax.set_xticklabels(plot_genes, fontsize=11)
        ax.set_yticks(range(len(pivot_mean.index)))
        ax.set_yticklabels(pivot_mean.index, fontsize=8)
        ax.set_xlabel("遺伝子", fontsize=11)
        ax.set_ylabel("クラスター", fontsize=11)
        ax.set_title("ドットプロット（サイズ=発現率%, 色=平均発現量）\n上位10クラスター", fontsize=11)

        # カラーバー
        sm = plt.cm.ScalarMappable(cmap="Reds",
                                    norm=plt.Normalize(0, pivot_mean.max().max()))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label="平均発現量 log(CP10K+1)", shrink=0.5)

        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.pdf"), bbox_inches="tight")
        plt.savefig(os.path.join(OUT_DIR, "dotplot_Gpr176_Gnaz_RGS16.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print("ドットプロット保存完了")

# ================================================================
# 最終サマリー
# ================================================================
print("\n" + "=" * 50)
print("         解析完了サマリー")
print("=" * 50)
print(f"データセット : {GEO_ID} (Kozareva et al. 2021)")
print(f"総細胞数     : {n_tot:,}")
print(f"Gpr176陽性   : {int(n_pos):,} 細胞 ({100*n_pos/n_tot:.1f}%)")
print()
print("Gpr176陽性クラスターでの発現:")
for r in results:
    print(f"  {r['gene']:6s}: 発現率 {r['pct_expressed_pos']:.1f}% "
          f"(vs {r['pct_expressed_neg']:.1f}%), "
          f"log2FC={r['log2FC']:.2f}, "
          f"p={r['mannwhitney_pval']:.2e} {r['significance']}")
print()
print(f"出力フォルダ : {OUT_DIR}/")
print("  - cluster_Gpr176_summary.csv   (クラスター別Gpr176陽性率)")
print("  - Gnaz_RGS16_in_Gpr176clusters.csv  (統計検定結果)")
print("  - violin_Gnaz_RGS16_Gpr176clusters.pdf/png")
print("  - barplot_expression_pct.pdf/png")
print("  - dotplot_Gpr176_Gnaz_RGS16.pdf/png")
print("=" * 50)
