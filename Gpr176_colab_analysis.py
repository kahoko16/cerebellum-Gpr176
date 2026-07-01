"""
Gpr176_colab_analysis.py
Google Colab用

論文: Kozareva et al. 2021, Nature
  "A transcriptomic atlas of mouse cerebellar cortex
   comprehensively defines cell types"
GEO: GSE165371

目的: Gpr176陽性細胞にGnaz（Gz）とRGS16が共発現するかを検証
     （全細胞 + プルキンエ細胞での解析）

【Colabでの実行手順】
1. Google DriveにGSE165371_cb_adult_mouse.tar.gzをアップロード
2. このスクリプトをColabにアップロードし、セルごとに実行
3. 結果はGoogle Driveに保存される
"""

# ==============================================================
# Cell 1: Google Driveのマウント
# ==============================================================
from google.colab import drive
drive.mount('/content/drive')

# ==============================================================
# Cell 2: パッケージインストール
# ==============================================================
# !pip install -q scanpy anndata matplotlib seaborn scipy

# ==============================================================
# Cell 3: 設定（ここだけ変更してください）
# ==============================================================
import os

# Google Drive内のファイルパス
# 例: マイドライブ直下に置いた場合
DRIVE_DIR  = "/content/drive/MyDrive"          # Driveのルート
DATA_FILE  = "GSE165371_cb_adult_mouse.tar.gz" # tarファイル名
# または解凍済みの場合（ファイルをそのままDriveに置いた場合）
MTX_FILE      = "cb_adult_mouse.mtx.gz"
BARCODES_FILE = "cb_adult_mouse_barcodes.txt"
GENES_FILE    = "cb_adult_mouse_genes.txt"

GENES_OF_INT = ["Gpr176", "Gnaz", "Rgs16"]
OUT_DIR      = os.path.join(DRIVE_DIR, "Gpr176_cluster_results")
os.makedirs(OUT_DIR, exist_ok=True)

# ==============================================================
# Cell 4: ファイルの準備（解凍）
# ==============================================================
import tarfile, gzip, shutil

data_dir = "/content/cb_data"
os.makedirs(data_dir, exist_ok=True)

tar_path = os.path.join(DRIVE_DIR, DATA_FILE)
mtx_path = os.path.join(DRIVE_DIR, MTX_FILE)

if os.path.exists(tar_path):
    print(f"解凍中: {DATA_FILE}")
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(data_dir)
    print("解凍完了")
    # 解凍されたファイルを探す
    for root, dirs, files in os.walk(data_dir):
        for f in files:
            print(f"  {f}")
elif os.path.exists(mtx_path):
    print("解凍済みファイルを使用")
    data_dir = DRIVE_DIR
else:
    raise FileNotFoundError(
        f"ファイルが見つかりません。\n"
        f"Google Driveの {DRIVE_DIR} に以下のどちらかを置いてください:\n"
        f"  - {DATA_FILE}（tar.gz）\n"
        f"  または\n"
        f"  - {MTX_FILE}, {BARCODES_FILE}, {GENES_FILE}（解凍済み）"
    )

# ==============================================================
# Cell 5: ファイルパスの自動検出
# ==============================================================
import re

def find_file(directory, patterns):
    for root, dirs, files in os.walk(directory):
        for f in files:
            for pat in patterns:
                if re.search(pat, f, re.IGNORECASE):
                    return os.path.join(root, f)
    return None

mtx_path  = find_file(data_dir, [r"\.mtx(\.gz)?$"])
bar_path  = find_file(data_dir, [r"barcodes(\.(txt|tsv))?$",
                                  r"barcodes(\.txt|\.tsv)?(\.gz)?$"])
gene_path = find_file(data_dir, [r"genes(\.(txt|tsv))?$",
                                  r"features(\.(txt|tsv))?$",
                                  r"genes(\.txt|\.tsv)?(\.gz)?$"])

print(f"mtx      : {mtx_path}")
print(f"barcodes : {bar_path}")
print(f"genes    : {gene_path}")

# ==============================================================
# Cell 6: 遺伝子リストとバーコードの読み込み（軽い）
# ==============================================================
import pandas as pd
import numpy as np

def read_txt(path):
    open_fn = gzip.open if path.endswith(".gz") else open
    with open_fn(path, "rt") as f:
        return [line.strip() for line in f]

genes    = read_txt(gene_path)
barcodes = read_txt(bar_path)

# タブ区切りの場合（features.tsvは2列目が遺伝子名）
if "\t" in genes[0]:
    genes = [g.split("\t")[1] if len(g.split("\t")) > 1
             else g.split("\t")[0] for g in genes]

print(f"総遺伝子数: {len(genes)}")
print(f"総細胞数  : {len(barcodes)}")

# 対象遺伝子のインデックスを特定（1-indexed, MTX形式に合わせる）
gene_lower = [g.lower() for g in genes]
target_idx = {}  # {gene_name: 1-based index}
for g in GENES_OF_INT:
    if g.lower() in gene_lower:
        idx = gene_lower.index(g.lower()) + 1  # 1-indexed
        actual = genes[idx - 1]
        target_idx[g] = (idx, actual)
        print(f"  {g} → {actual} (行 {idx})")
    else:
        print(f"  {g} → 未検出")

if "Gpr176" not in target_idx:
    raise ValueError("Gpr176が遺伝子リストに見つかりません")

# ==============================================================
# Cell 7: MTXファイルをスキャンして対象遺伝子だけ抽出
# ==============================================================
from scipy.sparse import csr_matrix

print("MTXファイルをスキャン中（数分かかります）...")

open_fn = gzip.open if mtx_path.endswith(".gz") else open

# 行(gene_idx) → {col: val} の辞書
data_store = {g: {"rows": [], "cols": [], "vals": []}
              for g in target_idx}

with open_fn(mtx_path, "rt") as f:
    # ヘッダーをスキップ
    for line in f:
        if not line.startswith("%"):
            n_genes_mtx, n_cells_mtx, nnz = map(int, line.strip().split())
            break

    print(f"  行列サイズ: {n_genes_mtx} × {n_cells_mtx}, nnz={nnz:,}")
    target_set = {info[0]: gname for gname, info in target_idx.items()}

    count = 0
    for line in f:
        row, col, val = line.strip().split()
        row = int(row)
        if row in target_set:
            gname = target_set[row]
            data_store[gname]["cols"].append(int(col) - 1)  # 0-indexed
            data_store[gname]["vals"].append(float(val))
        count += 1
        if count % 5_000_000 == 0:
            print(f"  {100*count/nnz:.1f}% スキャン済み...")

print("スキャン完了")

# スパース行列に変換（遺伝子 × 細胞）
n_cells = len(barcodes)
expr_dict = {}
for gname, info in target_idx.items():
    d = data_store[gname]
    mat = csr_matrix(
        (d["vals"], ([0] * len(d["cols"]), d["cols"])),
        shape=(1, n_cells)
    )
    expr_dict[gname] = np.asarray(mat.todense()).flatten()
    print(f"  {gname}: 発現細胞数 {(expr_dict[gname] > 0).sum():,} / {n_cells:,}")

# ==============================================================
# Cell 8: Gpr176陽性細胞の定義
# ==============================================================
gpr176_expr = expr_dict["Gpr176"]
gpr176_pos  = gpr176_expr > 0

n_pos = gpr176_pos.sum()
print(f"Gpr176陽性細胞: {n_pos:,} / {n_cells:,} ({100*n_pos/n_cells:.1f}%)")

meta = pd.DataFrame({
    "cell":             barcodes,
    "Gpr176_expr":      gpr176_expr,
    "Gpr176_positive":  gpr176_pos,
    "Gpr176_cluster":   np.where(gpr176_pos, "Gpr176_positive", "Gpr176_negative"),
})
for gname in ["Gnaz", "Rgs16"]:
    if gname in expr_dict:
        meta[f"{gname}_expr"] = expr_dict[gname]

# ==============================================================
# Cell 9: 統計検定（Mann-Whitney U）
# ==============================================================
from scipy import stats

results = []
for gname in ["Gnaz", "Rgs16"]:
    if gname not in expr_dict:
        print(f"{gname}: 未検出のためスキップ")
        continue

    expr = expr_dict[gname]
    pos_vals = expr[gpr176_pos]
    neg_vals = expr[~gpr176_pos]

    stat, pval = stats.mannwhitneyu(pos_vals, neg_vals, alternative="greater")
    mean_pos = pos_vals.mean()
    mean_neg = neg_vals.mean()
    pct_pos  = 100 * (pos_vals > 0).mean()
    pct_neg  = 100 * (neg_vals > 0).mean()
    log2fc   = np.log2((mean_pos + 1e-6) / (mean_neg + 1e-6))
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

    results.append({
        "gene": gname,
        "mean_Gpr176pos": round(mean_pos, 4),
        "mean_Gpr176neg": round(mean_neg, 4),
        "pct_expressed_pos": round(pct_pos, 2),
        "pct_expressed_neg": round(pct_neg, 2),
        "log2FC": round(log2fc, 4),
        "mannwhitney_pval": float(f"{pval:.4e}"),
        "significance": sig,
    })

    print(f"\n[{gname}]")
    print(f"  Gpr176+: 平均={mean_pos:.4f}, 発現率={pct_pos:.1f}%")
    print(f"  Gpr176-: 平均={mean_neg:.4f}, 発現率={pct_neg:.1f}%")
    print(f"  log2FC={log2fc:.3f}, p={pval:.3e} {sig}")

result_df = pd.DataFrame(results)
result_df.to_csv(os.path.join(OUT_DIR, "Gnaz_RGS16_in_Gpr176clusters.csv"), index=False)
print(f"\n結果CSV保存: {OUT_DIR}/Gnaz_RGS16_in_Gpr176clusters.csv")

# ==============================================================
# Cell 10: 可視化
# ==============================================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

colors = {"Gpr176_positive": "#E64B35", "Gpr176_negative": "#4DBBD5"}
plot_genes = ["Gpr176"] + [g for g in ["Gnaz", "Rgs16"] if g in expr_dict]

# --- バイオリンプロット ---
fig, axes = plt.subplots(1, len(plot_genes), figsize=(4 * len(plot_genes), 5))
if len(plot_genes) == 1:
    axes = [axes]

for ax, gname in zip(axes, plot_genes):
    plot_data = pd.DataFrame({
        "expr":  expr_dict[gname],
        "group": meta["Gpr176_cluster"],
    })
    sns.violinplot(data=plot_data, x="group", y="expr",
                   palette=colors, ax=ax, cut=0, inner="box")
    ax.set_title(gname, fontsize=12, fontweight="bold", fontstyle="italic")
    ax.set_xlabel("")
    ax.set_ylabel("UMI count")
    ax.set_xticklabels(["Gpr176+", "Gpr176-"], fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.suptitle("Kozareva et al. 2021 (GSE165371)\nGpr176陽性 vs 陰性細胞での遺伝子発現",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "violin_Gnaz_RGS16_Gpr176.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(OUT_DIR, "violin_Gnaz_RGS16_Gpr176.png"), dpi=150, bbox_inches="tight")
plt.show()
print("バイオリンプロット保存完了")

# --- 発現率棒グラフ ---
fig, ax = plt.subplots(figsize=(6, 5))
x    = np.arange(len(result_df))
w    = 0.35
for i, (grp, color) in enumerate([("pct_expressed_pos", "#E64B35"),
                                   ("pct_expressed_neg", "#4DBBD5")]):
    label = "Gpr176+ cells" if "pos" in grp else "Gpr176- cells"
    vals  = result_df[grp].tolist()
    bars  = ax.bar(x + i * w, vals, w, label=label,
                   color=color, alpha=0.85, edgecolor="black", linewidth=0.5)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5, f"{v:.1f}%",
                ha="center", va="bottom", fontsize=9)

for i, row in result_df.iterrows():
    ymax = max(row["pct_expressed_pos"], row["pct_expressed_neg"]) + 4
    ax.text(x[i] + w / 2, ymax, row["significance"],
            ha="center", fontsize=13, fontweight="bold")

ax.set_xticks(x + w / 2)
ax.set_xticklabels(result_df["gene"].tolist(), fontsize=12, fontstyle="italic")
ax.set_ylabel("発現細胞率 (%)", fontsize=11)
ax.set_title("Gpr176陽性細胞におけるGnaz・RGS16の発現率\nKozareva et al. 2021 (GSE165371)",
             fontsize=11)
ax.legend(fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "barplot_expression_pct.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(OUT_DIR, "barplot_expression_pct.png"), dpi=150, bbox_inches="tight")
plt.show()
print("棒グラフ保存完了")

# ==============================================================
# Cell 11: メタデータの読み込みとバーコード照合
# ==============================================================
META_FILE = os.path.join(DRIVE_DIR, "full_cb_metadata.csv")

full_meta = pd.read_csv(META_FILE, index_col=0)  # index = "IXa_M003_XXXX..." 形式
print(f"メタデータ: {len(full_meta):,} 細胞")
print(f"細胞タイプ分布:")
print(full_meta["final_annotation_cluster"].value_counts().head(20))

# バーコード照合
# メタデータ: "IXa_M003_CTTCCTTTCCGTATGA" → 最後の "_" 以降がバーコード本体
# MTXバーコード: "CTTCCTTTCCGTATGA-1" または "CTTCCTTTCCGTATGA" 形式
meta_bc      = full_meta.index.tolist()
meta_bc_trim = [b.rsplit("_", 1)[-1] for b in meta_bc]  # バーコード本体を抽出

# MTXバーコードから "-1" サフィックスを除去して照合
mtx_bc_trim  = [b.split("-")[0] for b in barcodes]

# 照合テーブルを作成
bc_map = pd.DataFrame({
    "mtx_barcode":  barcodes,
    "mtx_bc_trim":  mtx_bc_trim,
})
meta_df = full_meta.copy()
meta_df["meta_bc_trim"] = meta_bc_trim

merged = bc_map.merge(
    meta_df.reset_index().rename(columns={"index": "meta_barcode"}),
    left_on  = "mtx_bc_trim",
    right_on = "meta_bc_trim",
    how      = "left"
)
merged.index = range(len(merged))

n_matched = merged["final_annotation_cluster"].notna().sum()
print(f"\nMTX細胞数       : {len(barcodes):,}")
print(f"メタデータ一致数 : {n_matched:,} ({100*n_matched/len(barcodes):.1f}%)")

# 一致率が低い場合は別の照合方法を試みる
if n_matched / len(barcodes) < 0.5:
    print("⚠ 一致率が低い → サンプルIDなしで再照合を試みます")
    # メタデータのバーコードをそのまま使う
    merged2 = bc_map.merge(
        meta_df.reset_index().rename(columns={"index": "meta_barcode"}),
        left_on  = "mtx_barcode",
        right_on = "meta_barcode",
        how      = "left"
    )
    n_matched2 = merged2["final_annotation_cluster"].notna().sum()
    print(f"再照合一致数: {n_matched2:,}")
    if n_matched2 > n_matched:
        merged = merged2

# 発現データとメタデータを統合
cell_annotation = merged["final_annotation_cluster"].fillna("Unknown").tolist()
cell_subtype    = merged["final_annotation_subcluster"].fillna("Unknown").tolist()
cell_region     = merged["regions"].fillna("Unknown").tolist() if "regions" in merged.columns else ["Unknown"] * len(barcodes)

meta["cell_type"]    = cell_annotation
meta["cell_subtype"] = cell_subtype
meta["region"]       = cell_region

# プルキンエ細胞フラグ
meta["is_purkinje"] = meta["cell_type"].str.contains("Purkinje", case=False, na=False) & \
                      (meta["cell_type"] != "REMOVED")

print(f"\nプルキンエ細胞数: {meta['is_purkinje'].sum():,}")
print(f"細胞タイプ別Gpr176陽性率:")
for ct in meta["cell_type"].value_counts().head(10).index:
    ct_mask = (meta["cell_type"] == ct) & (meta["cell_type"] != "REMOVED")
    pct = 100 * meta.loc[ct_mask, "Gpr176_positive"].mean()
    n   = ct_mask.sum()
    print(f"  {ct:30s}: {pct:5.1f}% (n={n:,})")

# ==============================================================
# Cell 12: ドットプロット（細胞タイプ × 遺伝子）
# ==============================================================
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# REMOVED を除外
valid_mask  = meta["cell_type"] != "REMOVED"
valid_meta  = meta[valid_mask].copy()

# 表示する細胞タイプ（Gpr176陽性率上位15 + Purkinjeを必ず含む）
ct_gpr176 = (
    valid_meta.groupby("cell_type")["Gpr176_positive"]
    .agg(["mean", "count"])
    .rename(columns={"mean": "pct_gpr176", "count": "n_cells"})
    .query("n_cells >= 50")
    .sort_values("pct_gpr176", ascending=False)
)
purkinje_cts = [c for c in ct_gpr176.index if "Purkinje" in c]
other_cts    = [c for c in ct_gpr176.head(15).index if c not in purkinje_cts]
show_cts     = purkinje_cts + other_cts
show_cts     = list(dict.fromkeys(show_cts))[:20]  # 重複除去・上位20

print(f"ドットプロット対象クラスター数: {len(show_cts)}")

# 各クラスター × 遺伝子の平均発現量と発現率を計算
dot_data = []
for ct in show_cts:
    ct_mask = valid_meta["cell_type"] == ct
    ct_cells_idx = valid_meta.index[ct_mask]
    for gname in plot_genes:
        vals = expr_dict[gname][ct_cells_idx]
        dot_data.append({
            "cell_type": ct,
            "gene":      gname,
            "mean_expr": vals.mean(),
            "pct_expr":  100 * (vals > 0).mean(),
            "n_cells":   ct_mask.sum(),
        })

dot_df = pd.DataFrame(dot_data)

# Gpr176陽性率でクラスターを並べ替え
ct_order = ct_gpr176.loc[ct_gpr176.index.isin(show_cts)].sort_values("pct_gpr176").index.tolist()
dot_df["cell_type"] = pd.Categorical(dot_df["cell_type"], categories=ct_order, ordered=True)
dot_df = dot_df.sort_values("cell_type")

# プロット
fig, ax = plt.subplots(figsize=(len(plot_genes) * 1.8 + 2, len(show_cts) * 0.45 + 2))

for j, gene in enumerate(plot_genes):
    gene_df = dot_df[dot_df["gene"] == gene]
    for i, ct in enumerate(ct_order):
        row = gene_df[gene_df["cell_type"] == ct]
        if len(row) == 0:
            continue
        size  = row["pct_expr"].values[0] * 4   # 発現率 → ドットサイズ
        color = row["mean_expr"].values[0]
        ax.scatter(j, i, s=size, c=[[color]], cmap="Reds",
                   vmin=0, vmax=dot_df["mean_expr"].max(),
                   edgecolors="gray", linewidths=0.3)

ax.set_xticks(range(len(plot_genes)))
ax.set_xticklabels(plot_genes, fontsize=11, fontstyle="italic")
ax.set_yticks(range(len(ct_order)))
ax.set_yticklabels(ct_order, fontsize=8)
ax.set_xlabel("遺伝子", fontsize=11)
ax.set_ylabel("細胞タイプ（Gpr176陽性率昇順）", fontsize=10)
ax.set_title("Gpr176 / Gnaz / Rgs16 発現（細胞タイプ別）\nKozareva et al. 2021 (GSE165371)",
             fontsize=11)

# カラーバー（平均発現量）
sm = plt.cm.ScalarMappable(cmap="Reds",
     norm=plt.Normalize(0, dot_df["mean_expr"].max()))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
cbar.set_label("平均UMIカウント", fontsize=9)

# サイズ凡例
for pct_val in [10, 30, 50]:
    ax.scatter([], [], s=pct_val * 4, c="gray", alpha=0.5,
               label=f"{pct_val}%")
ax.legend(title="発現細胞率", loc="lower right", fontsize=8, title_fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "dotplot_celltype_Gpr176_Gnaz_RGS16.pdf"), bbox_inches="tight")
plt.savefig(os.path.join(OUT_DIR, "dotplot_celltype_Gpr176_Gnaz_RGS16.png"), dpi=150, bbox_inches="tight")
plt.show()
print("ドットプロット保存完了")

# ==============================================================
# Cell 13: プルキンエ細胞のみの詳細解析
# ==============================================================
purk_mask = meta["is_purkinje"].values
n_purk    = purk_mask.sum()
print(f"\nプルキンエ細胞解析 (n={n_purk:,})")

if n_purk > 0:
    purk_results = []
    for gname in ["Gnaz", "Rgs16"]:
        if gname not in expr_dict:
            continue
        expr    = expr_dict[gname]
        gpr_pos = gpr176_pos & purk_mask
        gpr_neg = (~gpr176_pos) & purk_mask

        pos_vals = expr[gpr_pos]
        neg_vals = expr[gpr_neg]

        if len(pos_vals) == 0 or len(neg_vals) == 0:
            continue

        stat, pval = stats.mannwhitneyu(pos_vals, neg_vals, alternative="greater")
        pct_pos = 100 * (pos_vals > 0).mean()
        pct_neg = 100 * (neg_vals > 0).mean()
        log2fc  = np.log2((pos_vals.mean() + 1e-6) / (neg_vals.mean() + 1e-6))
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"

        purk_results.append({
            "gene": gname,
            "n_Gpr176pos_purkinje": int(gpr_pos.sum()),
            "n_Gpr176neg_purkinje": int(gpr_neg.sum()),
            "pct_expressed_pos": round(pct_pos, 2),
            "pct_expressed_neg": round(pct_neg, 2),
            "log2FC": round(log2fc, 4),
            "mannwhitney_pval": float(f"{pval:.4e}"),
            "significance": sig,
        })
        print(f"  [{gname}] プルキンエ内 Gpr176+: {pct_pos:.1f}% vs Gpr176-: {pct_neg:.1f}%  {sig}")

    purk_df = pd.DataFrame(purk_results)
    purk_df.to_csv(os.path.join(OUT_DIR, "purkinje_Gnaz_RGS16_results.csv"), index=False)

    # プルキンエ細胞バイオリンプロット
    fig, axes = plt.subplots(1, len(plot_genes), figsize=(4 * len(plot_genes), 5))
    if len(plot_genes) == 1:
        axes = [axes]

    for ax, gname in zip(axes, plot_genes):
        plot_data = pd.DataFrame({
            "expr":  expr_dict[gname][purk_mask],
            "group": np.where(gpr176_pos[purk_mask], "Gpr176+", "Gpr176-"),
        })
        sns.violinplot(data=plot_data, x="group", y="expr",
                       palette={"Gpr176+": "#E64B35", "Gpr176-": "#4DBBD5"},
                       ax=ax, cut=0, inner="box")
        ax.set_title(gname, fontsize=12, fontweight="bold", fontstyle="italic")
        ax.set_xlabel("")
        ax.set_ylabel("UMI count")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.suptitle(f"プルキンエ細胞のみ (n={n_purk:,})\nGpr176陽性 vs 陰性での発現",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "violin_purkinje_Gpr176_Gnaz_RGS16.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "violin_purkinje_Gpr176_Gnaz_RGS16.png"), dpi=150, bbox_inches="tight")
    plt.show()
    print("プルキンエ細胞バイオリンプロット保存完了")

# ==============================================================
# Cell 14: 最終サマリー
# ==============================================================
print("\n" + "=" * 50)
print("         解析完了サマリー")
print("=" * 50)
print(f"データセット : GSE165371 (Kozareva et al. 2021)")
print(f"総細胞数     : {n_cells:,}")
print(f"Gpr176陽性   : {int(n_pos):,} 細胞 ({100*n_pos/n_cells:.1f}%)")
print()
print("Gpr176陽性細胞での発現:")
for r in results:
    print(f"  {r['gene']:6s}: 発現率 {r['pct_expressed_pos']:.1f}%"
          f" (vs {r['pct_expressed_neg']:.1f}%),"
          f" log2FC={r['log2FC']:.2f},"
          f" p={r['mannwhitney_pval']:.2e} {r['significance']}")
print()
print(f"出力フォルダ: {OUT_DIR}/")
print("  Gnaz_RGS16_in_Gpr176clusters.csv")
print("  violin_Gnaz_RGS16_Gpr176.pdf/png")
print("  barplot_expression_pct.pdf/png")
print("=" * 50)
