#!/usr/bin/env python3
"""
Immune-Focused NMF Niche Analysis (Tier 2) — Option A
======================================================
Drops Tumor and Epithelial columns from Cell2Location abundances,
runs NMF on immune + fibroblast cell types only to discover
immune-centric niches. Then cross-references with CT/IM/Stroma
regions and compares niche proportions across metastasis groups.

This is Step 2 of the niche-driven analysis framework for Aim 2.

Usage: Run cells sequentially in Jupyter or as a script.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations
import os
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

META_PATH = "/Users/adiallo/Desktop/Thesis/Data_Documents/data_all.csv"
SAMPLE_DIR = "/Users/adiallo/Desktop/Thesis/Data_Documents/All_Data/Visium_samples/Visium_Data/samples"
OUTPUT_DIR = "/Users/adiallo/Desktop/Thesis/Results/AIM_2_Results/Immune_Niche_Analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Abundance key (q05 = 5% quantile from Cell2Location)
ABUNDANCE_KEY = 'q05_cell_abundance_w_sf'

# Region definition (same as Tier 1)
TUMOR_LABELS = ['cancer_cancer']
STROMA_LABELS = ['benign_stroma', 'benign_fat', 'inter', 'stroma', 'benign_muscularis']
BOUNDARY_ROWS = 5  # 500µm ≈ 5 Visium spot spacings

# NMF parameters
K_RANGE = range(3, 16)        # sweep range for immune-focused (fewer types = fewer niches)
N_RESTARTS = 10                # more restarts for stability given smaller matrix
RANDOM_STATE = 42

# ==============================================================================
# 2. CELL TYPE CLASSIFICATION
#    - Matches the taxonomy used in CellChat and DE pipelines
#    - Returns (group, subgroup) for fine-grained immune niche discovery
# ==============================================================================

def classify_c2l_column(col_name):
    """
    Classify a Cell2Location column into biological category.
    Returns (category, keep_for_immune_nmf)

    Categories kept for immune NMF:
      - CD8_T, CD8_T_exhausted, CD4_T, Treg, gdT, NK, ILC, Other_T
      - B_cell, Plasma
      - Macrophage, DC, DC2, Granulocyte, Mast
      - CAF, Fibroblast, Myofibroblast, Pericyte, Smooth_muscle

    Categories DROPPED (Option A):
      - Tumor, Epithelial, Endothelial, Schwann
    """
    clean = col_name.replace("meanscell_abundance_w_sf_", "").replace("q05_cell_abundance_w_sf_", "")

    # --- TUMOR (DROP) ---
    if clean.startswith("Tumor"):
        return ("Tumor", False)

    # --- EPITHELIAL (DROP) ---
    if clean.startswith("cE"):
        return ("Epithelial", False)

    # --- STRUCTURAL ---
    if clean.startswith("cS"):
        import re
        num_match = re.search(r'\d+', clean)
        if num_match:
            num = int(num_match.group())
            if 1 <= num <= 14:
                return ("Endothelial", False)  # DROP — not relevant to immune question
            if num == 32:
                return ("Smooth_muscle", True)
            if num in (27, 28, 29, 30, 31):
                return ("CAF", True)
            if num == 26:
                return ("Myofibroblast", True)
            if num in (22, 24, 25):
                return ("Fibroblast", True)
            if num == 23:
                return ("Fibroblast", True)
            if 15 <= num <= 20:
                return ("Pericyte", True)
            if num == 33:
                return ("Schwann", False)  # DROP
        return ("Fibroblast", True)

    # --- B CELLS & PLASMA ---
    if clean.startswith("cB"):
        return ("B_cell", True)
    if clean.startswith("cP"):
        return ("Plasma", True)

    # --- MAST ---
    if clean.startswith("cMA"):
        return ("Mast", True)

    # --- MYELOID ---
    if clean.startswith("cM"):
        if "cM10" in clean:
            return ("Granulocyte", True)
        if "cM01" in clean or "cM02" in clean:
            return ("Macrophage", True)
        return ("DC", True)

    # --- T CELLS & NK ---
    if clean.startswith("cTNI"):
        import re
        num_match = re.search(r'\d+', clean)
        if num_match:
            num = int(num_match.group())
            if 1 <= num <= 7:
                return ("CD4_T", True)
            if num in (8, 9):
                return ("Treg", True)
            if 10 <= num <= 13:
                return ("CD8_T", True)
            if 14 <= num <= 16:
                return ("CD8_T_exhausted", True)
            if num in (17, 18):
                return ("gdT", True)
            if 19 <= num <= 22:
                return ("Other_T", True)
            if 23 <= num <= 25:
                return ("NK", True)
            if num == 26:
                return ("ILC", True)
        return ("Other_T", True)

    return ("Other", False)


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

def get_mets_group(row):
    """Assign metastasis group from metadata."""
    dist = str(row.get('Distant_Mets', 'False')).upper()
    ln = str(row.get('ln_only', 'False')).upper()
    if dist in ['TRUE', 'T']:
        return 'Distant_Mets'
    elif ln in ['TRUE', 'T']:
        return 'LN_Mets'
    return 'No_Mets'


def compute_regions(coords, histology, tumor_labels, stroma_labels, boundary_rows=5):
    """Compute CT/IM/Stroma regions using boundary distance method."""
    is_tumor = np.isin(histology, tumor_labels)
    is_stroma = np.isin(histology, stroma_labels)

    nn = NearestNeighbors(n_neighbors=2).fit(coords)
    dists, _ = nn.kneighbors(coords)
    spot_spacing = np.median(dists[:, 1])
    boundary_dist = boundary_rows * spot_spacing

    tumor_coords = coords[is_tumor]
    stroma_coords = coords[is_stroma]

    region = np.full(len(histology), 'Other', dtype=object)

    if len(tumor_coords) > 0 and len(stroma_coords) > 0:
        nn_tumor = NearestNeighbors(n_neighbors=1).fit(tumor_coords)
        dist_to_tumor, _ = nn_tumor.kneighbors(coords)
        dist_to_tumor = dist_to_tumor.ravel()

        nn_stroma = NearestNeighbors(n_neighbors=1).fit(stroma_coords)
        dist_to_stroma, _ = nn_stroma.kneighbors(coords)
        dist_to_stroma = dist_to_stroma.ravel()

        for i in range(len(histology)):
            if is_tumor[i]:
                region[i] = 'IM' if dist_to_stroma[i] <= boundary_dist else 'CT'
            elif is_stroma[i]:
                region[i] = 'IM' if dist_to_tumor[i] <= boundary_dist else 'Stroma'
    elif len(tumor_coords) > 0:
        region[is_tumor] = 'CT'
    elif len(stroma_coords) > 0:
        region[is_stroma] = 'Stroma'

    return region


# ==============================================================================
# 4. LOAD ALL SAMPLES & BUILD IMMUNE-FOCUSED MATRIX
# ==============================================================================

print("=" * 70)
print("STEP 1: Loading samples and building immune-focused abundance matrix")
print("=" * 70)

meta = pd.read_csv(META_PATH)
meta['deident'] = meta['deident'].astype(str)
meta['mets_group'] = meta.apply(get_mets_group, axis=1)

all_abundances_full = []   # full abundance (for DE later)
all_abundances_immune = [] # immune + fibroblast only
all_coords = []
all_histology = []
all_regions = []
all_sample_ids = []
all_mets = []
sample_spot_counts = {}
column_classification = None

for _, row in meta.iterrows():
    sid = str(row['deident'])
    path = os.path.join(SAMPLE_DIR, f"{sid}.h5ad")
    if not os.path.exists(path):
        continue

    try:
        adata = sc.read_h5ad(path)

        if ABUNDANCE_KEY not in adata.obsm:
            print(f"  {sid}: missing {ABUNDANCE_KEY}, skipping")
            continue
        if 'histology' not in adata.obs.columns:
            print(f"  {sid}: missing histology, skipping")
            continue

        abund = adata.obsm[ABUNDANCE_KEY]
        coords = adata.obsm['spatial']
        histology = adata.obs['histology'].values

        # Classify columns on first sample
        if column_classification is None:
            column_classification = {}
            for col in abund.columns:
                cat, keep = classify_c2l_column(col)
                column_classification[col] = {'category': cat, 'keep': keep}

            immune_cols = [c for c, v in column_classification.items() if v['keep']]
            dropped_cols = [c for c, v in column_classification.items() if not v['keep']]

            print(f"\n  Column classification:")
            print(f"    Total C2L columns: {len(column_classification)}")
            print(f"    Kept for immune NMF: {len(immune_cols)}")
            print(f"    Dropped: {len(dropped_cols)}")
            print(f"\n  Kept categories:")
            kept_cats = {}
            for c, v in column_classification.items():
                if v['keep']:
                    kept_cats.setdefault(v['category'], []).append(c)
            for cat, cols in sorted(kept_cats.items()):
                print(f"    {cat}: {len(cols)} columns")
            print(f"\n  Dropped categories:")
            drop_cats = {}
            for c, v in column_classification.items():
                if not v['keep']:
                    drop_cats.setdefault(v['category'], []).append(c)
            for cat, cols in sorted(drop_cats.items()):
                print(f"    {cat}: {len(cols)} columns")

        # Compute regions
        region = compute_regions(coords, histology, TUMOR_LABELS, STROMA_LABELS, BOUNDARY_ROWS)

        n_spots = len(histology)
        all_abundances_full.append(abund)
        all_abundances_immune.append(abund[immune_cols])
        all_coords.append(coords)
        all_histology.extend(histology)
        all_regions.extend(region)
        all_sample_ids.extend([sid] * n_spots)
        all_mets.extend([row['mets_group']] * n_spots)
        sample_spot_counts[sid] = n_spots

        n_ct = np.sum(region == 'CT')
        n_im = np.sum(region == 'IM')
        n_str = np.sum(region == 'Stroma')
        print(f"  {sid}: {n_spots} spots (CT={n_ct}, IM={n_im}, Stroma={n_str}) — {row['mets_group']}")

    except Exception as e:
        print(f"  {sid}: ERROR — {e}")

# Concatenate
print(f"\nConcatenating {len(sample_spot_counts)} samples...")
X_immune_df = pd.concat(all_abundances_immune, ignore_index=True)
X_full_df = pd.concat(all_abundances_full, ignore_index=True)

# Clean column names for display
immune_col_names = X_immune_df.columns.tolist()
immune_clean_names = [c.replace("meanscell_abundance_w_sf_", "").replace("q05_cell_abundance_w_sf_", "")
                      for c in immune_col_names]
immune_categories = [column_classification[c]['category'] for c in immune_col_names]

X_immune = X_immune_df.values

spot_meta = pd.DataFrame({
    'sample_id': all_sample_ids,
    'mets_group': all_mets,
    'histology': all_histology,
    'region': all_regions,
})

print(f"Immune-focused matrix: {X_immune.shape[0]} spots × {X_immune.shape[1]} cell types")
print(f"Samples loaded: {len(sample_spot_counts)}")
print(f"\nMets distribution:")
print(spot_meta.groupby(['sample_id', 'mets_group']).size().reset_index()
      .groupby('mets_group')['sample_id'].nunique())


# ==============================================================================
# 5. K SWEEP WITH MULTIPLE METRICS
# ==============================================================================

print(f"\n{'=' * 70}")
print("STEP 2: NMF K sweep (reconstruction error + silhouette)")
print("=" * 70)

errors = []
silhouettes = []
stability_scores = []

for k in K_RANGE:
    # Best of N_RESTARTS
    best_err = np.inf
    best_W = None
    all_assignments = []

    for restart in range(N_RESTARTS):
        m = NMF(n_components=k, init='nndsvda',
                random_state=RANDOM_STATE + restart, max_iter=2000, tol=1e-4)
        W_tmp = m.fit_transform(X_immune)
        if m.reconstruction_err_ < best_err:
            best_err = m.reconstruction_err_
            best_W = W_tmp
        all_assignments.append(np.argmax(W_tmp, axis=1))

    errors.append(best_err)

    # Silhouette on a subsample (full dataset too large)
    n_sub = min(10000, X_immune.shape[0])
    rng = np.random.RandomState(RANDOM_STATE)
    idx = rng.choice(X_immune.shape[0], n_sub, replace=False)
    labels_sub = np.argmax(best_W[idx], axis=1)
    if len(np.unique(labels_sub)) > 1:
        sil = silhouette_score(X_immune[idx], labels_sub, sample_size=n_sub, random_state=RANDOM_STATE)
    else:
        sil = -1
    silhouettes.append(sil)

    # Stability: mean pairwise agreement across restarts
    if len(all_assignments) > 1:
        agreements = []
        for a1, a2 in combinations(range(len(all_assignments)), 2):
            agree = np.mean(all_assignments[a1] == all_assignments[a2])
            agreements.append(agree)
        stability_scores.append(np.mean(agreements))
    else:
        stability_scores.append(1.0)

    print(f"  K={k}: error={best_err:.2f}, silhouette={sil:.4f}, stability={stability_scores[-1]:.4f}")

# Plot K selection metrics
fig, axes = plt.subplots(1, 3, figsize=(18, 5), facecolor='white')

axes[0].plot(list(K_RANGE), errors, 'o-', color='#4363d8', linewidth=2, markersize=6)
axes[0].set_xlabel('K', fontsize=13)
axes[0].set_ylabel('Reconstruction error', fontsize=13)
axes[0].set_title('Reconstruction Error', fontsize=14, fontweight='bold')
axes[0].grid(True, alpha=0.3)

axes[1].plot(list(K_RANGE), silhouettes, 'o-', color='#e6194b', linewidth=2, markersize=6)
axes[1].set_xlabel('K', fontsize=13)
axes[1].set_ylabel('Silhouette score', fontsize=13)
axes[1].set_title('Silhouette Score (higher = better)', fontsize=14, fontweight='bold')
axes[1].grid(True, alpha=0.3)

axes[2].plot(list(K_RANGE), stability_scores, 'o-', color='#3cb44b', linewidth=2, markersize=6)
axes[2].set_xlabel('K', fontsize=13)
axes[2].set_ylabel('Pairwise assignment agreement', fontsize=13)
axes[2].set_title('NMF Stability', fontsize=14, fontweight='bold')
axes[2].grid(True, alpha=0.3)

plt.suptitle('Immune-Focused NMF — K Selection Metrics', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'k_selection_metrics.png'), dpi=300, bbox_inches='tight')
plt.show()

# Save metrics for reference
k_metrics = pd.DataFrame({
    'K': list(K_RANGE),
    'reconstruction_error': errors,
    'silhouette': silhouettes,
    'stability': stability_scores,
})
k_metrics.to_csv(os.path.join(OUTPUT_DIR, 'k_selection_metrics.csv'), index=False)
print("\nReview the plots above and set K_FINAL below.")


# ==============================================================================
# 6. SET K_FINAL AND FIT MODEL
#    >>> REVIEW K SELECTION PLOTS ABOVE, THEN SET THIS VALUE <<<
# ==============================================================================

K_FINAL = 8  # ADJUST after reviewing K selection plots

print(f"\n{'=' * 70}")
print(f"STEP 3: Fitting final immune-focused NMF with K={K_FINAL}")
print("=" * 70)

model = NMF(n_components=K_FINAL, init='nndsvda',
            random_state=RANDOM_STATE, max_iter=2000, tol=1e-4)
W = model.fit_transform(X_immune)  # spots × K
H = model.components_               # K × immune_cell_types

print(f"Reconstruction error: {model.reconstruction_err_:.2f}")
print(f"Iterations: {model.n_iter_}")
print(f"W: {W.shape}, H: {H.shape}")

# Add niche info to spot metadata
spot_meta['immune_niche'] = np.argmax(W, axis=1) + 1
for k in range(K_FINAL):
    spot_meta[f'immune_niche_{k+1}_loading'] = W[:, k]


# ==============================================================================
# 7. NICHE CHARACTERIZATION — Name niches by their biology
# ==============================================================================

print(f"\n{'=' * 70}")
print("STEP 4: Characterizing immune niches")
print("=" * 70)

H_norm = H / H.sum(axis=1, keepdims=True)

# Map each column back to its category for aggregated view
cat_contributions = pd.DataFrame(0.0, index=[f"Niche_{i+1}" for i in range(K_FINAL)],
                                  columns=sorted(set(immune_categories)))

for i in range(K_FINAL):
    for j, col in enumerate(immune_col_names):
        cat = column_classification[col]['category']
        cat_contributions.loc[f"Niche_{i+1}", cat] += H_norm[i, j]

print("\n--- Category-level niche composition (normalized) ---")
print(cat_contributions.round(4).to_string())

# Print detailed top contributors per niche
print("\n--- Top 5 specific cell types per niche ---")
niche_labels = {}
for i in range(K_FINAL):
    top_idx = np.argsort(H[i])[::-1][:5]
    top_names = [immune_clean_names[j] for j in top_idx]
    top_vals = [H_norm[i, j] for j in top_idx]
    top_cats = [immune_categories[j] for j in top_idx]

    # Auto-label based on top 2 categories
    top2_cats = []
    seen = set()
    for c in top_cats:
        if c not in seen:
            top2_cats.append(c)
            seen.add(c)
        if len(top2_cats) == 2:
            break
    auto_label = " + ".join(top2_cats)
    niche_labels[i+1] = auto_label

    print(f"\nNiche {i+1} ({auto_label}):")
    for j, idx in enumerate(top_idx):
        print(f"  {immune_clean_names[idx]} [{immune_categories[idx]}]: {top_vals[j]:.4f}")

# Flag CD8-enriched niches
print("\n--- CD8 content per niche ---")
cd8_cols = [j for j, cat in enumerate(immune_categories) if cat in ('CD8_T', 'CD8_T_exhausted')]
for i in range(K_FINAL):
    cd8_total = sum(H_norm[i, j] for j in cd8_cols)
    cd8_cat = {immune_categories[j]: H_norm[i, j] for j in cd8_cols}
    cd8_active = sum(v for k, v in cd8_cat.items() if k == 'CD8_T')
    cd8_exh = sum(v for k, v in cd8_cat.items() if k == 'CD8_T_exhausted')
    flag = " *** CD8-ENRICHED ***" if cd8_total > 0.05 else ""
    print(f"  Niche {i+1} ({niche_labels[i+1]}): CD8_T={cd8_active:.4f}, "
          f"CD8_T_exhausted={cd8_exh:.4f}, total={cd8_total:.4f}{flag}")


# ==============================================================================
# 8. HEATMAP — Cell type contributions per immune niche
# ==============================================================================

print(f"\n{'=' * 70}")
print("STEP 5: Generating plots")
print("=" * 70)

# Detailed heatmap (individual C2L columns)
ct_df = pd.DataFrame(H_norm.T, index=immune_clean_names,
                     columns=[f"Niche {i+1}\n({niche_labels[i+1]})" for i in range(K_FINAL)])

g = sns.clustermap(
    ct_df, cmap='YlOrRd',
    figsize=(max(14, K_FINAL * 1.8), max(10, len(immune_clean_names) * 0.35)),
    row_cluster=True, col_cluster=False,
    linewidths=0.5, linecolor='white',
    dendrogram_ratio=(0.1, 0.0),
    cbar_kws={'label': 'Relative contribution', 'shrink': 0.5},
    yticklabels=True, xticklabels=True,
)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=10, rotation=45, ha='right')
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
g.fig.suptitle(f'Immune-Focused NMF Niche Composition (K={K_FINAL})',
               fontsize=14, fontweight='bold', y=1.02)
plt.savefig(os.path.join(OUTPUT_DIR, 'immune_niche_heatmap_detailed.png'), dpi=300, bbox_inches='tight')
plt.show()

# Category-level heatmap (aggregated)
cat_plot = cat_contributions.T
cat_plot.columns = [f"Niche {i+1}\n({niche_labels[i+1]})" for i in range(K_FINAL)]

fig, ax = plt.subplots(figsize=(max(12, K_FINAL * 1.5), 8), facecolor='white')
sns.heatmap(cat_plot, cmap='YlOrRd', annot=True, fmt='.3f',
            linewidths=0.5, linecolor='white', ax=ax,
            cbar_kws={'label': 'Relative contribution', 'shrink': 0.8})
ax.set_title(f'Immune Niche Composition — Category Level (K={K_FINAL})',
             fontsize=14, fontweight='bold')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'immune_niche_heatmap_categories.png'), dpi=300, bbox_inches='tight')
plt.show()


# ==============================================================================
# 9. NICHE × REGION CROSS-TABULATION
# ==============================================================================

cross_tab = pd.crosstab(spot_meta['immune_niche'], spot_meta['region'], normalize='index') * 100
col_order = [c for c in ['CT', 'IM', 'Stroma', 'Other'] if c in cross_tab.columns]
cross_tab = cross_tab[col_order]

REGION_COLORS = {'CT': '#E64B35', 'IM': '#F39B7F', 'Stroma': '#4DBBD5', 'Other': '#CCCCCC'}

fig, ax = plt.subplots(figsize=(12, 6), facecolor='white')
cross_tab.plot(kind='bar', stacked=True, ax=ax,
               color=[REGION_COLORS[c] for c in col_order],
               edgecolor='white', linewidth=0.5)
ax.set_xlabel('Immune Niche', fontsize=13)
ax.set_ylabel('% of spots', fontsize=13)
ax.set_title(f'Immune Niche × Region (K={K_FINAL})', fontsize=15, fontweight='bold')
ax.legend(title='Region', bbox_to_anchor=(1.02, 1), fontsize=11)
ax.set_xticklabels([f'Niche {i}\n({niche_labels.get(i,"")})' for i in cross_tab.index],
                   rotation=45, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'immune_niche_region_crosstab.png'), dpi=300, bbox_inches='tight')
plt.show()


# ==============================================================================
# 10. NICHE PROPORTIONS BY METS GROUP (SAMPLE-LEVEL) WITH STATS
# ==============================================================================

# Compute sample-level niche proportions
sample_niche = (spot_meta.groupby(['sample_id', 'mets_group', 'immune_niche'])
                .size().reset_index(name='count'))
sample_totals = spot_meta.groupby('sample_id').size().reset_index(name='total')
sample_niche = sample_niche.merge(sample_totals, on='sample_id')
sample_niche['proportion'] = sample_niche['count'] / sample_niche['total'] * 100

mets_order = ['No_Mets', 'LN_Mets', 'Distant_Mets']
mets_palette = {'No_Mets': '#3cb44b', 'LN_Mets': '#ffe119', 'Distant_Mets': '#e6194b'}

# Boxplots
n_cols = min(4, K_FINAL)
n_rows = int(np.ceil(K_FINAL / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 5 * n_rows), facecolor='white')
axes_flat = axes.flatten() if K_FINAL > 1 else [axes]

stats_results = []

for i, ax in enumerate(axes_flat):
    niche_num = i + 1
    if niche_num > K_FINAL:
        ax.axis('off')
        continue

    data = sample_niche[sample_niche['immune_niche'] == niche_num]

    sns.boxplot(data=data, x='mets_group', y='proportion', order=mets_order,
                palette=mets_palette, ax=ax, fliersize=3, width=0.6)
    sns.stripplot(data=data, x='mets_group', y='proportion', order=mets_order,
                  color='black', alpha=0.5, size=4, ax=ax)

    # Kruskal-Wallis test
    groups = [data[data['mets_group'] == mg]['proportion'].values for mg in mets_order]
    groups = [g for g in groups if len(g) >= 2]

    if len(groups) >= 2:
        try:
            kw_stat, kw_p = kruskal(*groups)
        except:
            kw_stat, kw_p = np.nan, 1.0
    else:
        kw_stat, kw_p = np.nan, 1.0

    # Pairwise Mann-Whitney
    pairs_result = {}
    for mg1, mg2 in [('No_Mets', 'LN_Mets'), ('No_Mets', 'Distant_Mets'), ('LN_Mets', 'Distant_Mets')]:
        g1 = data[data['mets_group'] == mg1]['proportion'].values
        g2 = data[data['mets_group'] == mg2]['proportion'].values
        if len(g1) >= 2 and len(g2) >= 2:
            try:
                _, pw_p = mannwhitneyu(g1, g2, alternative='two-sided')
            except:
                pw_p = 1.0
        else:
            pw_p = 1.0
        pairs_result[f'{mg1}_vs_{mg2}'] = pw_p

    stats_results.append({
        'niche': niche_num,
        'label': niche_labels.get(niche_num, ''),
        'kruskal_wallis_p': kw_p,
        **pairs_result,
    })

    # Title with stats
    sig_marker = " *" if kw_p < 0.05 else ""
    label = niche_labels.get(niche_num, '')
    ax.set_title(f'Niche {niche_num} ({label})\nKW p={kw_p:.4f}{sig_marker}',
                 fontsize=10, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('% spots' if i % n_cols == 0 else '')
    ax.tick_params(axis='x', rotation=30, labelsize=8)

plt.suptitle(f'Immune Niche Proportions by Metastasis Group (K={K_FINAL})',
             fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'immune_niche_proportions_by_mets.png'),
            dpi=300, bbox_inches='tight')
plt.show()

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv(os.path.join(OUTPUT_DIR, 'immune_niche_mets_stats.csv'), index=False)
print("\n--- Statistical tests: Niche proportions by mets group ---")
print(stats_df.to_string(index=False))


# ==============================================================================
# 11. NICHE × REGION × METS — What immune niches are in CT vs Stroma by mets?
# ==============================================================================

print(f"\n{'=' * 70}")
print("STEP 6: Immune niche distribution within regions by mets status")
print("=" * 70)

for region_name in ['CT', 'IM', 'Stroma']:
    region_spots = spot_meta[spot_meta['region'] == region_name]
    if len(region_spots) < 100:
        print(f"  {region_name}: too few spots ({len(region_spots)}), skipping")
        continue

    # Sample-level niche proportions within this region
    rn = (region_spots.groupby(['sample_id', 'mets_group', 'immune_niche'])
          .size().reset_index(name='count'))
    rt = region_spots.groupby('sample_id').size().reset_index(name='total')
    rn = rn.merge(rt, on='sample_id')
    rn['proportion'] = rn['count'] / rn['total'] * 100

    # Stacked bar: mean niche proportions per mets group
    pivot = rn.pivot_table(index='mets_group', columns='immune_niche',
                           values='proportion', aggfunc='mean').fillna(0)
    pivot = pivot.reindex(mets_order)

    VIVID_COLORS = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
        '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
        '#aaffc3', '#dcbeff', '#9A6324', '#800000', '#469990'
    ]

    fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')
    pivot.plot(kind='bar', stacked=True, ax=ax,
               color=VIVID_COLORS[:K_FINAL], edgecolor='white', linewidth=0.5)
    ax.set_title(f'{region_name} — Mean Immune Niche Composition by Mets Group',
                 fontsize=14, fontweight='bold')
    ax.set_ylabel('Mean % of spots', fontsize=12)
    ax.set_xlabel('')
    ax.legend(title='Niche', bbox_to_anchor=(1.02, 1), fontsize=9,
              labels=[f'N{i} ({niche_labels.get(i,"")})' for i in pivot.columns])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'immune_niche_{region_name}_by_mets.png'),
                dpi=300, bbox_inches='tight')
    plt.show()

    # Per-niche boxplots within this region
    print(f"\n  --- {region_name}: Niche proportions by mets (Kruskal-Wallis) ---")
    for niche_num in range(1, K_FINAL + 1):
        nd = rn[rn['immune_niche'] == niche_num]
        groups = [nd[nd['mets_group'] == mg]['proportion'].values for mg in mets_order]
        groups = [g for g in groups if len(g) >= 2]
        if len(groups) >= 2:
            try:
                _, p = kruskal(*groups)
                sig = " *" if p < 0.05 else ""
                print(f"    Niche {niche_num} ({niche_labels.get(niche_num,'')}): p={p:.4f}{sig}")
            except:
                print(f"    Niche {niche_num}: test failed")


# ==============================================================================
# 12. PER-SAMPLE SPATIAL MAPS (REPRESENTATIVE EXAMPLES)
# ==============================================================================

example_samples = {
    'No_Mets': '38',
    'LN_Mets': '107',
    'Distant_Mets': '80'
}

cumulative = 0
sample_slices = {}
for sid, n_spots in sample_spot_counts.items():
    sample_slices[sid] = (cumulative, cumulative + n_spots)
    cumulative += n_spots

VIVID_COLORS_NICHE = [
    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
    '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabed4',
    '#aaffc3', '#dcbeff', '#9A6324', '#800000', '#469990'
]

for mets_label, sid in example_samples.items():
    if sid not in sample_slices:
        print(f"  Sample {sid} not found, skipping")
        continue

    s, e = sample_slices[sid]
    dominant_sample = spot_meta.iloc[s:e]['immune_niche'].values
    region_sample = np.array(all_regions[s:e])

    adata_tmp = sc.read_h5ad(os.path.join(SAMPLE_DIR, f"{sid}.h5ad"))
    coords_sample = adata_tmp.obsm['spatial']
    y_sample = coords_sample[:, 1].max() - coords_sample[:, 1]

    fig, axes = plt.subplots(1, 2, figsize=(22, 9), facecolor='white')

    # Panel 1: CT/IM/Stroma
    ax = axes[0]
    for r_label in ['CT', 'IM', 'Stroma', 'Other']:
        mask = region_sample == r_label
        if mask.sum() > 0:
            ax.scatter(coords_sample[mask, 0], y_sample[mask],
                       c=REGION_COLORS[r_label], s=8, edgecolors='none', label=r_label)
    ax.set_title('CT / IM / Stroma', fontsize=16, fontweight='bold')
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='center left', bbox_to_anchor=(0, -0.1), fontsize=9, markerscale=3, ncol=2)

    # Panel 2: Immune niche
    ax = axes[1]
    for n_i in range(K_FINAL):
        mask = dominant_sample == (n_i + 1)
        if mask.sum() > 0:
            label_str = f'Niche {n_i+1} ({niche_labels.get(n_i+1, "")})'
            ax.scatter(coords_sample[mask, 0], y_sample[mask],
                       c=VIVID_COLORS_NICHE[n_i], s=8, edgecolors='none', label=label_str)
    ax.set_title(f'Immune Niche (K={K_FINAL})', fontsize=16, fontweight='bold')
    ax.set_aspect('equal')
    ax.axis('off')
    ax.legend(loc='center left', bbox_to_anchor=(0, -0.15), fontsize=8, markerscale=3, ncol=2)

    plt.suptitle(f'Sample {sid} — {mets_label.replace("_", " ")}',
                 fontsize=20, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f'spatial_immune_niche_{sid}_{mets_label}.png'),
                dpi=300, bbox_inches='tight')
    plt.show()

    # Per-niche loading heatmaps
    W_sample = W[s:e]
    fig, axes_hot = plt.subplots(2, int(np.ceil(K_FINAL/2)),
                                  figsize=(7 * int(np.ceil(K_FINAL/2)), 14), facecolor='black')
    fig.subplots_adjust(hspace=0.3, wspace=0.15)
    for i, ax in enumerate(axes_hot.flat):
        if i >= K_FINAL:
            ax.axis('off')
            continue
        ax.set_facecolor('black')
        vals = W_sample[:, i]
        vmax = np.percentile(vals, 98) if np.percentile(vals, 98) > 0 else 1

        sc_plot = ax.scatter(coords_sample[:, 0], y_sample, c=vals, cmap='hot',
                             s=10, edgecolors='none', vmin=0, vmax=vmax)
        ax.set_title(f'Niche {i+1}\n({niche_labels.get(i+1, "")})',
                     fontsize=11, fontweight='bold', color='white')
        ax.set_aspect('equal')
        ax.axis('off')
        cbar = plt.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
        cbar.ax.yaxis.set_tick_params(color='white')
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')

    plt.suptitle(f'Sample {sid} ({mets_label.replace("_", " ")}) — Immune Niche Loadings',
                 fontsize=16, fontweight='bold', color='white', y=1.02)
    plt.savefig(os.path.join(OUTPUT_DIR, f'immune_niche_loadings_{sid}_{mets_label}.png'),
                dpi=300, bbox_inches='tight', facecolor='black')
    plt.show()


# ==============================================================================
# 13. SAVE RESULTS
# ==============================================================================

print(f"\n{'=' * 70}")
print("STEP 7: Saving results")
print("=" * 70)

# Spot-level metadata with immune niche assignments
spot_meta.to_csv(os.path.join(OUTPUT_DIR, 'spot_immune_niche_assignments.csv'), index=False)
print(f"  Saved spot_immune_niche_assignments.csv ({len(spot_meta)} spots)")

# Niche compositions (H matrix) — both raw and normalized
H_df = pd.DataFrame(H, index=[f"Niche_{i+1}" for i in range(K_FINAL)], columns=immune_clean_names)
H_df.to_csv(os.path.join(OUTPUT_DIR, 'immune_niche_weights.csv'))
print(f"  Saved immune_niche_weights.csv")

H_norm_df = pd.DataFrame(H_norm, index=[f"Niche_{i+1}" for i in range(K_FINAL)], columns=immune_clean_names)
H_norm_df.to_csv(os.path.join(OUTPUT_DIR, 'immune_niche_weights_normalized.csv'))
print(f"  Saved immune_niche_weights_normalized.csv")

# Category-level contributions
cat_contributions.to_csv(os.path.join(OUTPUT_DIR, 'immune_niche_category_contributions.csv'))
print(f"  Saved immune_niche_category_contributions.csv")

# Niche labels
labels_df = pd.DataFrame([{'niche': k, 'label': v} for k, v in niche_labels.items()])
labels_df.to_csv(os.path.join(OUTPUT_DIR, 'immune_niche_labels.csv'), index=False)
print(f"  Saved immune_niche_labels.csv")

# Sample-level niche proportions
sample_niche_prop = (spot_meta.groupby(['sample_id', 'mets_group', 'immune_niche'])
                     .size().unstack(fill_value=0))
sample_niche_prop = sample_niche_prop.div(sample_niche_prop.sum(axis=1), axis=0) * 100
sample_niche_prop.columns = [f'ImmuneNiche_{c}_pct' for c in sample_niche_prop.columns]
sample_niche_prop.to_csv(os.path.join(OUTPUT_DIR, 'sample_immune_niche_proportions.csv'))
print(f"  Saved sample_immune_niche_proportions.csv")

# W matrix (spot loadings) for DE pipeline
np.save(os.path.join(OUTPUT_DIR, 'W_immune_niche_loadings.npy'), W)
print(f"  Saved W_immune_niche_loadings.npy")

# Column mapping for downstream use
col_map = pd.DataFrame({
    'original_column': immune_col_names,
    'clean_name': immune_clean_names,
    'category': immune_categories,
})
col_map.to_csv(os.path.join(OUTPUT_DIR, 'immune_column_mapping.csv'), index=False)
print(f"  Saved immune_column_mapping.csv")

print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("Done!")
print(f"\n>>> Next step: Use spot_immune_niche_assignments.csv to drive pseudobulk DE (Step 4)")
print(f">>> Key file for DE: 'immune_niche' column provides niche grouping per spot")
