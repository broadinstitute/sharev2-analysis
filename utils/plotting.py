import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from enum import Enum


def plot_knee_curve(adata, output_folder="./", prefix="sample"):
    """
    Plot knee curve separating lines by lab sample ID
    """
    output_dir = Path(output_folder) / prefix
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    fontsize = 18
    colors = plt.cm.tab10.colors

    for i, sample in enumerate(adata.obs["lab_sample_id"].unique()):
        adata_subset = adata[adata.obs["lab_sample_id"] == sample]
        umi_counts = np.array(adata_subset.X.sum(axis=1)).flatten()
        umi_counts.sort()
        umi_counts = umi_counts[::-1]

        ax.loglog(range(len(umi_counts)), umi_counts, linewidth=2, label=sample, color=colors[i % len(colors)])

    ax.set_ylabel("UMI Counts", fontsize=fontsize)
    ax.set_xlabel("Barcodes", fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.axhline(y=500, linewidth=2, color="#505050", linestyle='--', label='500 UMI')
    ax.axhline(y=200, linewidth=2, color="#A0A0A0", linestyle='--', label='200 UMI')
    ax.legend(fontsize=14)
    plt.grid(True, which="both")

    plt.savefig(os.path.join(output_dir, "sample_knee_plot.png"), format='png', bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(output_dir, "sample_knee_plot.svg"), format="svg", bbox_inches="tight")
    plt.show()


def plot_expression_heatmap(df, output_folder="./", prefix="sample", gene=None, title=""):
    """
    Plot pseudobulk expression heatmap for a gene across a 96-well plate
    """
    output_dir = Path(output_folder) / prefix
    os.makedirs(output_dir, exist_ok=True)

    if gene not in df.index:
        raise ValueError(f"'{gene}' not found in the pseudobulk DataFrame.")

    gene_data = df.loc[gene]
    rows = list('ABCDEFGH')
    cols = list(range(1, 13))
    plate_df = pd.DataFrame(0, index=rows, columns=cols)

    for well, expression in gene_data.items():
        row = well[0]
        col = int(well[1:])
        plate_df.at[row, col] = expression

    plt.figure(figsize=(8, 4))
    sns.heatmap(
        plate_df, annot=True, fmt='.1f', cmap='viridis', cbar=True, linewidths=0.5,
        annot_kws={"size": 10}, cbar_kws={'label': 'Gene Expression'}
    )

    plot_title = f'{gene} {title}' if title else f'{gene} pseudobulk counts'
    plt.title(plot_title, fontsize=16)
    plt.xlabel('Column', fontsize=14)
    plt.ylabel('Row', fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(rotation=0, fontsize=12)

    plt.savefig(os.path.join(output_dir, f"{gene}_plate_heatmap.png"), format='png', bbox_inches='tight', dpi=300)
    plt.show()


def stacked_barplot_proportions(
    adata, cluster_key, var_key, fsize=(5, 12), annotations=True,
    reverse_order=False, color_dict=None, cluster_order=None,
    output_folder="./", prefix="sample"
):
    """
    Creates a horizontal stacked bar plot of proportions for a categorical variable across clusters.
    """
    output_dir = Path(output_folder) / prefix
    os.makedirs(output_dir, exist_ok=True)

    grouped_data = adata.groupby([cluster_key, var_key]).size().unstack(fill_value=0)
    proportions = grouped_data.div(grouped_data.sum(axis=1), axis=0)

    if cluster_order:
        proportions = proportions.loc[cluster_order]
    elif reverse_order:
        proportions = proportions.iloc[::-1]

    cluster_sizes = grouped_data.sum(axis=1)
    cluster_sizes = cluster_sizes.loc[proportions.index]

    unique_categories = grouped_data.columns.tolist()
    if color_dict:
        colors = [color_dict.get(cat, "gray") for cat in unique_categories]
    else:
        colors = sns.color_palette("husl", n_colors=len(unique_categories))

    _, ax = plt.subplots(figsize=fsize)
    proportions.plot(kind="barh", stacked=True, color=colors, ax=ax, width=0.8, edgecolor=None)

    if annotations:
        for i, txt in enumerate(cluster_sizes):
            ax.text(1.02, i, f"{txt:,}", fontsize=12, va="center", transform=ax.get_yaxis_transform())

    ax.set_xlim(0, 1.1)
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xlabel("Proportion", fontsize=14)
    ax.set_ylabel(cluster_key, fontsize=14)
    ax.set_title(f'{var_key} by {cluster_key}', fontsize=16)

    if annotations:
        ax.legend(title=var_key, bbox_to_anchor=(1.2, 1), loc="upper left")
    else:
        ax.get_legend().remove()

    plt.grid(False)
    plt.savefig(os.path.join(output_dir, f"{var_key}_by_{cluster_key}_barplot.png"), format='png', bbox_inches='tight', dpi=300)
    plt.show()


