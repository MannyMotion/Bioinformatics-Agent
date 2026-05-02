"""
rnaseq_plot_runner.py

Purpose: Generate RNA-seq Plotly plots in a separate process.
         PCA and heatmap computed here to avoid scikit-learn
         crashing the uvicorn server process on Windows.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-05-01
"""

import sys
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
import plotly.graph_objects as go


def generate_plots(data: dict) -> list[str]:
    out_dir = Path(data["output_dir"])
    file_name = data["file_name"]
    plot_paths = []

    gene_results = data["gene_results"]
    control_label = data["control_label"]
    treatment_label = data["treatment_label"]
    control_cols = data["control_cols"]
    treatment_cols = data["treatment_cols"]
    cpm_data = data["cpm_data"]

    # --- Plot 1: Volcano Plot ---
    up_x, up_y, up_names = [], [], []
    down_x, down_y, down_names = [], [], []
    ns_x, ns_y = [], []

    for g in gene_results:
        neg_log10_p = -np.log10(g["p_value"] + 1e-10)
        if g["significant"] and g["log2_fold_change"] > 0:
            up_x.append(g["log2_fold_change"])
            up_y.append(neg_log10_p)
            up_names.append(g["gene_id"])
        elif g["significant"] and g["log2_fold_change"] < 0:
            down_x.append(g["log2_fold_change"])
            down_y.append(neg_log10_p)
            down_names.append(g["gene_id"])
        else:
            ns_x.append(g["log2_fold_change"])
            ns_y.append(neg_log10_p)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=ns_x, y=ns_y, mode="markers",
                             marker=dict(color="#9E9E9E", size=8, opacity=0.6),
                             name="Not significant"))
    fig.add_trace(go.Scatter(x=down_x, y=down_y, mode="markers+text",
                             marker=dict(color="#2196F3", size=10),
                             text=down_names, textposition="top center",
                             name=f"Down ({len(down_x)})"))
    fig.add_trace(go.Scatter(x=up_x, y=up_y, mode="markers+text",
                             marker=dict(color="#F44336", size=10),
                             text=up_names, textposition="top center",
                             name=f"Up ({len(up_x)})"))
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="white", opacity=0.5)
    fig.add_vline(x=1.0, line_dash="dot", line_color="white", opacity=0.5)
    fig.add_vline(x=-1.0, line_dash="dot", line_color="white", opacity=0.5)
    fig.update_layout(
        title=f"Volcano Plot: {treatment_label} vs {control_label}",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(p-value)",
        template="plotly_dark", height=500
    )
    p = str(out_dir / f"{file_name}_volcano.html")
    fig.write_html(p)
    plot_paths.append(p)

    # --- Plot 2: PCA ---
    try:
        all_cols = control_cols + treatment_cols
        genes = list(cpm_data[all_cols[0]].keys()) if all_cols else []

        # Build matrix: rows=samples, cols=genes
        matrix = []
        for col in all_cols:
            row = [cpm_data[col].get(g, 0.0) for g in genes]
            matrix.append(row)

        matrix = np.array(matrix, dtype=float)

        # Simple PCA without scikit-learn
        matrix_centered = matrix - matrix.mean(axis=0)
        cov = np.cov(matrix_centered.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvectors = eigenvectors[:, idx]
        eigenvalues = eigenvalues[idx]
        coords = matrix_centered @ eigenvectors[:, :2]

        total_var = eigenvalues.sum()
        var_explained = eigenvalues[:2] / total_var * 100 if total_var > 0 else [0, 0]

        fig = go.Figure()
        for i, col in enumerate(control_cols):
            fig.add_trace(go.Scatter(
                x=[coords[i, 0]], y=[coords[i, 1]],
                mode="markers+text",
                marker=dict(color="#4CAF50", size=14),
                text=[col], textposition="top center",
                name=control_label, legendgroup="control",
                showlegend=(i == 0)
            ))

        offset = len(control_cols)
        for i, col in enumerate(treatment_cols):
            fig.add_trace(go.Scatter(
                x=[coords[offset + i, 0]], y=[coords[offset + i, 1]],
                mode="markers+text",
                marker=dict(color="#F44336", size=14),
                text=[col], textposition="top center",
                name=treatment_label, legendgroup="treatment",
                showlegend=(i == 0)
            ))

        fig.update_layout(
            title="PCA Plot — Sample Clustering",
            xaxis_title=f"PC1 ({var_explained[0]:.1f}% variance)",
            yaxis_title=f"PC2 ({var_explained[1]:.1f}% variance)",
            template="plotly_dark", height=450
        )
        p = str(out_dir / f"{file_name}_pca.html")
        fig.write_html(p)
        plot_paths.append(p)

    except Exception as e:
        print(f"PCA failed: {e}", file=sys.stderr)

    # --- Plot 3: Heatmap ---
    try:
        sig_genes = [g for g in gene_results if g["significant"]]
        top_genes = sorted(sig_genes, key=lambda x: abs(x["log2_fold_change"]), reverse=True)[:10]
        top_ids = [g["gene_id"] for g in top_genes]
        all_cols = control_cols + treatment_cols

        if top_ids and all_cols:
            z_vals = []
            for gene_id in top_ids:
                row = [cpm_data[col].get(gene_id, 0.0) for col in all_cols]
                z_vals.append(row)

            fig = go.Figure(data=go.Heatmap(
                z=z_vals,
                x=all_cols,
                y=top_ids,
                colorscale="RdBu_r",
                colorbar=dict(title="CPM"),
                text=[[f"{v:.0f}" for v in row] for row in z_vals],
                texttemplate="%{text}"
            ))
            fig.update_layout(
                title="Top Differentially Expressed Genes (CPM)",
                template="plotly_dark", height=450
            )
            p = str(out_dir / f"{file_name}_heatmap.html")
            fig.write_html(p)
            plot_paths.append(p)

    except Exception as e:
        print(f"Heatmap failed: {e}", file=sys.stderr)

    return plot_paths


if __name__ == "__main__":
    data_file = Path(sys.argv[1])
    data = json.loads(data_file.read_text())
    paths = generate_plots(data)
    print(json.dumps(paths))