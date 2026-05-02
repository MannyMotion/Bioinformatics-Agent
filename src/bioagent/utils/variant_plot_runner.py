"""
variant_plot_runner.py

Purpose: Generate variant annotation Plotly plots in a separate process.
Author:  Emmanuel Ogbu (Manny)
Date:    2026-05-01
"""

import sys
import json
from pathlib import Path
from collections import Counter

import plotly.graph_objects as go


def generate_plots(data: dict) -> list[str]:
    out_dir = Path(data["output_dir"])
    file_name = data["file_name"]
    variants = data["variants"]
    plot_paths = []

    # --- Plot 1: Quality Scores ---
    fig = go.Figure()
    colors = ["#4CAF50" if v["passes_quality"] else "#F44336" for v in variants]
    fig.add_trace(go.Bar(
        x=list(range(len(variants))),
        y=[v["qual"] for v in variants],
        marker_color=colors,
        text=[v["variant_id"][:10] for v in variants],
        textangle=-45
    ))
    fig.add_hline(y=70, line_dash="dash", line_color="white",
                  annotation_text="Quality threshold (70)")
    fig.update_layout(
        title=f"Variant Quality Scores — {file_name}",
        xaxis_title="Variant",
        yaxis_title="Quality Score (PHRED)",
        template="plotly_dark", height=400
    )
    p = str(out_dir / f"{file_name}_quality.html")
    fig.write_html(p)
    plot_paths.append(p)

    # --- Plot 2: Genes ---
    annotated = [v for v in variants if v["is_annotated"]]
    if annotated:
        gene_counts = Counter(v["gene"] for v in annotated)
        fig = go.Figure(go.Bar(
            x=list(gene_counts.values()),
            y=list(gene_counts.keys()),
            orientation="h",
            marker_color="#4f9eff"
        ))
        fig.update_layout(
            title=f"Variants per Gene — {file_name}",
            xaxis_title="Number of Variants",
            yaxis_title="Gene",
            template="plotly_dark", height=400
        )
        p = str(out_dir / f"{file_name}_genes.html")
        fig.write_html(p)
        plot_paths.append(p)

    # --- Plot 3: Clinical Significance Pie ---
    sig_counts = Counter(v["clinical_significance"] for v in variants)
    colors_map = {
        "Pathogenic": "#F44336",
        "Likely pathogenic": "#FF9800",
        "Uncertain significance": "#FFC107",
        "Benign": "#4CAF50",
        "Unknown": "#9E9E9E"
    }
    fig = go.Figure(go.Pie(
        labels=list(sig_counts.keys()),
        values=list(sig_counts.values()),
        marker_colors=[colors_map.get(l, "#9E9E9E") for l in sig_counts.keys()]
    ))
    fig.update_layout(
        title=f"Clinical Significance — {file_name}",
        template="plotly_dark", height=400
    )
    p = str(out_dir / f"{file_name}_significance.html")
    fig.write_html(p)
    plot_paths.append(p)

    return plot_paths


if __name__ == "__main__":
    data_file = Path(sys.argv[1])
    data = json.loads(data_file.read_text())
    paths = generate_plots(data)
    print(json.dumps(paths))