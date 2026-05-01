"""
plot_runner.py

Purpose: Generate Plotly HTML plots in a separate process.
         Uses Plotly instead of matplotlib to avoid Windows
         segfault issues with matplotlib in subprocesses.

Author:  Emmanuel Ogbu (Manny)
Date:    2026-05-01
"""

import sys
import json
import os
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def generate_plots(data: dict) -> list[str]:
    """Generate all QC plots as HTML files."""
    out_dir = Path(data["output_dir"])
    file_name = data["file_name"]
    plot_paths = []

    gc_values = data["gc_values"]
    lengths = data["lengths"]
    mean_gc = data["mean_gc"]
    median_length = data["median_length"]
    nuc_counts = data["nucleotide_counts"]

    # Plot 1: GC Content Distribution
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=gc_values, nbinsx=20,
                               marker_color="#2196F3", opacity=0.8,
                               name="GC Content"))
    fig.add_vline(x=mean_gc, line_dash="dash", line_color="#F44336",
                  annotation_text=f"Mean: {mean_gc}%")
    fig.add_vline(x=35.0, line_dash="dot", line_color="#FF9800",
                  annotation_text="Low (35%)")
    fig.add_vline(x=65.0, line_dash="dot", line_color="#FF9800",
                  annotation_text="High (65%)")
    fig.update_layout(
        title=f"GC Content Distribution — {file_name}",
        xaxis_title="GC Content (%)",
        yaxis_title="Number of Sequences",
        template="plotly_dark",
        height=400
    )
    p = str(out_dir / f"{file_name}_gc_content.html")
    fig.write_html(p)
    plot_paths.append(p)

    # Plot 2: Sequence Length Distribution
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=lengths, nbinsx=20,
                               marker_color="#4CAF50", opacity=0.8,
                               name="Length"))
    fig.add_vline(x=median_length, line_dash="dash", line_color="#F44336",
                  annotation_text=f"Median: {median_length:.0f}bp")
    fig.update_layout(
        title=f"Sequence Length Distribution — {file_name}",
        xaxis_title="Sequence Length (bp)",
        yaxis_title="Number of Sequences",
        template="plotly_dark",
        height=400
    )
    p = str(out_dir / f"{file_name}_lengths.html")
    fig.write_html(p)
    plot_paths.append(p)

    # Plot 3: Nucleotide Composition
    total = sum(nuc_counts.values())
    nucs = list(nuc_counts.keys())
    pcts = [v/total*100 for v in nuc_counts.values()] if total > 0 else [0]*5
    colors = ["#4CAF50", "#F44336", "#FF9800", "#2196F3", "#9E9E9E"]

    fig = go.Figure()
    fig.add_trace(go.Bar(x=nucs, y=pcts, marker_color=colors,
                         text=[f"{p:.1f}%" for p in pcts],
                         textposition="outside"))
    fig.update_layout(
        title=f"Nucleotide Composition — {file_name}",
        xaxis_title="Nucleotide",
        yaxis_title="Percentage (%)",
        template="plotly_dark",
        height=400
    )
    p = str(out_dir / f"{file_name}_composition.html")
    fig.write_html(p)
    plot_paths.append(p)

    return plot_paths


if __name__ == "__main__":
    data_file = Path(sys.argv[1])
    data = json.loads(data_file.read_text())
    paths = generate_plots(data)
    print(json.dumps(paths))