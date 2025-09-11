"""
Network analysis and visualization.

Build simple bipartite graphs from results and provide a quick drawing helper.
"""

from __future__ import annotations

from typing import Tuple

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


def build_compound_target_graph(df: pd.DataFrame) -> nx.Graph:
    G = nx.Graph()
    if not {"compound_id", "target_uniprot"}.issubset(df.columns):
        return G
    for _, row in df.iterrows():
        c = f"C:{row['compound_id']}"
        t = f"T:{row['target_uniprot']}"
        G.add_node(c, bipartite=0)
        G.add_node(t, bipartite=1)
        G.add_edge(c, t, weight=float(row.get("fused_score", 0.0)))
    return G


def layout_bipartite(G: nx.Graph) -> dict:
    top = {n for n, d in G.nodes(data=True) if d.get("bipartite") == 0}
    bottom = set(G) - top
    return nx.bipartite_layout(G, top_nodes=top)


def draw_bipartite(G: nx.Graph):
    """Draw a bipartite graph using networkx + matplotlib and return the figure."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if G.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "No nodes to visualize", ha="center", va="center")
        ax.set_axis_off()
        return fig
    pos = layout_bipartite(G)
    top_nodes = [n for n, d in G.nodes(data=True) if d.get("bipartite") == 0]
    bottom_nodes = [n for n in G if n not in top_nodes]
    nx.draw_networkx_nodes(G, pos, nodelist=top_nodes, node_color="#4e79a7", node_size=120, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=bottom_nodes, node_color="#f28e2b", node_size=120, ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.3, ax=ax)
    # Keep labels minimal to avoid clutter
    if G.number_of_nodes() <= 40:
        labels = {n: n.split(":", 1)[1] for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)
    ax.set_title("Compound–Target Bipartite Network")
    ax.axis("off")
    return fig
