import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from network_data import load_graphs

COLORS = {"nurse": "seagreen", "untreated": "gold", "treated": "crimson"}
CENTERS = {"nurse": (-2.5, 0), "untreated": (2.0, 1.8), "treated": (2.0, -1.8)}  # one cluster per group


def grouped_layout(G):
    """One spring-layout cluster per group, placed at a fixed center."""
    pos = {}
    for g, center in CENTERS.items():
        members = [n for n, grp in G.nodes(data="group") if grp == g]
        radius = 0.18 * np.sqrt(len(members))  # bigger groups get a bigger cluster
        pos.update(nx.spring_layout(G.subgraph(members), weight="weight", seed=1, center=center, scale=radius))
    return pos


graphs = load_graphs("colony276")
pos = grouped_layout(graphs["before"])  # same positions in both panels, so only the edges change

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
for ax, (period, G) in zip(axes, graphs.items()):
    nx.draw_networkx_edges(G, pos, edge_color="gray", width=0.3, alpha=0.3, ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=[COLORS[g] for _, g in G.nodes(data="group")], node_size=70,
                           edgecolors="black", linewidths=0.4, ax=ax)
    ax.set_title(period.capitalize() + " treatment")
    ax.axis("off")

legend = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=10, label=l)
          for l, c in [("Nurses (incl. queen)", COLORS["nurse"]), ("Untreated foragers", COLORS["untreated"]),
                       ("Treated foragers", COLORS["treated"])]]
fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False)
plt.tight_layout(rect=(0, 0.05, 1, 1))
plt.savefig("containment.png", dpi=150)
plt.show()