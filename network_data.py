"""Loads ant contact networks into networkx graphs. No plotting, so other scripts can import it freely.

    from network_data import load_graphs
    graphs = load_graphs("colony276")      # {"before": G, "after": G}
    G = graphs["before"]
    G.nodes[n]["group"]                    # "nurse", "untreated" or "treated"
"""
from pathlib import Path
import numpy as np
import networkx as nx

DATA = Path("ant_data")


def file_names(colony):
    return {"before": f"{colony}_pathogen_PreTreatment_TH-27_TD9_interactions_observed.txt",
            "after": f"{colony}_pathogen_PostTreatment_TH21_TD9_interactions_observed.txt"}


def load_treated(colony):
    """Set of tags of the experimentally treated foragers."""
    tw = np.loadtxt(DATA / "treated_worker_list.txt", dtype=str, skiprows=1)  # colony, tag, survived
    return {int(t) for c, t, _ in tw if c == colony.replace("colony", "")}


def build_graph(colony, fname, treated):
    """Nodes = ants (attributes: task, treated, group). Edge = any contact, weight = total contact frames."""
    tg = np.loadtxt(DATA / "task_groups.txt", dtype=str, skiprows=1)  # colony, tag, task_group
    tg = tg[tg[:, 0] == colony]
    d = np.loadtxt(DATA / fname, dtype=int, skiprows=1, usecols=(1, 2, 4))  # Tag1, Tag2, duration

    G = nx.Graph(colony=colony)
    for _, tag, task in tg:
        tag, task = int(tag), str(task)
        if tag in treated:
            group = "treated"
        else:
            group = "untreated" if task == "forager" else "nurse"  # queen counts as nurse
        G.add_node(tag, task=task, treated=tag in treated, group=group)
    for a, b, frames in d:
        a, b, frames = int(a), int(b), int(frames)
        if G.has_edge(a, b):
            G[a][b]["weight"] += frames
        else:
            G.add_edge(a, b, weight=frames)
    return G


def load_graphs(colony="colony276"):
    treated = load_treated(colony)
    return {period: build_graph(colony, fname, treated) for period, fname in file_names(colony).items()}


if __name__ == "__main__":
    # export for other tools (Gephi, R, etc.)
    Path("graphs").mkdir(exist_ok=True)
    for period, G in load_graphs().items():
        nx.write_graphml(G, f"graphs/{G.graph['colony']}_{period}.graphml")
        print(period, G.number_of_nodes(), "nodes", G.number_of_edges(), "edges")