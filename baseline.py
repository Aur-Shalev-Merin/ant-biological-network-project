from pathlib import Path
import numpy as np
import networkx as nx

DATA = Path("ant_data")
colony = "colony276"
colony_num = colony.replace("colony", "")

FILES = [("before", f"{colony}_pathogen_PreTreatment_TH-27_TD9_interactions_observed.txt"),
         ("after", f"{colony}_pathogen_PostTreatment_TH21_TD9_interactions_observed.txt")]


def build_graph(fname):
    """Contact network: nodes = ants, edge = at least one contact, weight = total contact frames."""
    tg = np.loadtxt(DATA / "task_groups.txt", dtype=str, skiprows=1)  # colony, tag, task_group
    tg = tg[tg[:, 0] == colony]
    d = np.loadtxt(DATA / fname, dtype=int, skiprows=1, usecols=(1, 2, 4))  # Tag1, Tag2, duration

    G = nx.Graph()
    for _, tag, group in tg:
        G.add_node(int(tag), task=group)
    for a, b, frames in d:
        if G.has_edge(a, b):
            G[a][b]["weight"] += frames
        else:
            G.add_edge(a, b, weight=frames)
    return G


def mean_degrees(G, infected):
    """Return Kav (all ants) and Kinfected_av (infected ants only)."""
    Kav = np.mean([k for _, k in G.degree()])
    Kinf = np.mean([G.degree(n) for n in infected])
    return Kav, Kinf


def group_transitivity(G, infected):
    """Transitivity of the whole network and of the subgraphs induced by infected / uninfected ants."""
    uninfected = [n for n in G if n not in infected]
    return (nx.transitivity(G),
            nx.transitivity(G.subgraph(infected)),
            nx.transitivity(G.subgraph(uninfected)))


def block_matrix(G):
    """2x2 edge probabilities: [[N-N, N-F], [N-F, F-F]] using known task labels (queen counted as nurse)."""
    F = [n for n, t in G.nodes(data="task") if t == "forager"]
    N = [n for n, t in G.nodes(data="task") if t != "forager"]
    nn = G.subgraph(N).number_of_edges() / (len(N) * (len(N) - 1) / 2)
    ff = G.subgraph(F).number_of_edges() / (len(F) * (len(F) - 1) / 2)
    nf = len(list(nx.edge_boundary(G, N, F))) / (len(N) * len(F))
    return np.array([[nn, nf], [nf, ff]])


# infected = experimentally exposed foragers
tw = np.loadtxt(DATA / "treated_worker_list.txt", dtype=str, skiprows=1)  # colony, tag, survived
infected = [int(t) for c, t, _ in tw if c == colony_num]

for period, fname in FILES:
    G = build_graph(fname)
    Kav, Kinf = mean_degrees(G, infected)
    total, inf, uninf = group_transitivity(G, infected)
    print(f"{period}: Kav = {Kav:.1f}   Kinfected_av = {Kinf:.1f}")
    print(f"{period}: transitivity  total = {total:.3f}   infected = {inf:.3f}   uninfected = {uninf:.3f}")
    print(f"{period}: block matrix (rows/cols: nurses, foragers)")
    print(np.round(block_matrix(G), 3))

#run 3x3 sbm:
def block_matrix_3(G, infected):
    """3x3 edge probabilities between: nurses (incl. queen), untreated foragers, treated foragers."""
    inf = set(infected)
    nurses = [n for n, t in G.nodes(data="task") if t != "forager"]
    ufor = [n for n, t in G.nodes(data="task") if t == "forager" and n not in inf]
    ifor = [n for n in G if n in inf]
    groups = [nurses, ufor, ifor]

    M = np.zeros((3, 3))
    for i, A in enumerate(groups):
        for j, B in enumerate(groups):
            if i == j:
                possible = len(A) * (len(A) - 1) / 2
                M[i, j] = G.subgraph(A).number_of_edges() / possible if possible else np.nan
            else:
                M[i, j] = len(list(nx.edge_boundary(G, A, B))) / (len(A) * len(B))
    return M, [len(g) for g in groups]

for period, fname in FILES:
    G = build_graph(fname)
    M, sizes = block_matrix_3(G, infected)
    print(period, "sizes (nurses, untreated foragers, treated foragers):", sizes)
    print(np.round(M, 3))