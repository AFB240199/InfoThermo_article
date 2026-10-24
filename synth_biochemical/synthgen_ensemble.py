import numpy as np
import networkx as nx

# -----------------------
# Parameters
# -----------------------
N = 100          # number of nodes
k = 8            # average degree
c = 100          # number of realizations per network
base_seed = 123  # base seed
seeds = base_seed + np.arange(c)


# -----------------------
# Containers (vectors of adjacency matrices)
# -----------------------
A_er_vec  = []
A_ws_vec  = []
A_ba_vec  = []
A_sbm_vec = []
A_rg_vec  = []

# -----------------------
# Fixed parameters (same as notebook)
# -----------------------

# ER
p_er = k / (N - 1)

# WS
k_ws = k
p_ws = 0.1

# BA
m_ba = int(k / 2)

# SBM
m = 4
n_block = N // m
block_sizes = [n_block] * m
p_in = 0.23
p_out = (k / n_block - p_in) / (m - 1)

if p_out < 0 or p_out > 1:
    raise ValueError(f"Invalid p_out={p_out:.4f}. Try another p_in.")

P = [[p_in if i == j else p_out for j in range(m)] for i in range(m)]

# RG
r = np.sqrt(k / (np.pi * N))
r = r + 0.05 * r

# -----------------------
# Network generation loop
# -----------------------
for z in seeds:

    # ER
    G_er = nx.erdos_renyi_graph(N, p_er, seed=int(z))
    A_er_vec.append(nx.to_numpy_array(G_er))

    # WS
    G_ws = nx.watts_strogatz_graph(N, k_ws, p_ws, seed=int(z))
    A_ws_vec.append(nx.to_numpy_array(G_ws))

    # BA
    G_ba = nx.barabasi_albert_graph(N, m_ba, seed=int(z))
    A_ba_vec.append(nx.to_numpy_array(G_ba))

    # SBM
    G_sbm = nx.stochastic_block_model(block_sizes, P, seed=int(z))
    A_sbm_vec.append(nx.to_numpy_array(G_sbm))

    # RG
    G_rg = nx.random_geometric_graph(N, r, seed=int(z))
    A_rg_vec.append(nx.to_numpy_array(G_rg))

# -----------------------
# Convert to numpy object arrays (vector of matrices)
# -----------------------
A_er_vec  = np.array(A_er_vec,  dtype=object)
A_ws_vec  = np.array(A_ws_vec,  dtype=object)
A_ba_vec  = np.array(A_ba_vec,  dtype=object)
A_sbm_vec = np.array(A_sbm_vec, dtype=object)
A_rg_vec  = np.array(A_rg_vec,  dtype=object)

# -----------------------
# Save to disk
# -----------------------
np.savez(
    "network_ensembles.npz",
    A_er=A_er_vec,
    A_ws=A_ws_vec,
    A_ba=A_ba_vec,
    A_sbm=A_sbm_vec,
    A_rg=A_rg_vec,
    seeds=seeds,
)

print("Saved network ensembles:")
print(" - ER, WS, BA, SBM, RG")
print(f" - {c} realizations each")
