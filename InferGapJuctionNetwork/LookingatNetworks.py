# %%
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
G = nx.read_gml('0_0.15.gml')
G = nx.read_gml('401_0.27.gml')

# %%
av_conduct = []
tot_conduct = []
deg = []
for node in G.nodes():
    weights = []
    for edge in G.edges(node):
        weights.append(float(G[edge[0]][edge[1]]['weight']))
    av_conduct.append(np.mean(weights))
    tot_conduct.append(sum(weights))
    deg.append(len(weights))
# %%

Z = [x for _,x in sorted(zip(deg,av_conduct))]


# %%
