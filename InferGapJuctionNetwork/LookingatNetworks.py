# %%
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
G = nx.read_gml('Islet_Networks/size4171_it-270_er7.29.gml')
#G = nx.read_gml('401_0.27.gml')

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

edge_weights = nx.get_edge_attributes(G, 'weight')
# %%

Z = [x for _,x in sorted(zip(deg,av_conduct))]


# %%
plt.figure
options = {"node_size": 30, "alpha": 0.9, "edge_color": [0.9, 0.9,0.9]}
nx.draw(G, **options)
# %%
#calculate percent frequency
gjdist = pd.read_csv('TotGJConductDistribution.csv', index_col = 0)
gjconduct = [float(item)*203 for item in list(gjdist)]
gjfreq = list(gjdist.iloc[0,:])
gjfreq_s = list(gjdist.iloc[2,:])

conduct = nx.get_node_attributes(G,'Conduct')
x = list(conduct.values())
x = [float(i) for i in x]
hist = np.zeros(np.shape(gjfreq))
gjconduct[-1]=1000
for i in range(0, len(gjconduct)):
    if i == 0:
        hist[i] = len([k for k in x if k <= gjconduct[i]])/len(x)
    else:
        hist[i] = len([k for k in x if k <= gjconduct[i] and k > gjconduct[i-1]])/len(x)
    tot_err = (sum(((hist - gjfreq)/gjfreq_s)**2))**(1/2)
gjconduct[-1]=2.639



# %%
plt.figure
plt.bar(gjconduct, gjfreq)
plt.bar(gjconduct, hist)

# %%
# %%
plt.figure
plt.bar(gjconduct, hist)

plt.figure
plt.bar(gjconduct, gjfreq)

# %%
