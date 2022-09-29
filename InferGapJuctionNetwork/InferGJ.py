## InferGJ.py
# This code uses FRAP data from Farnsworth et al. 2014 to infer many different configurartions of GJ weights and networks in the islet
# %%
from asyncio import base_tasks
from itertools import product
import math
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy.random as npr
from scipy import optimize as op
import scipy.stats as stats 
import random
import seaborn as sns

# %% Define skew
def skew_norm_pdf(x,e=0,w=1,a=0):
    # adapated from:
    # http://stackoverflow.com/questions/5884768/skew-normal-distribution-in-scipy
    t = (x-e) / w
    return 2.0 * w * stats.norm.pdf(t) * stats.norm.cdf(a*t)

# %% Define optimization function
def calc_conduct_optim(edgeweights, G, gjconduct, gjfreq):
    for i, x in enumerate(G.edges):
        G.edges[x[0],x[1]]['weight']=edgeweights[i]

    for node in G.nodes:
        tot_cond = 0
        for i in list(G.neighbors(node)):
            tot_cond = tot_cond+G.get_edge_data(node,i)['weight']
        G.nodes[node]['Conduct'] = tot_cond

    conduct = nx.get_node_attributes(G,'Conduct')
    x = list(conduct.values())
    tot_error = cost(x, gjconduct, gjfreq)
    #print(tot_error)
    return tot_error

def calc_conduct(edgeweights, G, gjconduct, gjfreq):
    for i, x in enumerate(G.edges):
        G.edges[x[0],x[1]]['weight']=str(edgeweights[i])

    for node in G.nodes:
        tot_cond = 0
        for i in list(G.neighbors(node)):
            tot_cond = tot_cond+float(G.get_edge_data(node,i)['weight'])
        G.nodes[node]['Conduct'] = str(tot_cond)
    conduct = nx.get_node_attributes(G,'Conduct')
    x = [float(i) for i in list(conduct.values())]
    return G


def cost(x,gjconduct,gjfreq):
    #calculate percent frequency
    hist = np.zeros(np.shape(gjfreq))
    gjconduct.append(1000)
    for i in range(0, len(gjconduct)-1):
        if i == 0:
            hist[i] = len([k for k in x if k <= gjconduct[i]])/len(x)
        else:
            hist[i] = len([k for k in x if k <= gjconduct[i] and k > gjconduct[i-1]])/len(x)
    gjconduct.pop()
    tot_err = np.linalg.norm(np.subtract(hist, gjfreq))
    return tot_err

# %%
def run_networkbuild(itter):
    random.seed(itter)
    # First build a 3D network with sphere packing:
    b_diameter = 10 #diameter of beta cell: mu m
    i_area = 10000 #diameter of islet: mu m^2
    perc_b = .80 #percent of beta cells in the islet
    i_radius = (i_area/math.pi)**(1/2)

    beta_num = int(perc_b*i_area/b_diameter)
    #% Add edges
    #1. Randomly pick edge distrutions
    edgedist = np.random.normal(6.38+0.6, 1.35, beta_num) 

    #% Pack 3D grid with cells
    grid_rad = int(i_radius/b_diameter)
    points_inside_grid = product(range(-grid_rad,grid_rad),range(-grid_rad,grid_rad),range(-grid_rad,grid_rad))
    points_inside_grid = list(points_inside_grid)
    positions = [columns for columns in points_inside_grid if (columns[1]**2+columns[2]**2+columns[0]**2)**(1/2)<i_radius]
    beta_positions = np.array(positions)[np.random.choice(len(positions), size=int(perc_b*len(positions)), replace=False)]

    # Fill a network of radius i_radius with betacells
    G = nx.Graph()
    G.add_nodes_from(range(0,len(beta_positions)))
    positions_dict = {i:list(beta_positions[i]) for i in G.nodes}
    for n in G.nodes:
        G.nodes[n]["pos"] = positions_dict[n]
    thr = 6 #How far away the edge can be.
    node_sort = np.argsort(-edgedist)
    used_nodes=[]
    for i in node_sort: 
        #check if requirement is already satisfied
        if int(edgedist[i]) - G.degree[i] > 0:
            nodepos = G.nodes[i]['pos']
            avaliable_edges = [k for (k,v) in positions_dict.items() if 
                (v[0]-nodepos[0])**2+(v[1]-nodepos[1])**2+(v[2]-nodepos[2])**2**(1/2)<thr 
                and k not in used_nodes and k != i and G.degree[k]<edgedist[k]-1]
            #choose which nodes to connect proximal edges based on degree
            
            if len(avaliable_edges) > int(edgedist[i])-G.degree[i]:
                edges_to_connect = np.random.choice(avaliable_edges, int(edgedist[i])-G.degree[i])
                [G.add_edge(i, ed) for ed in edges_to_connect] #add edges
                used_nodes.append(i)
            else:
                avaliable_edges = [k for (k,v) in positions_dict.items() if 
                    (v[0]-nodepos[0])**2+(v[1]-nodepos[1])**2+(v[2]-nodepos[2])**2**(1/2)<thr]
                edges_to_connect = np.random.choice(avaliable_edges, int(edgedist[i])-G.degree[i])
                [G.add_edge(i, ed) for ed in edges_to_connect] #add edges
        else:
            used_nodes.append(i)
    G.remove_edges_from(nx.selfloop_edges(G))
    degree_sequence = [d for n, d in G.degree()]
    print(np.mean(degree_sequence))
    print(np.std(degree_sequence))

    #Assign GJ weighting such that the histogram looks like GJdist
    #
    gjdist = pd.read_csv('TotGJConductDistribution.csv')
    gj_conduct = [float(item)*203 for item in list(gjdist)]
    gj_freq = list(gjdist.iloc[0,:])
    #
    x=np.linspace(0,1.22,G.number_of_edges())
    edgeweights_0 = skew_norm_pdf(x,1.22/8,0.5,-0.1)
    #random.shuffle(edgeweights_0)
    #edgeweights_0 = npr.poisson(122/4, size=G.number_of_edges())/100
    edgeweights_dist = npr.weibull(5, size=G.number_of_edges())/100
    edgeweights_0 = np.interp(edgeweights_0, (edgeweights_0.min(),edgeweights_0.max()), (0, 1.22/3.5))
    #plt.hist(edgeweights_0)
    #print(edgeweights_0.mean())
    edgeweights_0 = [1e-5 if x<=0 else x for x in edgeweights_0] 
    final_weights = op.fmin(calc_conduct_optim, edgeweights_0, args = (G, gj_conduct, gj_freq),ftol=1e-15,xtol=1e-15,maxiter=1e10)

    
    err = calc_conduct_optim(final_weights, G, gj_conduct, gj_freq)
    print(err)
    G2 = calc_conduct(final_weights, G, gj_conduct, gj_freq)
    
    for i in range(0, len(G2.nodes)):
        G2.nodes[i]['pos'] = str(G2.nodes[i]['pos'])
        G2.nodes[i]['Conduct'] = str(G2.nodes[i]['Conduct'])
    
    if err < 0.16:
        nx.write_gml(G2, str(itter) + '_' + str(round(err,2))+'.gml')

if __name__ == "__main__":
    for i in range(1,150):
        run_networkbuild(i)

# %%
