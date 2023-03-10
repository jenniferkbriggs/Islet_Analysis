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

# %%
class MyBounds:
    def __init__(self, xmax=[1.1,1.1], xmin=[-1.1,-1.1] ):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)
    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


# %% Define skew
def skew_norm_pdf(x,e=0,w=1,a=0):
    # adapated from:
    # http://stackoverflow.com/questions/5884768/skew-normal-distribution-in-scipy
    t = (x-e) / w
    return 2.0 * w * stats.norm.pdf(t) * stats.norm.cdf(a*t)

# %% Define optimization function
def calc_conduct_optim(edgeweights, G, gjconduct, gjfreq, gj_freq_s, options=0):
    for i, x in enumerate(G.edges):
        G.edges[x[0],x[1]]['weight']=edgeweights[i]

    for node in G.nodes:
        tot_cond = 0
        for i in list(G.neighbors(node)):
            tot_cond = tot_cond+G.get_edge_data(node,i)['weight']
        G.nodes[node]['Conduct'] = tot_cond

    conduct = nx.get_node_attributes(G,'Conduct')
    x = list(conduct.values())
    tot_error, Hist = cost(x, gjconduct, gjfreq, gj_freq_s)
    #print(tot_error)
    if options == 1:
        return tot_error, Hist 
    else:
        return tot_error

def calc_conduct(edgeweights, G, gjconduct):
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


def cost(x,gjconduct,gjfreq, gj_freq_s):
    #calculate percent frequency
    hist = np.zeros(np.shape(gjfreq))
    gjconduct[-1]=1000
    for i in range(0, len(gjconduct)-1):
        if i == 0:
            hist[i] = len([k for k in x if k <= gjconduct[i]])/len(x)
        else:
            hist[i] = len([k for k in x if k <= gjconduct[i] and k > gjconduct[i-1]])/len(x)
    tot_err = (sum(((hist - gjfreq)/gj_freq_s)**2))**(1/2)
    #print(tot_err)
    return tot_err, hist

# %%
def run_networkbuild(itter, area):
    random.seed(itter)
    # First build a 3D network with sphere packing:
    b_diameter = 10+random.randint(0,9)/10 #diameter of beta cell: mu m
    i_area = area #10000 #diameter of islet: mu m^2
    perc_b = random.uniform(0.6,0.9)+random.uniform(0.0, 0.1) #percent of beta cells in the islet
    i_radius = (i_area/math.pi)**(1/2)

    #% Add edges

    #% Pack 3D grid with cells
    grid_rad = round(i_radius/b_diameter)
    points_inside_grid = product(range(-grid_rad,grid_rad),range(-grid_rad,grid_rad),range(-grid_rad,grid_rad))
    points_inside_grid = list(points_inside_grid)
    positions = [columns for columns in points_inside_grid if (columns[1]**2+columns[2]**2+columns[0]**2)**(1/2)<i_radius]
    
    beta_positions = np.array(positions)[np.random.choice(len(positions), size=int(perc_b*len(positions)), replace=False)]

    #1. Randomly pick edge distrutions
    #edgedist = np.random.normal(6.38+0.6, 1.35, len(beta_positions))
    edgedist = np.random.normal(6.38+0.6, 1.35, len(beta_positions)) 
 
    print('Number of Beta Cells' + str(len(edgedist)))
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
                ((v[0]-nodepos[0])**2+(v[1]-nodepos[1])**2+(v[2]-nodepos[2])**2)**(1/2)<thr 
                and k not in used_nodes and k != i and G.degree[k]<edgedist[k]-1]
            #choose which nodes to connect proximal edges based on degree
            
            if len(avaliable_edges) > int(edgedist[i])-G.degree[i]:
                edges_to_connect = np.random.choice(avaliable_edges, int(edgedist[i])-G.degree[i])
                [G.add_edge(i, ed) for ed in edges_to_connect] #add edges
                used_nodes.append(i)
            else:
                avaliable_edges = [k for (k,v) in positions_dict.items() if 
                    ((v[0]-nodepos[0])**2+(v[1]-nodepos[1])**2+(v[2]-nodepos[2])**2)**(1/2)<thr]
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
    gjdist = pd.read_csv('TotGJConductDistribution.csv', index_col = 0)
    gj_conduct = [float(item)*203 for item in list(gjdist)]
    gj_freq = list(gjdist.iloc[0,:])
    gj_freq_s = list(gjdist.iloc[2,:])
    #
    #x=np.linspace(1e-5,2,G.number_of_edges())
    #edgeweights_0 =  skew_norm_pdf(x,1.22/8, 0.5/4,-1.2)
    edgeweights_0 = np.random.normal(1.22/8, 0.5/3.5,G.number_of_edges())
    #edgeweights_0 = np.random.uniform(0, 1.22, np.size(x))
    #edgeweights_0 = npr.poisson(122/4, size=G.number_of_edges())/100
    #edgeweights_0 = npr.weibull(5, size=G.number_of_edges())/100
    #plt.hist(edgeweights_0)
    #print(edgeweights_0.mean())
    edgeweights_0 = [1e-5 if x<=0 else x for x in edgeweights_0] 
    
    # --- particle swarm ----- #
    #options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}

    #optimizer = GlobalBestPSO(n_particles=3000, dimensions=len(edgeweights_0), options=options, bounds=bounds)
    #bounds = ((1e-10, 3),)
    kwargs = {"G": G, "gjconduct" : gj_conduct, "gjfreq": gj_freq, "gjfreq_s": gj_freq_s}
    #cost, pos = optimizer.optimize(calc_conduct_optim, 1000, **kwargs)
    
    
    # -- final -- #
    bounds = op.Bounds(lb = 1e-5, ub = 1.22, keep_feasible=True)
    minimizer_kwargs = {"method":"L-BFGS-B", "bounds":bounds, "args":(G,  gj_conduct, gj_freq, gj_freq_s)}
    final_weights = op.basinhopping(calc_conduct_optim, edgeweights_0, minimizer_kwargs=minimizer_kwargs, niter=1000)#ftol=1e-20,xtol=1e-20,maxiter=1e10)
   # final_weights = op.minimize(calc_conduct_optim, edgeweights_0, args=(G,  gj_conduct, gj_freq, gj_freq_s), method='Nelder-Mead',bounds=bounds)#, ftol=1e-20,xtol=1e-20,maxiter=1e10)

    final_weights = final_weights.x
    err, Hist = calc_conduct_optim(final_weights, G, gj_conduct, gj_freq, gj_freq_s, options = 1)
    print(err)
    G2 = calc_conduct(final_weights, G, gj_conduct)
    
    for i in range(0, len(G2.nodes)):
        G2.nodes[i]['pos'] = str(G2.nodes[i]['pos'])
        G2.nodes[i]['Conduct'] = str(G2.nodes[i]['Conduct'])
    
    if err <= 6 + len(edgedist)/500:
        nx.write_gml(G2, 'size' + str(round(len(G.nodes))) + '_it' + str(itter) + '_er' + str(round(err,2))+'.gml')

if __name__ == "__main__":
    i = 1
    # for j in range(5000,1000,-1):
    #     i = i-1 
    #     run_networkbuild(i, j)
    for j in range(6000,500000,10):
        i = i+1 
        run_networkbuild(i, j)
# %%
