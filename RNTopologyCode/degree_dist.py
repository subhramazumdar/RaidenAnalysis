import json
import logging
import networkx as nx
import sys
from allied import readGraph,capacity_node,degree_avg
import heapq
from networkx.algorithms import approximation
import collections
from scipy import special as spcl
import random
import numpy as np
from random import randint
import operator

from graph_visualize import plot_graph, set_graph_color, set_node_color

import matplotlib.pyplot as plt
def degreeDistribution(G):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    
    
    plt.plot(np.array(deg), np.array(cnt), color='red', label = 'Empirical degree distribution')
   

    plt.title("Degree Distribution")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)
    leg = ax.legend();


    plt.yscale("log")
    plt.xscale("log")
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
    plt.savefig("plot_degree/"+sys.argv[1]+".png")

    
def main():
    with open(sys.argv[1]) as f:
        source = json.load(f)
        
    f1=open(sys.argv[2],"a")    
    G=readGraph(source)
    print('total edges:', G.number_of_edges())
    print('total nodes:', G.number_of_nodes())
    set_graph_color(G)
    
    cap=capacity_node(G)
    
    deg=degree_avg(G)
    
    f1.write(sys.argv[1]+" "+str(G.number_of_nodes())+" "+str(G.number_of_edges())+" "+str(cap)+" "+str(deg)+" "+str(nx.classes.function.density(G))+" "+str(len(list(nx.connected_components(G))))+" "+str(len(nx.algorithms.maximal_independent_set(G)))+" "+ str(len(list(nx.algorithms.bridges(G))))+" "+str(len(nx.algorithms.dominating_set(G)))+" "+str(nx.algorithms.chordal.is_chordal(G))+" "+str(nx.algorithms.assortativity.degree_assortativity_coefficient(G))+" "+str(approximation.clustering_coefficient.average_clustering(G))+" "+str(nx.algorithms.cluster.transitivity(G))+" ")
    #print("Density of RN: ",nx.classes.function.density(G))
    #print("Number of connected components: ", len(list(nx.connected_components(G))))
    
    #print("Maximal independent set: ", len(nx.algorithms.maximal_independent_set(G)))
    #print("Number of bridges: ", len(list(nx.algorithms.bridges(G))))
    #print("Size of the dominating set: ", len(nx.algorithms.dominating_set(G)))
    #print("Raiden Network is Chordal graph: ", nx.algorithms.chordal.is_chordal(G))
    #print("Raiden Network degree assortativity",nx.algorithms.assortativity.degree_assortativity_coefficient(G))
    
    #print("LN Wiener index", nx.algorithms.wiener_index(G)) #7686159.0
    #print("LN is Eulerian: ",nx.algorithms.is_eulerian(G))
    #print("LN is planar: ", nx.algorithms.planarity.check_planarity(G))
    #print("Number of isolates in LN: ", list(nx.isolates(G)))
   # print("LN's S-metric: ", smetric(G)) #0.6879664061934981
   # print("RN average clustering coefficient", approximation.clustering_coefficient.average_clustering(G))
    #print("RN's transitivity: ", nx.algorithms.cluster.transitivity(G))
    degreeDistribution(G)
    c = max(nx.connected_components(G), key=len)
    G=G.subgraph(c).copy() 
    f1.write(str(nx.algorithms.distance_measures.diameter(G))+" "+str(nx.algorithms.distance_measures.radius(G))+" "+str(nx.algorithms.shortest_paths.generic.average_shortest_path_length(G))+"\n")
    
    #centrality(G)
    #print("RN's largest clique size: ", nx.algorithms.approximation.clique.max_clique(G))
    #print("Adjacency spectrum of RN: ", nx.linalg.spectrum.adjacency_spectrum(G))
    #G=G_tmp.copy()
   
    
    #print(deg)
    #print(cnt)
    #goodnessOfFit(deg, cnt)
    #percolationThresholdPrediction(deg)
    #attackingBetweenness(G,deg,cnt)
    #plot_graph(G)
        
if __name__ == '__main__':
    main()
    
    
