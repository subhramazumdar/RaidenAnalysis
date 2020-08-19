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
    return deg,cnt
    

def remove_random_node(G,x):
    
    x=int(G.number_of_nodes()*x/100)
    
    set_nodes=[y for y in G.nodes][:x]
    for node in set_nodes:
        G.remove_node(node)

def remove_capacity_node(G,x):
    
    x=int(G.number_of_nodes()*x/100)
    capacity_node=sorted(G.nodes(data=True), key=lambda x:x[1]['capacity'],reverse=True)[:x]
    
    for node in capacity_node:
        G.remove_node(node[0])
    
def remove_betweeness_node(G,x):
    x=int(G.number_of_nodes()*x/100)
    measure_of_centrality = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    nx.set_node_attributes(G, measure_of_centrality, 'betweenness_centrality')
    betweenness_centrality_node=sorted(G.nodes(data=True), key=lambda x:x[1]['betweenness_centrality'],reverse=True)[:x]
    
    for node in betweenness_centrality_node:
        G.remove_node(node[0])
    
def remove_degree_node(G,x):
    x=int(G.number_of_nodes()*x/100)
    degrees = dict(G.degree())
    nx.set_node_attributes(G, degrees, 'degree') 
    
    degree_sequence=sorted(G.nodes(data=True), key=lambda x:x[1]['degree'],reverse=True)[:x]
    for node in degree_sequence:
        G.remove_node(node[0])

def powerList(my_list, exponent):
    
    list_deg=[]
    
    for x in my_list:
        if x>0:
            list_deg.append(x**exponent)
        else:
            list_deg.append(0)
        
        
    return list_deg    
    
def main():
    with open(sys.argv[1]) as f:
        source = json.load(f)
        
        
    G=readGraph(source)
    print('total edges:', G.number_of_edges())
    print('total nodes:', G.number_of_nodes())
    set_graph_color(G)
    
    capacity_node(G)
    
    degree_avg(G)
    print("Density of RN: ",nx.classes.function.density(G))
    print("Number of connected components: ", len(list(nx.connected_components(G))))
    
    print("Maximal independent set: ", len(nx.algorithms.maximal_independent_set(G)))
    print("Number of bridges: ", len(list(nx.algorithms.bridges(G))))
    print("Size of the dominating set: ", len(nx.algorithms.dominating_set(G)))
    print("Raiden Network is Chordal graph: ", nx.algorithms.chordal.is_chordal(G))
    print("Raiden Network degree assortativity",nx.algorithms.assortativity.degree_assortativity_coefficient(G))
    
    #print("LN Wiener index", nx.algorithms.wiener_index(G)) #7686159.0
    #print("LN is Eulerian: ",nx.algorithms.is_eulerian(G))
    #print("LN is planar: ", nx.algorithms.planarity.check_planarity(G))
    #print("Number of isolates in LN: ", list(nx.isolates(G)))
   # print("LN's S-metric: ", smetric(G)) #0.6879664061934981
    print("RN average clustering coefficient", approximation.clustering_coefficient.average_clustering(G))
    print("RN's transitivity: ", nx.algorithms.cluster.transitivity(G))
    G_tmp=G.copy()
    #vertexConnectivity(G)
    
    c = max(nx.connected_components(G), key=len)
    G=G.subgraph(c).copy() 
    
    print("Raiden Network diameter: ", nx.algorithms.distance_measures.diameter(G)) #6
    print("Raiden Network radius", nx.algorithms.distance_measures.radius(G)) #3
    print("Average shortest paths: ",nx.algorithms.shortest_paths.generic.average_shortest_path_length(G)) # 2.806222412074612
    #G=G_tmp.copy()
    #edgeConnectivity(G)
    G=G_tmp.copy()
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
    print("Number of connected components\n")
    print(nx.number_connected_components(G))
    
    x0=[]
    y0=[]
    
    
    
    
    for i in range(5,45,5):
        G=G_tmp.copy()
        print("Removing "+str(i)+" high capacity nodes\n")
        remove_random_node(G,i)
        #plot_graph(G)
        a=nx.number_connected_components(G)
        print("Number of connected components\n")
        print(a)
        x0.append(int(i/100*G.number_of_nodes()))
        y0.append(a)
     
        
        
        
    x1=[]
    y1=[]
    
    
    
    
    for i in range(5,45,5):
        G=G_tmp.copy()
        print("Removing "+str(i)+" high capacity nodes\n")
        remove_capacity_node(G,i)
        #plot_graph(G)
        a=nx.number_connected_components(G)
        print("Number of connected components\n")
        print(a)
        x1.append(int(i/100*G.number_of_nodes()))
        y1.append(a)
     
        
        
    
    
    x2=[]
    y2=[]
    for i in range(5,45,5):
        G=G_tmp.copy()
        print("Removing "+str(i)+" high bc nodes\n")
        remove_betweeness_node(G,i)
        a=nx.number_connected_components(G)
        print("Number of connected components\n")
        print(a)
        x2.append(int(i/100*G.number_of_nodes()))
        y2.append(a)
     
        
        
        
        
        
        
        #plot_graph(G)
        #a=nx.number_connected_components(G)
        #print("Number of connected components\n")
        #print(a)
        #x2.append(i)
        #y2.append(a)
        #try:
         #   c = max(nx.connected_components(G), key=len)
          #  G=G.subgraph(c).copy() 
            
           # print("Raiden Network diameter: ", nx.algorithms.distance_measures.diameter(G)) #6
            #print("Raiden Network radius", nx.algorithms.distance_measures.radius(G)) #3
            #print("Average shortest paths: ",nx.algorithms.shortest_paths.generic.average_shortest_path_length(G))
        
        #except Exception:
         #   print("aborting")
               
    x3=[]
    y3=[]
    for i in range(5,45,5):
        G=G_tmp.copy()
        print("Removing "+str(i)+" high degree nodes\n")
        remove_degree_node(G,i)
        a=nx.number_connected_components(G)
        print("Number of connected components\n")
        print(a)
        x3.append(int(i/100*G.number_of_nodes()))
        y3.append(a)
     
        
        
        
        
        
         #   c = max(nx.connected_components(G), key=len)
          #  G=G.subgraph(c).copy() 
    
           # print("Raiden Network diameter: ", nx.algorithms.distance_measures.diameter(G)) #6
           # print("Raiden Network radius", nx.algorithms.distance_measures.radius(G)) #3
           # print("Average shortest paths: ",nx.algorithms.shortest_paths.generic.average_shortest_path_length(G))
        #except Exception:
         #  print("aborting")
    
    plt.plot(x0, y0, color='purple', markersize=3, label="Random ")
    plt.plot(x0, y1, color='red', markersize=3, label="High Capacity ")
    plt.plot(x0, y2, color='blue', markersize=3, label="High BC")
    plt.plot(x0, y3, color='green', markersize=3, label="High Degree")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    plt.title("Robustness of Network")
    plt.ylabel("Number of connected component")
    plt.xlabel("Number of nodes removed")
    
    
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
#    plt.show()

    plt.legend(loc="upper left")
    #plt.show()

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.savefig("plot_connected_components/"+sys.argv[1]+".png")

    f.close()
    
    #degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    #degreeCount = collections.Counter(degree_sequence)
    #deg, cnt = zip(*degreeCount.items())
   # print(deg)
   # print(cnt)
    #goodnessOfFit(deg, cnt)
    #percolationThresholdPrediction(deg)
    # print "Degree sequence", degree_sequence
    
    #fittingGamma(deg, cnt)
    
if __name__ == '__main__':
    main()
    
    
