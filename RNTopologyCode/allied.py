import networkx as nx
from statistics import mean 
import logging

def readGraph(source):
    
    G = nx.Graph()
    
        
    #    print('Adding node with pubkey', node['pub_key'])

    # print(G.nodes)
    # 

    for edge in source['channels']:
        G.add_edge(edge['participant1'],edge['participant2'], node1_deposit=edge['deposit1'],node2_deposit=edge['deposit2'])

            
    return G

def capacity_node(G):
    
    nx.set_node_attributes(G,0,'capacity')
    for edge in G.edges:
        G.nodes[edge[0]]['capacity']+=G.edges[edge]['node1_deposit']
        G.nodes[edge[1]]['capacity']+=G.edges[edge]['node2_deposit']
        G.edges[edge]['capacity']=G.edges[edge]['node1_deposit']+G.edges[edge]['node2_deposit']
        
    #print("Total DAI Stablecoin held in the network")
    return sum([G.nodes[node]['capacity'] for node in G.nodes])
    

def degree_avg(G):
    
    degrees = dict(G.degree())
    nx.set_node_attributes(G, degrees, 'degree')
    
    for n in G.degree:
        if n[1]>=10:
            print(n)
#    print("Average Degree: ")
    return mean([n[1] for n in G.degree])
