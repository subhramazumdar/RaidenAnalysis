import networkx as nx
from statistics import mean 
import logging

def readGraph(source):
    
    G = nx.Graph()
    
        
    #    print('Adding node with pubkey', node['pub_key'])

    # print(G.nodes)
    # 
#    for node in source['nodes']:
 #       G.add_node(node)
  #      print(node)

        

            
    for edge in source['channels']:
        G.add_edge(edge['participant1'].lower(),edge['participant2'].lower(), node1_deposit=edge['deposit1'],node2_deposit=edge['deposit2'])
    
    for node in source['nodes']:
      if node.lower() not in G.nodes:
	        G.add_node(node.lower())
        	

    print(G.number_of_nodes())  
    for node in G.nodes:
        print(node)
              
    return G

def capacity_node(G):
    
    nx.set_node_attributes(G,0,'capacity')
    for edge in G.edges:
        
        G.nodes[edge[0]]['capacity']+=G.edges[edge]['node1_deposit']
        G.nodes[edge[1]]['capacity']+=G.edges[edge]['node2_deposit']
        G.edges[edge]['capacity']=G.edges[edge]['node1_deposit']+G.edges[edge]['node2_deposit']
        G.edges[edge]['start']=edge[0]
        G.edges[edge]['end']=edge[1]
        
    #print("Total DAI Stablecoin held in the network")
    return sum([G.nodes[node]['capacity'] for node in G.nodes])
    
def capacity_node_main(G,G_part,num_arg):
    
    nx.set_node_attributes(G,0,'capacity')
    included=[]
    for i in range(1,num_arg):
        for edge in G.edges:
            if edge not in included and edge in G_part[num_arg-i].edges:
                G.nodes[edge[0]]['capacity']+=G_part[num_arg-i].edges[edge]['node1_deposit']
                G.nodes[edge[1]]['capacity']+=G_part[num_arg-i].edges[edge]['node2_deposit']
                included.append(edge)
               
        
    #print("Total DAI Stablecoin held in the network")
    return sum([G.nodes[node]['capacity'] for node in G.nodes])

def degree_avg(G):
    
    degrees = dict(G.degree())
    nx.set_node_attributes(G, degrees, 'degree')
    
    #for n in G.degree:
        #if n[1]>=10:
         #   print(n)
#    print("Average Degree: ")
    return mean([n[1] for n in G.degree])
