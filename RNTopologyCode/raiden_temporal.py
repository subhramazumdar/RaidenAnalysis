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

def bc_aggregated(G):
    
    measure_of_centrality = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    nx.set_node_attributes(G, measure_of_centrality, 'betweenness_centrality')
    
    return [(v,G.nodes[v]['betweenness_centrality']) for v in G.nodes]
    
def bc_average(G,num_arg):
    
    for i in range(1,num_arg):
        measure_of_centrality = nx.betweenness_centrality(G[i], k=None, normalized=True, weight=None, endpoints=False, seed=None)
        nx.set_node_attributes(G[i], measure_of_centrality, 'betweenness_centrality')
        
        
    

    
    
def bc_temporal(G,G_main,val,shortest_path_matrix):
    
    bc={}
    set_return=[]
    for k in G_main.nodes:
        bc[k]=0
        
    for i in range(1,val):
            #print(i)
            for k in G_main.nodes:
                #print(k)
                for u in G_main.nodes:
                    for v in G_main.nodes:
                        if u!=v and v!=k and u!=k:
                            
                            a=shortest_path_matrix[val][(u,v)][i]
                            if a>0:
                                #print(a)
                                for j in range(i+1,val+1):
                                
                                    bc[k]=bc[k]+(shortest_path_matrix[j-1][(u,k)][i]*shortest_path_matrix[val][(k,v)][j])/a
                                    #print(u,end=' ')
                                    #print(v,end=' ')
                                    #print(shortest_path_matrix[j-1][(u,k)][i]*shortest_path_matrix[val][(k,v)][j],end=' ')
                                    #print(a)
                            
                            
    n=G_main.number_of_nodes()                        
    for k in G_main.nodes:
        set_return.append((k,2*bc[k]/(val*(n-1)*(n-2))))
                            
    return set_return

def compute(mat,shortest,G_main,G,val):
    
    for u in G_main.nodes:
        for v in G_main.nodes:
           mat[(u,v)]=dict.fromkeys(range(1,val),99999)
           shortest[(u,v)]=dict.fromkeys(range(1,val),0)                  
        
    for v in G_main.nodes:
        
           mat[(v,v)][val-1]=1
           shortest[(v,v)][val-1]=1
           
        
    for e in G[val-1].edges:
        
        shortest[(e[0],e[1])][val-1]=1
        shortest[(e[1],e[0])][val-1]=1
        mat[(e[0],e[1])][val-1]=1
        mat[(e[1],e[0])][val-1]=1
        
    
    for i in range(val-2,0,-1):
        
        for u in G_main.nodes:
            
            for v in G_main.nodes:
                    
                #if (u,v) in G[i].edges:
                 #   matrix[(u,v)][i]=1
                #else:    
                    for k in G_main.nodes:
                        if (u==k or (u,k) in G[i].edges) and mat[(u,v)][i]>=mat[(k,v)][i+1]+1:
                            if mat[(u,v)][i]>mat[(k,v)][i+1]+1:
                                mat[(u,v)][i]=mat[(k,v)][i+1]+1
                            #print("u",end=' ')
                            #print(u,end=' ')                    
                            #print("k",end=' ')                
                            #print(k,end=' ')
                            #print("v",end=' ')
                            #print(v)
                            
                            shortest[(u,v)][i]=shortest[(u,v)][i]+shortest[(k,v)][i+1]
                            
                        
                 
    return mat,shortest

def temporal_degree(G,num_arg,G_main):
    
    return_set=[]
    for i in range(1,num_arg):
        degrees = dict(G[i].degree())
        nx.set_node_attributes(G[i], degrees, 'degree') 
        degrees.clear()
    
    for v in G_main.nodes:
        val=0
        for i in range(1,num_arg):
            if v in G[i].nodes:
                val=val+G[i].nodes[v]['degree']
                
        return_set.append((v,val/((num_arg-1)*(G_main.number_of_nodes()-1))))    
    return return_set
    
        
def aggregated_degree(G,num_arg):
    
    return_set=[]
    degrees = dict(G.degree())
    nx.set_node_attributes(G, degrees, 'degree') 
    
    for v in G.nodes:
        return_set.append((v,G.nodes[v]['degree']/(G.number_of_nodes()-1)))
        
    return return_set    


def plot_graph(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,str1,str2):
    
    x0=[]
    y0=[]
    
    
    
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print("Number of connected components\n")
        print(a)
        x0.append(k)
        y0.append(a)
     
        
        
        
    x1=[]
    y1=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print("Number of connected components\n")
        print(a)
        x1.append(k)
        y1.append(a)
     
     
    x2=[]
    y2=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print("Number of connected components\n")
        print(a)
        x2.append(k)
        y2.append(a)
     

    x3=[]
    y3=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print("Number of connected components\n")
        print(a)
        x3.append(k)
        y3.append(a)
     

    plt.plot(x0, y0, 'r^--',markersize=8, label="Temporal ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="Aggregated ")
    plt.plot(x0, y2, 'yo--', label="Average ")
    plt.plot(x0, y3, 'gs--', label="Normal")
    
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

    plt.savefig("plot_temporal_connected_components/test"+str1+"-"+str2+".png")
    plt.close()

def plot_capacity(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,str1,str2):
    
    x0=[]
    y0=[]
    
    
    capacity_init=sum([G.nodes[x]['capacity'] for x in G.nodes])
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.nodes[x]['capacity'] for x in G_tmp.nodes])
        percentage_dec=(capacity_init-capacity_now)*100/capacity_init    
        x0.append(k)
        y0.append(percentage_dec)
     
        
        
        
    x1=[]
    y1=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.nodes[x]['capacity'] for x in G_tmp.nodes])
        percentage_dec=(capacity_init-capacity_now)*100/capacity_init    
        x1.append(k)
        y1.append(percentage_dec)
     
     
    x2=[]
    y2=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.nodes[x]['capacity'] for x in G_tmp.nodes])
        percentage_dec=(capacity_init-capacity_now)*100/capacity_init    
        x2.append(k)
        y2.append(percentage_dec)
     
    x3=[]
    y3=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.nodes[x]['capacity'] for x in G_tmp.nodes])
        percentage_dec=(capacity_init-capacity_now)*100/capacity_init    
        x3.append(k)
        y3.append(percentage_dec)
     

    
    plt.plot(x0, y0, 'r^--', markersize=8,label="Temporal ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="Aggregated ")
    plt.plot(x0, y2, 'yo--', label="Average ")
    plt.plot(x0, y3, 'gs--', label="Normal ")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    plt.title("Robustness of Network")
    plt.ylabel("Loss in Capacity")
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

    plt.savefig("plot_capacity_connected_components/test"+str1+"-"+str2+".png")
    plt.close()
    
    
def remove_degree_node(G,setlist,count):
    

    x=int(G.number_of_nodes()*count/100)
    
    #print(setlist)
    degree_sequence=sorted(setlist, key=lambda x:x[1],reverse=True)[:x]
    print(degree_sequence)
    
    for node in degree_sequence:
        
        G.remove_node(node[0])

    
    
def plot_shortest(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,str1,str2):
    x0=[]
    y0=[]
    
    
    
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        
        
        x0.append(k)
        
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y0.append(a)
        
        except Exception:
            y0.append(0)
            print("aborting") 
        
        
        
    x1=[]
    y1=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        
        
        x1.append(k)
        
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y1.append(a)
        
        except Exception:
            y1.append(0)
            print("aborting")
    x2=[]
    y2=[]
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        
        
        x2.append(k)
        
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y2.append(a)
        
        except Exception:
            y2.append(0)
            print("aborting")
    x3=[]
    y3=[]
    
    
    
    
    for i in range(5,45,5):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        
        
        x3.append(k)
        
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y3.append(a)
        
        except Exception:
            y3.append(0)
            print("aborting") 


    
    plt.plot(x0, y0, 'r^--', markersize=8,label="Temporal ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="Aggregated ")
    plt.plot(x0, y2, 'yo--', label="Average ")
    plt.plot(x0, y3, 'gs--', label="Normal ")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    plt.title("Robustness of Network")
    plt.ylabel("Average Shortest Path length")
    plt.xlabel("Number of nodes removed")
    
    
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
#    plt.show()

    plt.legend(loc="upper right")
    #plt.show()

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.savefig("plot_shortest_connected_components/test"+str1+"-"+str2+".png")
    plt.close()
def Distribution(set_temp_deg,set_agg_deg,set_avg_deg,str1,str2):
    
    
    # print "Degree sequence", degree_sequence
    
    tempCount = collections.Counter(set_temp_deg)
    deg, cnt = zip(*tempCount.items())

    aggCount=collections.Counter(set_agg_deg)
    deg1, cnt1 = zip(*aggCount.items())
    
    avgCount=collections.Counter(set_avg_deg)
    deg2, cnt2 = zip(*avgCount.items())
    
    
    plt.plot(np.array(deg), np.array(cnt), 'r^-',markersize=7, label = 'Temporal distribution')
    plt.plot(np.array(deg1), np.array(cnt1), 'bs-', markersize=8, label = 'Aggregated distribution')
    plt.plot(np.array(deg2), np.array(cnt2), 'go-', label = 'Average distribution')
    
    plt.title(str1)
    plt.ylabel("Count")
    plt.xlabel("measure_of_centrality")
    #ax.set_xticks([d + 0.4 for d in deg])
    #ax.set_xticklabels(deg)
    #leg = ax.legend();

    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
    plt.legend(loc="upper right")
    plt.savefig(str1+"-"+str2+".png")
    plt.close()


        
def main():
    
    G_main=nx.Graph()
    num_arg=len(sys.argv)-1
    source={}
    set_nodes=[]
    G={}
    for i in range(1,num_arg):
    
        with open(sys.argv[i]) as f:
            #print(sys.argv[i])
            source[i] = json.load(f)
        
            G[i]=readGraph(source[i])
            print(G[i].number_of_nodes())
            for v in G[i].nodes:
                if v not in set_nodes: 
                    set_nodes.append(v)
                    
            capacity_node(G[i])        
        f.close()            
    
    bet_cen={}
    cb={}
    cb_avg={}
    cb_new=dict()
    
    
    matrix = dict()
    shortest_path_matrix=dict()
    
    
    
    for v in set_nodes:
        cb_new[v]=dict.fromkeys(range(1,num_arg),0)
        
    for v in set_nodes:
        cb[v]=0
        cb_avg[v]=0
        bet_cen[v]=0
        G_main.add_node(v)

    set_edge=[]        
    
    for i in range(1,num_arg):
        for e in G[i].edges:
            if e not in G_main.edges:
                
                G_main.add_edge(e[0],e[1])
                
    
    for i in range(1,num_arg):
        matrix[i]={}
        shortest_path_matrix[i]={}
        
        matrix[i],shortest_path_matrix[i]=compute(matrix[i],shortest_path_matrix[i],G_main,G,i+1)    
    
    #k=max([G[i].number_of_nodes() for i in range(1,7)])
    #print(k)
    
    
    
                    
    #for u in G_main.nodes:
     #   print(u)
      #  for v in G_main.nodes:
       #     for i in range(1,num_arg):
        #        print(shortest_path_matrix[num_arg-1][(u,v)][i],end=' ')
         #   print("\n")    
    
    for v in G_main.nodes:                    
        for i in range(1,num_arg):
            
            for k in G_main.nodes:
            
                if k!=v and matrix[num_arg-1][(v,k)][i]<99999:
                    
                    cb[v]=cb[v]+(1/matrix[num_arg-1][(v,k)][i])
        
    cc_temporal=[]
    for v in G_main.nodes:
        #print(cb[v])
        cb[v]=cb[v]/((num_arg-1)*(G_main.number_of_nodes()-1))
        cc_temporal.append((v,cb[v]))
    
    print(cc_temporal)
    for v in G_main.nodes:
        for k in G_main.nodes:
            
            if v!=k and nx.has_path(G_main,v,k):
                
                cb_avg[v]=cb_avg[v]+(1/(len(nx.shortest_path(G_main,v,k))-1))
    
    print("\n\n")      
    #print(G_main.number_of_nodes())
    #print(G_main.number_of_edges())
    
    cc_agg=[]
    for v in G_main.nodes:
     
        cb_avg[v]=cb_avg[v]/(G_main.number_of_nodes()-1)         
        cc_agg.append((v,cb_avg[v]))
    
    
    print(cc_agg)
    for i in range(1,num_arg):
        
        for v in G[i].nodes:
            for k in G[i].nodes:
            
                if v!=k and nx.has_path(G[i],v,k):
                    lenpath=len(nx.shortest_path(G[i],v,k))-1
                    cb_new[v][i]=cb_new[v][i]+(1/lenpath)
                    
            
            
            cb_new[v][i]=cb_new[v][i]/(G[i].number_of_nodes()-1)
            #print(cb_new[v][i])
        
    print("\n\n")        
    cc_average=[]
    for v in G_main.nodes:
        val=0
        for i in range(1,num_arg):
            val=val+cb_new[v][i]
            
        val=val/(num_arg-1)
        cc_average.append((v,val))
        
    
    print(cc_average)
    print("\n\n") 
    bc_ag=bc_aggregated(G_main)
    for node in G_main.nodes:
        print(G_main.nodes[node]['betweenness_centrality'])
    print("\n\n")         
    bc_average(G,num_arg)
    
    set_bc=[]
    for node in G_main.nodes:
        val=0
        for i in range(1,num_arg):
            if node in G[i].nodes:
                val=val+G[i].nodes[node]['betweenness_centrality']
        set_bc.append((node,val/(num_arg-1)))        
        
    
    
    
        
    bc_tmp=bc_temporal(G,G_main,num_arg-1,shortest_path_matrix)
    
    set_temp_deg=temporal_degree(G,num_arg,G_main)
    set_agg_deg=aggregated_degree(G_main,num_arg)
    
    #plotDeg(set_temp_deg,set_agg_deg,set_temp_deg)
    
  #  print(bc_tmp)
  #  print(bc_ag)
  #  print(set_bc)
    
    #plot_graph(G[num_arg-1],set_temp_deg,set_agg_deg,set_temp_deg)
    #print(bc_tmp)
    #print(bc_ag)
    #plot_graph(G[num_arg-1],bc_tmp,bc_ag,set_bc)
    
    measure_of_centrality = nx.betweenness_centrality(G[num_arg-1], k=None, normalized=True, weight=None, endpoints=False, seed=None)
    nx.set_node_attributes(G[num_arg-1], measure_of_centrality, 'betweenness_centrality')
    betweenness_centrality_node=[]
    
    for v in G[num_arg-1].nodes:
        betweenness_centrality_node.append((v,G[num_arg-1].nodes[v]['betweenness_centrality']))
    
    measure_of_centrality1 = nx.closeness_centrality(G[num_arg-1],u=None, distance=None, wf_improved=True)
    nx.set_node_attributes(G[num_arg-1], measure_of_centrality1, 'closeness_centrality')
    closeness_centrality_node=[]
    
    for v in G[num_arg-1].nodes:
        closeness_centrality_node.append((v,G[num_arg-1].nodes[v]['closeness_centrality']))
        
    
    degree_node=[]
    
    for v in G[num_arg-1].nodes:
        degree_node.append((v,G[num_arg-1].nodes[v]['degree']))
    
    
    plot_graph(G[num_arg-1],set_temp_deg,set_agg_deg,set_temp_deg,degree_node,"degree",sys.argv[num_arg])
    plot_graph(G[num_arg-1],bc_tmp,bc_ag,set_bc,betweenness_centrality_node,"betweenness",sys.argv[num_arg])
    plot_graph(G[num_arg-1],cc_temporal,cc_agg,cc_average,closeness_centrality_node,"closeness",sys.argv[num_arg])
    
    #plotDeg(bc_tmp,bc_ag,set_bc)
    plot_shortest(G[num_arg-1],set_temp_deg,set_agg_deg,set_temp_deg,degree_node,"degree",sys.argv[num_arg])
    plot_shortest(G[num_arg-1],bc_tmp,bc_ag,set_bc,betweenness_centrality_node,"betweenness",sys.argv[num_arg])
    plot_shortest(G[num_arg-1],cc_temporal,cc_agg,cc_average,closeness_centrality_node,"closeness",sys.argv[num_arg])
    
        
        
    plot_capacity(G[num_arg-1],set_temp_deg,set_agg_deg,set_temp_deg,degree_node,"degree",sys.argv[num_arg])
    
    
    
    plot_capacity(G[num_arg-1],bc_tmp,bc_ag,set_bc,betweenness_centrality_node,"betweenness",sys.argv[num_arg])
    plot_capacity(G[num_arg-1],cc_temporal,cc_agg,cc_average,closeness_centrality_node,"closeness",sys.argv[num_arg])
    
      # degree sequence
    Distribution(sorted([d[1] for d in set_temp_deg], reverse=True),sorted([d[1] for d in set_agg_deg], reverse=True),sorted([d[1] for d in set_temp_deg], reverse=True),"Degree Distribution",sys.argv[num_arg])
    Distribution(sorted([d[1] for d in bc_tmp], reverse=True),sorted([d[1] for d in bc_ag], reverse=True),sorted([d[1] for d in set_bc], reverse=True),"Betweeness Centrality Distribution",sys.argv[num_arg])
    Distribution(sorted([d[1] for d in cc_temporal], reverse=True),sorted([d[1] for d in cc_agg], reverse=True),sorted([d[1] for d in cc_average], reverse=True),"Closeness Centrality Distribution",sys.argv[num_arg])
    
if __name__ == '__main__':
    main()
    
    
