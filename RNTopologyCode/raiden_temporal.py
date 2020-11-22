import json
import logging
import numpy as np
from bokeh.palettes import plasma
import networkx as nx
import sys
from allied import readGraph,capacity_node,degree_avg,capacity_node_main
import heapq
from networkx.algorithms import approximation,s_metric
import collections
from scipy import special as spcl
import random
import numpy as np
from random import randint
import operator

from random import seed

from graph_visualize import plot_graph1, set_graph_color, set_node_color

import math

import types

from mpmath import mp
from decimal import *
import string

import matplotlib.pyplot as plt

def aggregated_capacity(G):
    
    
    return [(v,G.nodes[v]['capacity']) for v in G.nodes]

def temporal_capacity(G,G_main,num_arg):
    return_set=[]
    for v in G_main.nodes:
        val=0
        for i in range(1,num_arg):
            if v in G[i].nodes:
                val=val+G[i].nodes[v]['capacity']
                
        return_set.append((v,val/((num_arg-1)*(G_main.number_of_nodes()-1))))    
    return return_set

def bc_aggregated(G):
    
    measure_of_centrality = nx.betweenness_centrality(G, k=None, normalized=True, weight=None, endpoints=False, seed=None)
    nx.set_node_attributes(G, measure_of_centrality, 'betweenness_centrality')
    
    return [(v,G.nodes[v]['betweenness_centrality']) for v in G.nodes]
    
def bc_average(G,num_arg):
    
    for i in range(1,num_arg):
        measure_of_centrality = nx.betweenness_centrality(G[i], k=None, normalized=True, weight=None, endpoints=False, seed=None)
        nx.set_node_attributes(G[i], measure_of_centrality, 'betweenness_centrality')
        
        
    

def li_smax_graph(degree_seq):
    """Generates a graph based with a given degree sequence and maximizing
    the s-metric.  Experimental implementation.
    Maximum s-metrix  means that high degree nodes are connected to high
    degree nodes.
    - `degree_seq`: degree sequence, a list of integers with each entry
       corresponding to the degree of a node.
       A non-graphical degree sequence raises an Exception.
    Reference::
      @unpublished{li-2005,
       author = {Lun Li and David Alderson and Reiko Tanaka
                and John C. Doyle and Walter Willinger},
       title = {Towards a Theory of Scale-Free Graphs:
               Definition, Properties, and  Implications (Extended Version)},
       url = {http://arxiv.org/abs/cond-mat/0501169},
       year = {2005}
      }
    The algorithm::
     STEP 0 - Initialization
     A = {0}
     B = {1, 2, 3, ..., n}
     O = {(i; j), ..., (k, l),...} where i < j, i <= k < l and
             d_i * d_j >= d_k *d_l
     wA = d_1
     dB = sum(degrees)
     STEP 1 - Link selection
     (a) If |O| = 0 TERMINATE. Return graph A.
     (b) Select element(s) (i, j) in O having the largest d_i * d_j , if for
             any i or j either w_i = 0 or w_j = 0 delete (i, j) from O
     (c) If there are no elements selected go to (a).
     (d) Select the link (i, j) having the largest value w_i (where for each
             (i, j) w_i is the smaller of w_i and w_j ), and proceed to STEP 2.
     STEP 2 - Link addition
     Type 1: i in A and j in B.
             Add j to the graph A and remove it from the set B add a link
             (i, j) to the graph A. Update variables:
             wA = wA + d_j -2 and dB = dB - d_j
             Decrement w_i and w_j with one. Delete (i, j) from O
     Type 2: i and j in A.
         Check Tree Condition: If dB = 2 * |B| - wA.
             Delete (i, j) from O, continue to STEP 3
         Check Disconnected Cluster Condition: If wA = 2.
             Delete (i, j) from O, continue to STEP 3
         Add the link (i, j) to the graph A
         Decrement w_i and w_j with one, and wA = wA -2
     STEP 3
         Go to STEP 1
    The article states that the algorithm will result in a maximal s-metric.
    This implementation can not guarantee such maximality. I may have
    misunderstood the algorithm, but I can not see how it can be anything
    but a heuristic. Please contact me at sundsdal@gmail.com if you can
    provide python code that can guarantee maximality.
    Several optimizations are included in this code and it may be hard to read.
    Commented code to come.
    """

    if not is_valid_degree_sequence(degree_seq):
        raise nx.NetworkXError('Invalid degree sequence')
    degree_seq.sort()  # make sure it's sorted
    degree_seq.reverse()
    degrees_left = degree_seq[:]
    A_graph = nx.Graph()
    A_graph.add_node(0)
    a_list = [False] * len(degree_seq)
    b_set = set(range(1, len(degree_seq)))
    a_open = set([0])
    O = []
    for j in b_set:
        heapq.heappush(O, (-degree_seq[0] * degree_seq[j], (0, j)))
    wa = degrees_left[0]  # stubs in a_graph
    db = sum(degree_seq) - degree_seq[0]  # stubs in b-graph
    a_list[0] = True  # node 0 is now in a_Graph
    bsize = len(degree_seq) - 1  # size of b_graph
    selected = []
    weight = 0
    while O or selected:
        if len(selected) < 1:
            firstrun = True
            while O:
                (newweight, (i, j)) = heapq.heappop(O)
                if degrees_left[i] < 1 or degrees_left[j] < 1:
                    continue
                if firstrun:
                    firstrun = False
                    weight = newweight
                if not newweight == weight:
                    break
                heapq.heappush(selected, [-degrees_left[i], \
                                          -degrees_left[j], (i, j)])
            if not weight == newweight:
                heapq.heappush(O, (newweight, (i, j)))
            weight *= -1
        if len(selected) < 1:
            break

        [w1, w2, (i, j)] = heapq.heappop(selected)
        if degrees_left[i] < 1 or degrees_left[j] < 1:
            continue
        if a_list[i] and j in b_set:
            # TYPE1
            a_list[j] = True
            b_set.remove(j)
            A_graph.add_node(j)
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            wa += degree_seq[j] - 2
            db -= degree_seq[j]
            bsize -= 1
            newweight = weight
            if not degrees_left[j] == 0:
                a_open.add(j)
                for k in b_set:
                    if A_graph.has_edge(j, k): continue
                    w = degree_seq[j] * degree_seq[k]
                    if w > newweight:
                        newweight = w
                    if weight == w and not newweight > weight:
                        heapq.heappush(selected, [-degrees_left[j], \
                                                  -degrees_left[k], (j, k)])
                    else:
                        heapq.heappush(O, (-w, (j, k)))
                if not weight == newweight:
                    while selected:
                        [w1, w2, (i, j)] = heapq.heappop(selected)
                        if degrees_left[i] * degrees_left[j] > 0:
                            heapq.heappush(O, [-degree_seq[i] * degree_seq[j], (i, j)])
            if degrees_left[i] == 0:
                a_open.discard(i)

        else:
            # TYPE2
            if db == (2 * bsize - wa):
                # tree condition
                # print "removing because tree condition    "
                continue
            elif db < 2 * bsize - wa:
                raise nx.NetworkXError("THIS SHOULD NOT HAPPEN! - not graphable")
                continue
            elif wa == 2 and bsize > 0:
                # print "removing because disconnected  cluster"
                # disconnected cluster condition
                continue
            elif wa == db - (bsize) * (bsize - 1):
                # print "MYOWN removing because disconnected  cluster"
                continue
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            if degrees_left[i] < 1:
                a_open.discard(i)
            if degrees_left[j] < 1:
                a_open.discard(j)
            wa -= 2
            if not degrees_left[i] < 0 and not degrees_left[j] < 0:
                selected2 = (selected)
                selected = []
                while selected2:
                    [w1, w1, (i, j)] = heapq.heappop(selected2)
                    if degrees_left[i] * degrees_left[j] > 0:
                        heapq.heappush(selected, [-degrees_left[i], \
                                                  -degrees_left[j], (i, j)])
    return A_graph    

def is_valid_degree_sequence(deg_sequence):
    """Return True if deg_sequence is a valid sequence of integer degrees
    equal to the degree sequence of some simple graph.
      - `deg_sequence`: degree sequence, a list of integers with each entry
         corresponding to the degree of a node (need not be sorted).
         A non-graphical degree sequence (i.e. one not realizable by some
         simple graph) will raise an exception.
    See Theorem 1.4 in [chartrand-graphs-1996]. This algorithm is also used
    in havel_hakimi_graph()
    References:
    [chartrand-graphs-1996] G. Chartrand and L. Lesniak, "Graphs and Digraphs",
                            Chapman and Hall/CRC, 1996.
    """
    # some simple tests
    if deg_sequence == []:
        return True  # empty sequence = empty graph
    if not nx.utils.is_list_of_ints(deg_sequence):
        return False  # list of ints
    if min(deg_sequence) < 0:
        return False  # each int not negative
    if sum(deg_sequence) % 2:
        return False  # must be even

    # successively reduce degree sequence by removing node of maximum degree
    # as in Havel-Hakimi algorithm

    s = deg_sequence[:]  # copy to s
    while s:
        s.sort()  # sort in non-increasing order
        if s[0] < 0:
            return False  # check if removed too many from some node

        d = s.pop()  # pop largest degree
        if d == 0: return True  # done! rest must be zero due to ordering

        # degree must be <= number of available nodes
        if d > len(s):   return False

        # remove edges to nodes of next higher degrees
        s.reverse()  # to make it easy to get at higher degree nodes.
        for i in range(d):
            s[i] -= 1

    # should never get here b/c either d==0, d>len(s) or d<0 before s=[]
    return False    
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

def bc_edge_temporal(G,G_main,val,shortest_path_matrix):
    
    bc={}
    set_return=[]
    for e in G_main.edges:
        bc[e]=0
        
    for i in range(1,val):
            #print(i)
            for k in G_main.edges:
                #print(k)
                for u in G_main.nodes:
                    for v in G_main.nodes:
                        if u!=v and u!=k[0] and v!=k[1]:
                            
                            a=shortest_path_matrix[val][(u,v)][i]
                            if a>0:
                                #print(a)
                                for j in range(i+1,val):
                                
                                    bc[e]=bc[e]+(shortest_path_matrix[j-1][(u,k[0])][i]*shortest_path_matrix[val][(k[1],v)][j+1])/a
                                    #print(u,end=' ')
                                    #print(v,end=' ')
                                    #print(shortest_path_matrix[j-1][(u,k)][i]*shortest_path_matrix[val][(k,v)][j],end=' ')
                                    #print(a)
                            
                            
    n=G_main.number_of_nodes()                        
    for k in G_main.edges:
        set_return.append((k,2*bc[k]/(val*(n-1)*(n))))
                            
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


    
def plot_graph(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,randomlist,str1,str2):
    
    x0=[]
    y0=[]
    z0=[]
    
    

    print("temporal\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        #print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        #a=nx.number_connected_components(G_tmp)
        print(a)
        x0.append(k)
        y0.append(a)
        z0.append(val)
        
        
        
    x1=[]
    y1=[]
    z1=[]
    print("aggregated\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        #print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x1.append(k)
        y1.append(a)
        z1.append(val)
        """
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x1.append(k)
        y1.append(a)
        """
     
     
    x2=[]
    y2=[]
    z2=[]
    
    print("average\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        #a=nx.number_connected_components(G_tmp)
        #print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x2.append(k)
        y2.append(a)
        z2.append(val)
        """
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x2.append(k)
        y2.append(a)
        """

    x3=[]
    y3=[]
    z3=[]
    print("normal\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        #a=nx.number_connected_components(G_tmp)
        #print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x3.append(k)
        y3.append(a)
        z3.append(val)
        """
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        a=nx.number_connected_components(G_tmp)
        print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x3.append(k)
        y3.append(a)
        """
    x4=[]
    y4=[]
    z4=[]
    print("random\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        
        
        node_list = np.array(list(G_tmp.nodes()))
          
        k=int(i/100*G_tmp.number_of_nodes())
        val=0
        for j in range(1,k+1):
              for node1 in G_tmp.neighbors(node_list[randomlist[j-1]]):
                  val=val+G_tmp.edges[node_list[randomlist[j-1]],node1]['capacity']
              G_tmp.remove_node(node_list[randomlist[j-1]])
              
        #plot_graph(G)
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        #a=nx.number_connected_components(G_tmp)
        #print([len(c) for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True)])
        print("Number of connected components\n")
        print(a)
        x4.append(k)
        y4.append(a)
        z4.append(val)
    plt.plot(x0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(x0, y2, 'yo--', label="closeness_centrality")
    plt.plot(x0, y3, 'gs--', label="capacity")
    plt.plot(x0, y4, '-p',markersize=4, color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    #plt.title("Robustness of Network")
    plt.ylabel("$\Delta_r$")
    plt.xlabel("Number of nodes removed")
    plt.legend(loc="upper left")
    #plt.show()

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.savefig("plot_temporal_connected_components/test"+str1+"-"+str2+".png")
    plt.close()
    
    plt.plot(z0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(z1, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(z2, y2, 'yo--', label="closeness_centrality")
    plt.plot(z3, y3, 'gs--', label="capacity")
    plt.plot(z4, y4, '-p', markersize=4,color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    #plt.title("Robustness of Network")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_r$")
    plt.legend(loc="upper left")
    plt.savefig("plot_temporal_connected_budget/test"+str1+"-"+str2+".png")
    plt.close()
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
#    plt.show()

    

def gini_plot(G,k,*args):
    x=[]
    y={}
    count=-1
    for i in range(1,k):
       x.append(i)

    for stri in args:
        count=count+1
        y[count]=[]
       
       
        for i in range(1,k):
            gini=0
            sum_gini=sum([G[i].nodes[node1][stri] for node1 in G[i].nodes])
            if sum_gini>0:
                for node in G[i].nodes:
                    for node1 in G[i].nodes:
                        gini=gini+abs(G[i].nodes[node][stri]-G[i].nodes[node1][stri])                    

            if sum_gini>0:
                val=gini/(2*sum_gini*G[i].number_of_nodes())
                y[count].append(val)
            else:
                y[count].append(1)
			#y[count].append(gini/(2*sum_gini*G[i].number_of_nodes()))    
			
			
			#
				    
				
			
    plt.plot(x, y[0], 'r^--',markersize=8, label="BC")
    plt.plot(x, y[1],  'bs--', markersize=10, label="CC")
    plt.plot(x, y[2], 'yo--', label="Degree")
    plt.plot(x, y[3], 'gs--', label="Capacity")
    
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    plt.title("Evolution of Gini Coefficient")
    plt.ylabel("Gini Coefficient")
    plt.xlabel("Timestamp")
    
    
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):
    #        label.set_visible(False)
#    plt.show()
    plt.xticks(rotation=10)
    plt.legend(loc="upper right")
    #plt.show()

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.savefig("Gini.png")
    plt.close()
	    
    
  #  gini=0
  #  sum_gini=sum([node1[1] for node1 in bc_tmp])   
  #  for node in bc_tmp:
  #     for node1 in bc_tmp:
   #         gini=gini+abs(node[1]-node1[1])
#       gini=gini/node[1]
   # print(gini/(2*sum_gini*G[num_arg-1].number_of_nodes()))   
    


def plot_capacity(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,randomlist,str1,str2):
    
    x0=[]
    y0=[]
    z0=[]
    
    capacity_init=sum([G.nodes[x]['capacity'] for x in G.nodes])
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        percentage_dec=(capacity_init-capacity_now)/capacity_init    
        x0.append(k)
        y0.append(percentage_dec)
        z0.append(val)
     
        
        
        
    x1=[]
    y1=[]
    z1=[]
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        percentage_dec=(capacity_init-capacity_now)/capacity_init    
        x1.append(k)
        y1.append(percentage_dec)
        z1.append(val)
     
     
    x2=[]
    y2=[]
    z2=[]
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        percentage_dec=(capacity_init-capacity_now)/capacity_init    
        x2.append(k)
        y2.append(percentage_dec)
        z2.append(val)
     
    x3=[]
    y3=[]
    z3=[]
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        
        percentage_dec=(capacity_init-capacity_now)/capacity_init    
        x3.append(k)
        y3.append(percentage_dec)
        z3.append(val)
     

    x4=[]
    y4=[]
    z4=[]
    print("random\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        
        
        node_list = np.array(list(G_tmp.nodes()))
          
        k=int(i/100*G_tmp.number_of_nodes())
        val=0
        for j in range(1,k+1):
              print(j)
              for node1 in G_tmp.neighbors(node_list[randomlist[j-1]]):
                  val=val+G_tmp.edges[node_list[randomlist[j-1]],node1]['capacity']
              #val=val+G_tmp.nodes[]['capacity']
              G_tmp.remove_node(node_list[randomlist[j-1]])
              
        #plot_graph(G)
        capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        percentage_dec=(capacity_init-capacity_now)/capacity_init    
        x4.append(k)
        y4.append(percentage_dec)
        z4.append(val)
        
    plt.plot(x0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(x0, y2, 'yo--', label="closeness_centrality")
    plt.plot(x0, y3, 'gs--', label="capacity")
    plt.plot(x0, y4, '-p',markersize=4, color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    
    plt.ylabel("$\Delta_c$")
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
    plt.plot(z0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(z1, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(z2, y2, 'yo--', label="closeness_centrality")
    plt.plot(z3, y3, 'gs--', label="capacity")
    plt.plot(z4, y4, '-p', markersize=4, color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    #plt.title("Robustness of Network")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_c$")
    plt.legend(loc="upper left")
    plt.savefig("plot_capacity_connected_budget/test"+str1+"-"+str2+".png")
    plt.close()
    
    
def remove_degree_node(G,setlist,count):
    

    x=int(G.number_of_nodes()*count/100)
    
    #print(setlist)
    degree_sequence=sorted(setlist, key=lambda x:x[1],reverse=True)[:x]
   # print(degree_sequence)
    val=0
    for node in degree_sequence:
        
        #print("Node")
       # print(node[0])	
        for node1 in G.neighbors(node[0]):
                  val=val+G.edges[node[0],node1]['capacity']
        
        G.remove_node(node[0])
        """
        list_edge=sorted([G.edges[node[0],v] for v in G.neighbors(node[0])], key=lambda x:x['capacity'],reverse=False)
        #print(list_edge)
        
        
        
        
        
        count=0
        for e in list_edge:
            
            G.nodes[e['start']]['capacity']=G.nodes[e['start']]['capacity']-G.edges[(e['start'],e['end'])]['node1_deposit']
            G.nodes[e['end']]['capacity']=G.nodes[e['end']]['capacity']-G.edges[(e['start'],e['end'])]['node2_deposit']
            G.remove_edge(e['start'],e['end'])
            count=count+1
            if count==2:
                break
        """    
    return val
    
def plot_shortest(G,set_temp_deg,set_agg_deg,set_avg_deg,set_normal,randomlist,str1,str2):
    x0=[]
    y0=[]
    z0=[]
    
    
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_temp_deg,i)
        #plot_graph(G)
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                
                                
                                
                                count=count+1
                    
        n=G.number_of_nodes()
        x0.append(k)
        y0.append(1-(count/(n*(n-1))))
        z0.append(val)
        """
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y0.append(a)
        
        except Exception:
            y0.append(0)
            print("aborting") 
        
        """
        
    x1=[]
    y1=[]
    z1=[]
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_agg_deg,i)
        #plot_graph(G)
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                count=count+1
                    
        n=G.number_of_nodes()
        x1.append(k)
        y1.append(1-(count/(n*(n-1))))
        z1.append(val)
        """
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
        """    
    x2=[]
    y2=[]
    z2=[]
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_avg_deg,i)
        #plot_graph(G)
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                count=count+1
                    
        n=G.number_of_nodes()
        x2.append(k)
        y2.append(1-(count/(n*(n-1))))
        z2.append(val)
        """
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
        """    
 
    x3=[]
    y3=[]
    z3=[]
    
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        val=remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                count=count+1
                    
        n=G.number_of_nodes()
        x3.append(k)
        y3.append(1-(count/(n*(n-1))))
        z3.append(val)
        """
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())
        remove_degree_node(G_tmp,set_normal,i)
        #plot_graph(G)
        
        
        x3.append(k)
        
        try:
            c = max(nx.connected_components(G_tmp), key=len)
            G_tmp=G_tmp.subgraph(c).copy() 
            
            a=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G_tmp)
            y3.append(a)y3.append(1-(count/pow((n-1),2)))
        
        except Exception:
            y3.append(0)
            print("aborting")
        """
    x4=[]
    y4=[]
    z4=[]
    print("random\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        
        
        node_list = np.array(list(G_tmp.nodes()))
        val=0  
        k=int(i/100*G_tmp.number_of_nodes())
        for j in range(1,k+1):
              for node1 in G_tmp.neighbors(node_list[randomlist[j-1]]):
                  val=val+G_tmp.edges[node_list[randomlist[j-1]],node1]['capacity']
              G_tmp.remove_node(node_list[randomlist[j-1]])
              
        #plot_graph(G)
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                count=count+1
                    
        n=G.number_of_nodes()
        x4.append(k)
        y4.append(1-(count/(n*(n-1))))
        z4.append(val)
        
        
    
        
        
        
        
    
    plt.plot(x0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(x0, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(x0, y2, 'yo--', label="closeness_centrality")
    plt.plot(x0, y3, 'gs--', label="capacity")
    plt.plot(x0, y4, '-p', markersize=4,color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    
    plt.ylabel("$\Delta_s$")
    plt.xlabel("Number of nodes removed")
    
    
    #n = 20
    #invisible = [202,213,226,291,323,407,415,267,269] ##disturbing degrees
    #for index, label in enumerate(ax.xaxis.get_ticklabels()):
    #    if deg[index] % n != 0 and deg[index]<200 or (deg[index] in invisible):right
    #        label.set_visible(False)
#    plt.show()

    plt.legend(loc="upper left")
    #plt.show()

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.savefig("plot_shortest_connected_components/test"+str1+"-"+str2+".png")
    plt.close()
    
    
    
    plt.plot(z0, y0, 'r^--', markersize=8,label="Degree ")
    plt.plot(z1, y1,  'bs--', markersize=10, label="betweenness_centrality ")
    plt.plot(z2, y2, 'yo--', label="closeness_centrality")
    plt.plot(z3, y3, 'gs--', label="capacity")
    plt.plot(z4, y4,'-p',markersize=4, color='black', label="random")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    #plt.title("Robustness of Network")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_s$")
    plt.legend(loc="upper left")
    plt.savefig("plot_shortest_connected_budget/test"+str1+"-"+str2+".png")
    plt.close()
    
    
    
def plot_edges_shortest(G,randomlist):
    x5=[]
    y5=[]
    z5=[]
    print("cut\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        k=int(i/100*G_tmp.number_of_nodes())

        node_list = np.array(list(G_tmp.nodes()))
        listcut=[]
        print(len(randomlist))
        j=0
        val=0
        while j < 2*k:
            
            source=node_list[randomlist[j]]
            sink=node_list[randomlist[j+1]]
            
            print(source)
            print(sink)
            cut_value, partition = nx.minimum_cut(G_tmp,source, sink,capacity='capacity')

            reachable, non_reachable = partition
            cutset = set()

            for u, nbrs in ((n, G_tmp[n]) for n in reachable):

                cutset.update((u, v) for v in nbrs if v in non_reachable and (v,u) not in cutset)

            print(sorted(cutset))
            for e in cutset:
                if e not in listcut and e in G_tmp.edges:
                    listcut.append(e)
            j=j+2
            
            
        for e in listcut:
            
            if e in G_tmp.edges:
	            val=val+G_tmp.edges[e]['capacity']
	            G_tmp.remove_edge(e[0],e[1])


            
        count=0
        for c in sorted(nx.connected_components(G_tmp), key=len, reverse=True):
            S=G_tmp.subgraph(c).copy()   
            for node in S.nodes:
                for node1 in S.nodes:                      
                            if node in S.nodes and node1 in S.nodes and node!=node1 and nx.has_path(S,node,node1)==True:
                                count=count+1
                    
        n=G.number_of_nodes()    
        print(count)  
        x5.append(len(listcut))
        y5.append(1-(count/(n*(n-1))))
        z5.append(val)    
    plt.plot(x5, y5, 'r^--', markersize=8,label="Cut")
    plt.xlabel("Number of edges")
    plt.ylabel("$\Delta_s$")
    plt.legend(loc="upper left")
    plt.savefig("shortest.png")
    plt.close()
    
    plt.plot(z5, y5, 'b^--', markersize=8,label="Cut")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_s$")
    plt.legend(loc="upper left")
    plt.savefig("shortestbudget.png")
    plt.close()
 #   plt.show()

def plot_edges_capacity(G,randomlist):
    x5=[]
    y5=[]
    z5=[]
    print("cut\n")
    for i in range(1,11,1):
        G_tmp=G.copy()
        capacity_init=sum([G.nodes[x]['capacity'] for x in G.nodes])
        k=int(i/100*G_tmp.number_of_nodes())

        node_list = np.array(list(G_tmp.nodes()))
        listcut=[]
        print(len(randomlist))
        j=0
        val=0
        while j < 2*k:
            
            source=node_list[randomlist[j]]
            sink=node_list[randomlist[j+1]]
            
            print(source)
            print(sink)
            cut_value, partition = nx.minimum_cut(G_tmp,source, sink,capacity='capacity')

            reachable, non_reachable = partition
            cutset = set()

            for u, nbrs in ((n, G_tmp[n]) for n in reachable):

                cutset.update((u, v) for v in nbrs if v in non_reachable and (v,u) not in cutset)
                


            print(sorted(cutset))
            for e in cutset:
                if e not in listcut and e in G_tmp.edges:
                    listcut.append(e)
            j=j+2
            
            
        for e in listcut:
            if e in G_tmp.edges:
	            val=val+G_tmp.edges[e]['capacity']
	            G_tmp.remove_edge(e[0],e[1])
            	           	       

            
            
   #     capacity_now=sum([G_tmp.edges[e]['capacity'] for e in G_tmp.edges])
        percentage_dec=val/capacity_init    
        
        x5.append(len(listcut))
        y5.append(percentage_dec)
        z5.append(val)    
    plt.plot(x5, y5, 'r^--', markersize=8,label="Cut")
    plt.xlabel("Number of edges")
    plt.ylabel("$\Delta_c$")
    plt.legend(loc="upper left")
    plt.savefig("capacity.png")
    plt.close()
    #plt.show()
    plt.plot(z5, y5, 'b^--', markersize=8,label="Cut")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_c$")
    plt.legend(loc="upper left")
    plt.savefig("capacitybudget.png")
    plt.close()
    


def plot_edges_reach(G,randomlist):
    x5=[]
    y5=[]
    z5=[]
    print("cut\n")
    print(randomlist)
    for i in range(1,11,1):
        G_tmp=G.copy()
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        r=len(largest_cc)
        k=int(i/100*G_tmp.number_of_nodes())
        
        node_list = np.array(list(G_tmp.nodes()))
        listcut=[]
        print(len(randomlist))
        j=0
        val=0
        while j < 2*k:
            
            source=node_list[randomlist[j]]
            sink=node_list[randomlist[j+1]]
            
            print(source)
            print(sink)
            cut_value, partition = nx.minimum_cut(G_tmp,source, sink,capacity='capacity')

            reachable, non_reachable = partition
            cutset = set()
            print("reach")
            print(reachable)
            print("unreach")
            print(non_reachable)		
            for u, nbrs in ((n, G_tmp[n]) for n in reachable):

                cutset.update((u, v) for v in nbrs if v in non_reachable and (v,u) not in cutset)
                



            print(sorted(cutset))
            for e in cutset:
                print(e in G_tmp.edges)
                if e not in listcut and e in G_tmp.edges:
                    listcut.append(e)
            j=j+2
            
            
        for e in listcut:
            print(e)
            print(e in G_tmp.edges)
            if e in G_tmp.edges:
	            val=val+G_tmp.edges[e]['capacity']
	            G_tmp.remove_edge(e[0],e[1])


            
        largest_cc = max(nx.connected_components(G_tmp), key=len)
        
        #a=nx.number_connected_components(G_tmp)
        a=(r-len(largest_cc))/r
        
        
        x5.append(len(listcut))
        print("listcut")
        print(len(listcut))

        y5.append(a)
        z5.append(val)    
    plt.plot(x5, y5, 'r^--', markersize=8,label="Cut")
    plt.xlabel("Number of edges")
    plt.ylabel("$\Delta_r$")
    plt.legend(loc="upper left")
    #plt.savefig("plot_shortest_connected_budget/test"+str1+"-"+str2+".png")
    #plt.close()
    plt.savefig("reach.png")
    plt.close()
    plt.plot(z5, y5, 'b^--', markersize=8,label="Cut")
    plt.xlabel("Adversarial Budget (Raiden Energy Token)")
    plt.ylabel("$\Delta_r$")
    plt.legend(loc="upper left")
    plt.savefig("reachbudget.png")
    plt.close()
    #plt.show()

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

def goodnessOfFit(deg, cnt):
    #calculated p-value: 0.8172
    (Kmin, gamma, nmin, baseDist) = fittingGamma(deg, cnt)
    ntotal = np.sum(cnt)

    powerLawProbs = generatePowerLawDist(gamma, Kmin, np.amax(deg))
    counter = 0
    # we chose 2500 iterations to have 2 digits precision for the p-value
    for i in range(2500):
        syntDeg, syntCnt = generateSyntheticData(ntotal, nmin, gamma, Kmin, np.amax(deg), powerLawProbs)
        optK, approxGamma, nmin, minDist=fittingGamma(syntDeg,syntCnt)
        print(i, approxGamma, minDist)
        if(baseDist < minDist):
            counter+=1
    print("p-value:", counter/2500)
    return syntDeg,syntCnt

def generatePowerLawDist(gamma, Kmin, maxDeg):
    probI = []
    for q in range(Kmin, 2000, 1):
        probI.append(pow(q, -gamma) / spcl.zeta(gamma, Kmin))
    probI.append(1-np.sum(probI))
    return probI

def gammaF(deg, cnt, Kmin):
    sum = 0
    counter = 0
    index = 0
    for x in deg:
        if(x>=Kmin):
            print(Kmin)
            sum += cnt[counter]*np.log(x / (Kmin - 0.5))
            index+=cnt[counter]
        counter+=1
    gamma = 1 + index * (1 / sum)
    return (gamma, index)

#degreeCentrality
#2.1387317708757214
def fittingGamma(deg, cnt):
    minmaxD = 100
    optimalK = 0
    bestgamma = 0
    nmin = 0
    maxD = []
    for Kmin in deg:
        gamma, ind = gammaF(deg, cnt, Kmin)

        cdf = [] #cumulative distribution func
        cdfOrdered = {}
        for y in deg:
            if(spcl.zeta(gamma,Kmin)!=0.0):
                CDF = 1-(spcl.zeta(gamma,y))/(spcl.zeta(gamma,Kmin))
            else:
                CDF = 0
            cdf.append(CDF)
            cdfOrdered[y] = CDF

        #Kolmogorov-Smirnov test
        maxDistance = 0
        for z in deg:
            if (z >= Kmin):
                if abs(cdfOrdered[z]-empiricalCDF(deg,cnt, Kmin,z))>maxDistance:
                    maxDistance = abs(cdfOrdered[z]-empiricalCDF(deg,cnt,Kmin,z))
        if maxDistance < minmaxD:
            minmaxD = maxDistance
            optimalK = Kmin
            nmin = ind
            bestgamma = gamma
        maxD.append(maxDistance)

    #deviation=1/math.sqrt(np.sum(cnt)*((mp.zeta(gamma,optimalK,2)/spcl.zeta(gamma,optimalK))-math.pow((mp.zeta(gamma,optimalK,1)/spcl.zeta(gamma,optimalK)),2)))

    #print(optimalK)
    #print("Approximated exponent in the power-law distribution: ",gammaF(deg, cnt, optimalK))
    #print("Deviation of the approximation: ", deviation)
    #print("Max distance between CDF and empirical data: ", minmaxD)
    #plt.plot(np.array(deg), np.array(maxD), '.', color='green')
    #plt.title("Kolmogorov-Smirnov test")
    #plt.ylabel("D")
    #plt.xlabel("K")
    #plt.yscale("log")

    return (optimalK, bestgamma, nmin, minmaxD)
def empiricalCDF(deg, cnt, Kmin, z):
    sum = 0
    counter = 0
    total = 0
    for w in deg:
        if (w>=Kmin):
            total+=cnt[counter]
            if (w<=z):
                sum+=cnt[counter]
        counter+=1
    return sum/total

def generateSyntheticData(ntotal, nmin, gamma, Kmin, Kmax, powerLawProbs):
    syntheticDegrees = []
    for x in range(ntotal):
        rand=random.randint(1,ntotal)
        if rand < ntotal-nmin:
            syntheticDegrees.append(random.randint(1, Kmin - 1))
        else:
            syntheticDegrees.append(np.random.choice(np.arange(Kmin, 2000+1), p=powerLawProbs))

    syntheticDegreeCount = collections.Counter(syntheticDegrees)
    syntDeg, syntCnt = zip(*syntheticDegreeCount.items())     
    return syntDeg,syntCnt

def li_smax_graph(degree_seq):
    """Generates a graph based with a given degree sequence and maximizing
    the s-metric.  Experimental implementation.
    Maximum s-metrix  means that high degree nodes are connected to high
    degree nodes.
    - `degree_seq`: degree sequence, a list of integers with each entry
       corresponding to the degree of a node.
       A non-graphical degree sequence raises an Exception.
    Reference::
      @unpublished{li-2005,
       author = {Lun Li and David Alderson and Reiko Tanaka
                and John C. Doyle and Walter Willinger},
       title = {Towards a Theory of Scale-Free Graphs:
               Definition, Properties, and  Implications (Extended Version)},
       url = {http://arxiv.org/abs/cond-mat/0501169},
       year = {2005}
      }
    The algorithm::
     STEP 0 - Initialization
     A = {0}
     B = {1, 2, 3, ..., n}
     O = {(i; j), ..., (k, l),...} where i < j, i <= k < l and
             d_i * d_j >= d_k *d_l
     wA = d_1
     dB = sum(degrees)
     STEP 1 - Link selection
     (a) If |O| = 0 TERMINATE. Return graph A.
     (b) Select element(s) (i, j) in O having the largest d_i * d_j , if for
             any i or j either w_i = 0 or w_j = 0 delete (i, j) from O
     (c) If there are no elements selected go to (a).
     (d) Select the link (i, j) having the largest value w_i (where for each
             (i, j) w_i is the smaller of w_i and w_j ), and proceed to STEP 2.
     STEP 2 - Link addition
     Type 1: i in A and j in B.
             Add j to the graph A and remove it from the set B add a link
             (i, j) to the graph A. Update variables:
             wA = wA + d_j -2 and dB = dB - d_j
             Decrement w_i and w_j with one. Delete (i, j) from O
     Type 2: i and j in A.
         Check Tree Condition: If dB = 2 * |B| - wA.
             Delete (i, j) from O, continue to STEP 3
         Check Disconnected Cluster Condition: If wA = 2.
             Delete (i, j) from O, continue to STEP 3
         Add the link (i, j) to the graph A
         Decrement w_i and w_j with one, and wA = wA -2
     STEP 3
         Go to STEP 1
    The article states that the algorithm will result in a maximal s-metric.
    This implementation can not guarantee such maximality. I may have
    misunderstood the algorithm, but I can not see how it can be anything
    but a heuristic. Please contact me at sundsdal@gmail.com if you can
    provide python code that can guarantee maximality.
    Several optimizations are included in this code and it may be hard to read.
    Commented code to come.
    """

    if not is_valid_degree_sequence(degree_seq):
        raise nx.NetworkXError('Invalid degree sequence')
    degree_seq.sort()  # make sure it's sorted
    degree_seq.reverse()
    degrees_left = degree_seq[:]
    A_graph = nx.Graph()
    A_graph.add_node(0)
    a_list = [False] * len(degree_seq)
    b_set = set(range(1, len(degree_seq)))
    a_open = set([0])
    O = []
    for j in b_set:
        heapq.heappush(O, (-degree_seq[0] * degree_seq[j], (0, j)))
    wa = degrees_left[0]  # stubs in a_graph
    db = sum(degree_seq) - degree_seq[0]  # stubs in b-graph
    a_list[0] = True  # node 0 is now in a_Graph
    bsize = len(degree_seq) - 1  # size of b_graph
    selected = []
    weight = 0
    while O or selected:
        if len(selected) < 1:
            firstrun = True
            while O:
                (newweight, (i, j)) = heapq.heappop(O)
                if degrees_left[i] < 1 or degrees_left[j] < 1:
                    continue
                if firstrun:
                    firstrun = False
                    weight = newweight
                if not newweight == weight:
                    break
                heapq.heappush(selected, [-degrees_left[i], \
                                          -degrees_left[j], (i, j)])
            if not weight == newweight:
                heapq.heappush(O, (newweight, (i, j)))
            weight *= -1
        if len(selected) < 1:
            break

        [w1, w2, (i, j)] = heapq.heappop(selected)
        if degrees_left[i] < 1 or degrees_left[j] < 1:
            continue
        if a_list[i] and j in b_set:
            # TYPE1
            a_list[j] = True
            b_set.remove(j)
            A_graph.add_node(j)
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            wa += degree_seq[j] - 2
            db -= degree_seq[j]
            bsize -= 1
            newweight = weight
            if not degrees_left[j] == 0:
                a_open.add(j)
                for k in b_set:
                    if A_graph.has_edge(j, k): continue
                    w = degree_seq[j] * degree_seq[k]
                    if w > newweight:
                        newweight = w
                    if weight == w and not newweight > weight:
                        heapq.heappush(selected, [-degrees_left[j], \
                                                  -degrees_left[k], (j, k)])
                    else:
                        heapq.heappush(O, (-w, (j, k)))
                if not weight == newweight:
                    while selected:
                        [w1, w2, (i, j)] = heapq.heappop(selected)
                        if degrees_left[i] * degrees_left[j] > 0:
                            heapq.heappush(O, [-degree_seq[i] * degree_seq[j], (i, j)])
            if degrees_left[i] == 0:
                a_open.discard(i)

        else:
            # TYPE2
            if db == (2 * bsize - wa):
                # tree condition
                # print "removing because tree condition    "
                continue
            elif db < 2 * bsize - wa:
                raise nx.NetworkXError("THIS SHOULD NOT HAPPEN! - not graphable")
                continue
            elif wa == 2 and bsize > 0:
                # print "removing because disconnected  cluster"
                # disconnected cluster condition
                continue
            elif wa == db - (bsize) * (bsize - 1):
                # print "MYOWN removing because disconnected  cluster"
                continue
            A_graph.add_edge(i, j)
            degrees_left[i] -= 1
            degrees_left[j] -= 1
            if degrees_left[i] < 1:
                a_open.discard(i)
            if degrees_left[j] < 1:
                a_open.discard(j)
            wa -= 2
            if not degrees_left[i] < 0 and not degrees_left[j] < 0:
                selected2 = (selected)
                selected = []
                while selected2:
                    [w1, w1, (i, j)] = heapq.heappop(selected2)
                    if degrees_left[i] * degrees_left[j] > 0:
                        heapq.heappush(selected, [-degrees_left[i], \
                                                  -degrees_left[j], (i, j)])
    return A_graph
       
def main():
    
    G_main=nx.Graph()
    num_arg=len(sys.argv)-1
    source={}
    set_nodes=[]
    G={}
    for i in range(1,num_arg):
    
        with open(sys.argv[i]) as f:
            print(sys.argv[i])
            source[i] = json.load(f)
        
            G[i]=readGraph(source[i])
            G[i].remove_nodes_from(list(nx.isolates(G[i])))
            set_graph_color(G[i])
            #plot_graph1(G[i])
            
            
            
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
                
    capacity_node_main(G_main,G,num_arg)
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
    
    bc_edge_tmp=bc_edge_temporal(G,G_main,num_arg-1,shortest_path_matrix)
    
    set_temp_deg=temporal_degree(G,num_arg,G_main)
    set_agg_deg=aggregated_degree(G_main,num_arg)
    cap_agg=aggregated_capacity(G_main)
    cap_temp=temporal_capacity(G,G_main,num_arg)
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
    cap_node=[]    
    for v in G[num_arg-1].nodes:
        closeness_centrality_node.append((v,G[num_arg-1].nodes[v]['closeness_centrality']))
    for v in G[num_arg-1].nodes:
        cap_node.append((v,G[num_arg-1].nodes[v]['capacity']))    
    
    #plot_graph1(G[num_arg-1])
    degree_node=[]
    
    for v in G[num_arg-1].nodes:
        degree_node.append((v,G[num_arg-1].nodes[v]['degree']))
    
    
   
    
    
#    cap_cen=sorted(G[num_arg-1].nodes(data=True), key=lambda x:x[1]['capacity'],reverse=True)
    cap_cen=sorted(G[num_arg-1].nodes(data=True), key=lambda x:x[1]['capacity'],reverse=True)
    
    for node in cap_cen:
        print(node)

    cap_node=[]
    for v in G[num_arg-1].nodes:
        cap_node.append((v,G[num_arg-1].nodes[v]['capacity']))
           
     #  cap_node=cap_node/(num_arg-1)
      # G[num_arg-1].nodes[node]['capacity']=cap_node
       
       #print(str(node)+" "+str(cap_node))
    
    cap_cen=sorted(G[num_arg-1].nodes(data=True), key=lambda x:x[1]['capacity'],reverse=True)
    for i in range(1,num_arg):    
	    measure_of_centrality1 = nx.closeness_centrality(G[i],u=None, distance=None, wf_improved=True)
	    nx.set_node_attributes(G[i], measure_of_centrality1, 'closeness_centrality')
    gini_plot(G,num_arg,'betweenness_centrality','closeness_centrality','degree','capacity')
    
        
    
    
   # goodnessOfFit(deg,cnt)
    
    for i in range(1,num_arg):
	    print("Degree assortativity",nx.algorithms.assortativity.degree_assortativity_coefficient(G[i]))

	   
    degree_sequence = sorted([d for n, d in G[num_arg-1].degree()], reverse=True)  # degree sequence
    G_max = li_smax_graph(list(degree_sequence))
    print("smetric",float(sum([G[num_arg-1].degree(u) * G[num_arg-1].degree(v) for (u, v) in G[num_arg-1].edges()]))/float(sum([G_max.degree(u) * G_max.degree(v) for (u, v) in G_max.edges()])))

        
    #Generate 5 random numbers between 10 and 30
    randomlist = []
    randomlist=random.sample(range(G[num_arg-1].number_of_nodes()), int(20*G[num_arg-1].number_of_nodes()/100))
    """
    for i in range(0,30):
        n = random.randint(1,G[num_arg-1].number_of_nodes())
        randomlist.append(n)
    #randomlist = random.sample(range(1,G[num_arg-1].number_of_nodes() ), 30)
    """
    print(randomlist)
    #plotDeg(bc_tmp,bc_ag,set_bc)
    plot_graph(G[num_arg-1],set_temp_deg,bc_tmp,cc_temporal,cap_temp,randomlist,"Temporal",sys.argv[num_arg])
    plot_graph(G[num_arg-1],set_agg_deg,bc_ag,cc_agg,cap_agg,randomlist,"Aggregated",sys.argv[num_arg])
    plot_graph(G[num_arg-1],set_temp_deg,set_bc,cc_average,cap_temp,randomlist,"Average",sys.argv[num_arg])
    plot_graph(G[num_arg-1],degree_node,betweenness_centrality_node,closeness_centrality_node,cap_node,randomlist,"Normal",sys.argv[num_arg])
    
    plot_shortest(G[num_arg-1],set_temp_deg,bc_tmp,cc_temporal,cap_temp,randomlist,"Temporal",sys.argv[num_arg])
    plot_shortest(G[num_arg-1],set_agg_deg,bc_ag,cc_agg,cap_agg,randomlist,"Aggregated",sys.argv[num_arg])
    plot_shortest(G[num_arg-1],set_temp_deg,set_bc,cc_average,cap_temp,randomlist,"Average",sys.argv[num_arg])
    plot_shortest(G[num_arg-1],degree_node,betweenness_centrality_node,closeness_centrality_node,cap_node,randomlist,"Normal",sys.argv[num_arg])
    
   
    plot_capacity(G[num_arg-1],set_temp_deg,bc_tmp,cc_temporal,cap_temp,randomlist,"Temporal",sys.argv[num_arg])
    
    plot_capacity(G[num_arg-1],set_agg_deg,bc_ag,cc_agg,cap_agg,randomlist,"Aggregated",sys.argv[num_arg])
    plot_capacity(G[num_arg-1],set_temp_deg,set_bc,cc_average,cap_temp,randomlist,"Average",sys.argv[num_arg])
    plot_capacity(G[num_arg-1],degree_node,betweenness_centrality_node,closeness_centrality_node,cap_node,randomlist,"Normal",sys.argv[num_arg])
    randomlist=random.sample(range(G[num_arg-1].number_of_nodes()), 2*int(10/100*G[num_arg-1].number_of_nodes()))
    plot_edges_reach(G[num_arg-1],randomlist)

    plot_edges_shortest(G[num_arg-1],randomlist)
    plot_edges_capacity(G[num_arg-1],randomlist)
    
      # degree sequence
    Distribution(sorted([d[1] for d in set_temp_deg], reverse=True),sorted([d[1] for d in set_agg_deg], reverse=True),sorted([d[1] for d in set_temp_deg], reverse=True),"Degree Distribution",sys.argv[num_arg])
    Distribution(sorted([d[1] for d in bc_tmp], reverse=True),sorted([d[1] for d in bc_ag], reverse=True),sorted([d[1] for d in set_bc], reverse=True),"Betweeness Centrality Distribution",sys.argv[num_arg])
    Distribution(sorted([d[1] for d in cc_temporal], reverse=True),sorted([d[1] for d in cc_agg], reverse=True),sorted([d[1] for d in cc_average], reverse=True),"Closeness Centrality Distribution",sys.argv[num_arg])
    
    deg={}
    cnt={}
    
    print(G[num_arg-1].number_of_nodes())
    
   # num=len([d for n,d in G[num_arg-1].degree() if d>0])
    #alpha=1+(num/sum([np.log(d/(mindeg-0.5)) for n,d in G[num_arg-1].degree() if d>0]))
    #print(alpha)
    alpha=0    
    #print(alpha)    
    
if __name__ == '__main__':
    main()
    
    
