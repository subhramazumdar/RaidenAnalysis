import networkx as nx
import json
from pprint import pprint
import math
from matplotlib import mlab
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import random
import types
import operator
from matplotlib.ticker import MaxNLocator
import collections
import numpy as np
from mpmath import mp
from decimal import *
import string
from random import randint
from random import seed
from scipy import special as spcl
import heapq
from networkx.algorithms import approximation


allMoneyStaked = 0
allWeights = []

def defineGraph(data) -> object:
    global allMoneyStaked
    G = nx.Graph()
    for x in range(len(data)):
        G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], capacity=data[x]['capacity'])
        allWeights.append(int(data[x]['capacity']))
        allMoneyStaked+=int(data[x]['capacity'])
    return G


##https://graph.lndexplorer.com/api/graph
def readFile() -> object:
    with open('graph20190418.json') as f:
        data = json.load(f)
    return data['edges']

def main():
    G = defineGraph(readFile())
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    #basicStatistics(G)
    degreeDistribution(G) #percentile graph needs to be added
    weightsDistribution(G)
    shortestPaths(G)
    drawGraph(G)
    vertexConnectivity(G)
    edgeConnectivity(G)
    randomVertexConnectivity(G)
    centrality(G)
    clusteringCoefficient(G)
    simpleStatistics(G)
    fittingGamma(deg, cnt)
    goodnessOfFit(deg, cnt)
    percolationThresholdPrediction(deg)
    attackingBetweenness(G,deg, cnt)
    attackingHighDegrees(G)
    drawHistogram(G)
    improvingRobustness(G)
    edgeCapacityDistributions(G)

def edgeCapacityDistributions(G):
    richestStakers = {}
    for n in G.nodes():
        capacityOfn = 0
        for m in G.neighbors(n):
            capacityOfn += int(G[n][m]['capacity'])
        richestStakers[n]=capacityOfn
    sortedStakers = list(sorted(richestStakers.items(), key=operator.itemgetter(1), reverse=True))
    sortedCapacities = []
    for i in sortedStakers:
        sortedCapacities.append(i[1])

    lostCapacityPercentage = []
    lostCapacities = []
    lostCapacity = 0
    counter = 0
    for k in sortedStakers:
        counter+=1
        lostCapacity+=k[1]
        lostCapacities.append(k[1])
        lostCapacityPercentage.append(lostCapacity/allMoneyStaked)
        print(counter, lostCapacity/allMoneyStaked)
        for l in G.neighbors(k[0]):
            index = 0
            previousCapacity= 0
            for m in range(len(sortedStakers)):
                if sortedStakers[m][0]==l:
                    previousCapacity = sortedStakers[m][1]
                    index = m
            del sortedStakers[index]
            sortedStakers.insert(index,(l,previousCapacity-int(G[k[0]][l]['capacity'])))
        sortedStakers.remove(k)
        G.remove_node(k[0])
        sorted(sortedStakers, key=operator.itemgetter(1), reverse=True)
        if counter == 50:
            break

    fig, ax1 = plt.subplots()
    t = np.arange(0, 50, 1)
    lns1 = ax1.plot(t, lostCapacityPercentage, 'b-', label='Capacity Loss %')

    ax1.set_xlabel('Number of removed nodes')
    # Make the y-axis label, ticks and tick labels match the line color.
    # ax1.set_ylabel('MKR Price', color='b')
    ax1.set_ylabel('Ratio of lost capacity and balk capacity', color='b')
    ax1.tick_params('y', colors='b')

    ax2 = ax1.twinx()
    lns3 = ax2.plot(t, lostCapacities, 'r-', label='Capacity loss in BTC')
    ax2.set_ylabel('Capacity loss by removing a node in 10 BTC', color='r')
    ax2.tick_params('y', colors='r')

    # added these three lines
    lns = lns3 + lns1  # lns2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=5)

    # Put a nicer background color on the legend.
    # legend.get_frame().set_facecolor('C0')

    fig.tight_layout()
    plt.show()


def improvingRobustness(G):
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    originalRobustness = 0
    #G.remove_nodes_from(list(nx.connected_components(G))[1])
    percolationThreshold = 0
    for x in range(1000):
        highestDegrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
        G.remove_node(next(iter(highestDegrees))[0])
        percolationThreshold += 1
        largest_cc = max(nx.connected_components(G), key=len)
        originalRobustness+=(len(largest_cc)/(originalGiantComponentSize*originalGiantComponentSize))
        print(x, originalRobustness)
        if (len(largest_cc) < originalGiantComponentSize / 100):
            print("Percolation Threshold: ", percolationThreshold, percolationThreshold / originalGiantComponentSize)
            break


def drawHistogram(G): #attack effect on avg shortest path length
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    a=[2.806222412074612, 2.871217085442275, 2.91560045459128, 2.9437603442568685, 2.9740025229103972, 3.0045886322304565, 3.037044687225701, 3.0553851019952716, 3.075570909398242, 3.0890175949054375, 3.1157498785629887, 3.131813039281431, 3.146865049134942, 3.163353070623199, 3.1801390239729423, 3.1856303960443033, 3.2041163947454585, 3.2125176424347455, 3.226031083245262, 3.2354480642086663, 3.2458709321869597, 3.2528784561124127, 3.2657448307356463, 3.275644221485353, 3.2882485866327467, 3.3025002028102657, 3.3139604454666687, 3.325104794319505, 3.3376269239216745, 3.3465795761157207, 3.354999742940903]
    b = list(range(0, 31))
    c=[2.806222412074612, 2.806840895554021, 2.8062060695075983, 2.806240810145709, 2.8054799649474504, 2.8050401865595025, 2.805183957841351, 2.80548635115336, 2.8071260429750997, 2.8071066774929436, 2.806607488708337, 2.806631682697912, 2.8063720821020826, 2.8065669709319105, 2.8056500936474484, 2.806809639127011, 2.806961234777178, 2.806665944865557, 2.806500923697486, 2.8066547251933036, 2.8064298690301666, 2.8060923181460953, 2.8060408625703954, 2.8057784569538504, 2.8057073212363317, 2.8056447519017835, 2.805837125196155, 2.805820866731267, 2.8060240140210753, 2.8058008903104583, 2.8054164087648368]
    plt.plot(b, a, '.', color='red', markersize=12, label="HDR")
    plt.plot(b,c,'.',color='blue',markersize=12, label="Random")
    plt.title("Attack effect on path length")
    plt.ylabel("Average shortest path length")
    plt.xlabel("Number of removed nodes")

    plt.legend(loc="upper left")

    plt.show()

#Percolation Threshold: 381 0.1627509611277232
#Percolation Threshold: 2258 0.964545066211021
#Percolation Threshold: 330 0.14096539940196498
def attackingHighDegrees(G):
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    G.remove_nodes_from(list(nx.connected_components(G))[1])
    percolationThreshold = 0
    PP = []  # probability that a random node belongs to the giant component
    p = []  # percolation threshold
    diameter = []
    for x in range(30):
        highestDegrees = sorted(G.degree, key=lambda x: x[1], reverse=True)
        G.remove_node(next(iter(highestDegrees))[0])
        percolationThreshold+=1
        largest_cc = max(nx.connected_components(G), key=len)
        PP.append(len(largest_cc) / originalGiantComponentSize)
        p.append(percolationThreshold)
        shortpahts=nx.algorithms.shortest_paths.generic.average_shortest_path_length(G.subgraph(largest_cc))
        diameter.append(shortpahts)
        print(percolationThreshold, len(largest_cc), shortpahts)
        if (len(largest_cc) < originalGiantComponentSize / 100):
            print("Percolation Threshold: ", percolationThreshold, percolationThreshold/originalGiantComponentSize)
            break
    print(diameter)

    return (p,PP, diameter)



#Percolation Threshold:  561 (betweenness centrality)
#Percolation Threshold: 1171 (closeness centrality)
#Percolation Threshold:  144 0.061512174284493806 (betweenness)
def attackingBetweenness(G, deg, cnt):
    gor = G.copy()
    GorHighDegreeAttack = G.copy()
    (pHigh, PPHigh) = attackingHighDegrees(GorHighDegreeAttack)
    (pRandom, PPrandom) = randomVertexConnectivity(gor)
    #betweennessC=nx.algorithms.centrality.closeness_centrality(G)
    betweennessC = nx.algorithms.centrality.betweenness_centrality(G, k=300, normalized=True) #higher k yields better approx of betweenness

    PP=[] #probability that a random node belongs to the giant component
    p=[] #percolation threshold
    sortedB = sorted(betweennessC.items(), key=operator.itemgetter(1), reverse=True)
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    percolationThreshold = 0
    for x in range(1000):
        #betweennessC = nx.algorithms.centrality.closeness_centrality(G)
        betweennessC = nx.algorithms.centrality.betweenness_centrality(G, k=50,normalized=True)  # higher k yields better approx of betweenness
        sortedB = sorted(betweennessC.items(), key=operator.itemgetter(1), reverse=True)
        try:
            G.remove_node(sortedB[0][0])
        except nx.exception.NetworkXError:
            continue
        G.remove_nodes_from(list(nx.isolates(G)))
        percolationThreshold+=1
        largest_cc = max(nx.connected_components(G), key=len)
        PP.append(len(largest_cc)/originalGiantComponentSize)
        p.append(percolationThreshold)
        if (len(largest_cc)/originalGiantComponentSize < 0.01):
            print("Percolation Threshold: ", percolationThreshold, percolationThreshold/originalGiantComponentSize)
            break
    for y in range(len(PPrandom)-len(PP)):
        PP.append(0)
    for z in range(len(PPrandom)-len(PPHigh)):
        PPHigh.append(0)
    plt.plot(pRandom, PP, '.', color='red', markersize=3, label="Attacking betweenness")
    plt.plot(pRandom, PPrandom, '.', color='blue', markersize=3, label="Random failures")
    plt.plot(pRandom, PPHigh, '.', color='green', markersize=3, label="Attacking high-degree")
    #plt.plot(p, PP, '.', color='red')
    #plt.plot(pRandom, PPrandom, '.', color='green')
    plt.title("Robustness of Lightning Network")
    plt.ylabel("Percentage of the giant component")
    plt.xlabel("Percentage of nodes")


    plt.legend(loc="upper right")

    #plt.legend(['Attack', 'Random failures'], loc='upper right')

    plt.show()


def percolationThresholdPrediction(deg):
    gamma = 2.1387317708757214
    kmin=np.amin(deg)
    kmax=np.amax(deg)
    fc = 1- 1/(((gamma-2)*pow(kmin,gamma-2)*pow(kmax,3-gamma)/(3-gamma))-1)
    print(kmin, kmax, fc)
    print("Estimated network size: ",pow((kmax/kmin),gamma-1))
    print("Avalanche exponent: ",gamma/(gamma-1))
    print("First moment: ",moments(gamma,kmin,1))
    print("Fc error: ",1-1/moments(gamma,kmin,1))

    #K=(2-gamma)*(pow(kmax,3-gamma)-pow(kmin,3-gamma))/((3-gamma)*(pow(kmax,2-gamma)-pow(kmin,2-gamma)))
    K = abs((2-gamma)/(3 - gamma)) * (pow(kmax, 3 - gamma)* pow(kmin, gamma - 2))
    Fc=1-1/(K-1)
    print("Fc for scale-free: ", Fc)

def moments(gamma, xmin, m):
    return (gamma-1)*pow(xmin,m)/(gamma-1-m)

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

    return syntDeg, syntCnt

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

def simpleStatistics(G):
    a=[2234, 2117, 1493, 2196, 2160, 1767, 1908, 1796, 685, 650, 1934, 2234, 1362, 1872, 1731, 493, 459, 1060, 486, 471, 1546, 1439, 861, 1766, 1604, 910, 2229, 878, 938, 819]
    print(np.mean(a)/G.order())

def basicStatistics(G):
    print("Number of LN nodes : ", G.order())
    print("Number of LN payment channels: ", G.size())
    print("Density of LN: ",nx.classes.function.density(G))
    print("Amount of BTC, denominated in Satoshis, in all payment channels:", allMoneyStaked)
    print("Number of connected components: ", len(list(nx.connected_components(G))))
    print("Maximal independent set: ", len(nx.algorithms.maximal_independent_set(G)))
    print("Number of bridges: ", len(list(nx.algorithms.bridges(G))))
    print("Size of the dominating set: ", len(nx.algorithms.dominating_set(G)))
    print("LN is Chordal graph: ", nx.algorithms.chordal.is_chordal(G))
    print("LN degree assortativity",nx.algorithms.assortativity.degree_assortativity_coefficient(G))
    #print("LN rich-club coefficient: ", nx.algorithms.richclub.rich_club_coefficient(G))
    #print("LN rich-club normalized coefficient: ", nx.algorithms.richclub.rich_club_coefficient(G, normalized=True))
    G.remove_nodes_from(list(nx.connected_components(G))[1]) #there is a small second component
    G.remove_nodes_from(list(nx.connected_components(G))[2])  # there is a small second component
    G.remove_nodes_from(list(nx.connected_components(G))[1])  # there is a small second component
    #print("LN diameter: ", nx.algorithms.distance_measures.diameter(G)) #6
    #print("LN radius", nx.algorithms.distance_measures.radius(G)) #3
    #print("LN Wiener index", nx.algorithms.wiener_index(G)) #7686159.0
    #print("LN is Eulerian: ",nx.algorithms.is_eulerian(G))
    #print("LN is planar: ", nx.algorithms.planarity.check_planarity(G))
    #print("Number of isolates in LN: ", list(nx.isolates(G)))
    #print("LN's S-metric: ", smetric(G)) #0.6879664061934981
    print("LN average clustering coefficient", approximation.clustering_coefficient.average_clustering(G))
    print("LN's transitivity: ", nx.algorithms.cluster.transitivity(G))
    print("Average shortest paths: ",nx.algorithms.shortest_paths.generic.average_shortest_path_length(G)) # 2.806222412074612
    print("LN's largest clique size: ", nx.algorithms.approximation.clique.max_clique(G))
    print("Adjacency spectrum of LN: ", nx.linalg.spectrum.adjacency_spectrum(G))
    #only 81 onion node :(

def clusteringCoefficient(G):
    sortedNodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    print(sortedNodes[0][0])
    highDegreeNodes = []
    for x in range(30):
        highDegreeNodes.append(sortedNodes[x][0])

    clusteringCoefficients=nx.algorithms.cluster.clustering(G)
    valueArr = []
    for key,value in clusteringCoefficients.items():
        valueArr.append(value)
        print(value)
    binwidth = 0.01
    plt.hist(valueArr, bins=np.arange(min(valueArr), max(valueArr) + binwidth, binwidth))

    plt.title("Clustering-coefficient histogram")
    plt.ylabel("Number of nodes")
    plt.xlabel("Clustering-coefficient")

    left, right = plt.xlim()  # return the current xlim
    plt.xlim(-0,1, 1.1)
    print(left, right)

    plt.yscale("log")

    plt.show()

def centrality(G):
    G.remove_nodes_from(list(nx.connected_components(G))[1])  # there is a small second component
    btwCentrality = nx.algorithms.centrality.closeness_centrality(G) #k=300, normalized=True, seed=123
    btwList=btwCentrality.items()
    valueArr = []
    for x in btwList:
        valueArr.append(x[1])
    binwidth=0.002
    plt.hist(valueArr, bins=np.arange(min(valueArr), max(valueArr) + binwidth, binwidth))

    plt.title("Closeness-centrality histogram")
    plt.ylabel("Number of nodes")
    plt.xlabel("Closeness-centrality")

    plt.yscale("log")

    plt.show()

#calculating critical percolation threshold for random failures and attacks
def randomVertexConnectivity(G):
    originalGiantComponentSize = len(list(nx.connected_components(G))[0])
    percolationThreshold = 0
    G.remove_nodes_from(list(nx.connected_components(G))[1])  # there is a small second component
    PP = []  # probability that a random node belongs to the giant component
    p = []  # percolation threshold
    diameter = []
    for x in range(30):
        sortedNodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)
        toremove=randint(0,len(sortedNodes)-1)
        G.remove_node(next(iter(sortedNodes[toremove])))
        del sortedNodes[toremove]
        largest_cc = max(nx.connected_components(G), key=len)
        shortpahts = nx.algorithms.shortest_paths.generic.average_shortest_path_length(G.subgraph(largest_cc))
        diameter.append(shortpahts)
        percolationThreshold += 1
        print(percolationThreshold, len(largest_cc), shortpahts)
        p.append(percolationThreshold/originalGiantComponentSize)
        PP.append(len(largest_cc)/originalGiantComponentSize)
        if(len(largest_cc)/originalGiantComponentSize<0.01):
            print("Percolation Threshold: ", percolationThreshold, percolationThreshold / originalGiantComponentSize)
            break

    print(diameter)
    #plt.plot(p, PP, '.', color='red')
    #plt.title("Robustness of Lightning Network")
    #plt.ylabel("Size of giant component")
    #plt.xlabel("Number of removed nodes")

    return (p,PP)

def vertexConnectivity(G):
    sortedNodes = sorted(G.degree(), key=lambda x: x[1], reverse=True)
    comp=[0]
    cnt=[2]
    Gor=[]
    for x in range(30):
        Gor.append(G.copy())
    print("what?",len(list(nx.connected_components(G))))
    for x in range(200):
        #print(sortedNodes[x][0])
        G.remove_node(sortedNodes[x][0])
        #Gor[x].remove_node(sortedNodes[x][0])
        comp.append(x+1)
        cnt.append(len(list(nx.connected_components(G))[0]))
        print(x+1,len(list(nx.connected_components(G))[0]))
        #cnt.append(len(list(nx.connected_components(Gor[x]))))
        #print(x+1,len(list(nx.connected_components(Gor[x]))),len(list(nx.connected_components(Gor[x]))[0]))
        #print(Gor[x].order())

    fig, ax = plt.subplots()
    plt.bar(comp, cnt, width=0.50, color='r')

    plt.title("Vertex connectivity")
    plt.ylabel("Number of connected components")
    plt.xlabel("Rank in the list of highest degree nodes")
    ax.set_xticks([d  for d in comp])
    ax.set_xticklabels(comp)

    left, right = plt.xlim()  # return the current xlim
    plt.xlim(0, 30.35)
    print(left,right)

    plt.show()

#4,286,775
def edgeConnectivity(G):
    edgeList=list(G.edges())
    edgeData={}
    for x in edgeList:
        capacity = G.get_edge_data(x[0],x[1])
        edgeData[x] = int(capacity['capacity'])

    sortedEdgeList=reversed(sorted(edgeData.items(), key=lambda kv: kv[1])) #sort edges by capacity
    counter = [0]
    cnt = 0
    target = 1000
    comp=[2]
    for i in sortedEdgeList:
        G.remove_edge(i[0][0],i[0][1])
        if(comp[-1] < len(list(nx.connected_components(G)))):
            print(cnt+1)
        comp.append(len(list(nx.connected_components(G))))
        counter.append(cnt + 1)
        if cnt == target:
            break
        cnt+=1



    fig, ax = plt.subplots()
    plt.bar(counter, comp, width=1, color='r')

    plt.title("Edge connectivity")
    plt.ylabel("Number of connected components")
    plt.xlabel("Number of removed high-capacity edges")
    ax.set_xticks([d for d in counter])
    ax.set_xticklabels(counter)

    left, right = plt.xlim()  # return the current xlim
    plt.xlim(-10, 1010)

    n = 100
    for index, label in enumerate(ax.xaxis.get_ticklabels()):
        if index % n != 0 :
            label.set_visible(False)

    plt.show()

def weightsDistribution(G):
    global allWeights
    allWeights.sort(reverse=True)
    edgeCount = collections.Counter(allWeights)
    weight, cnt = zip(*edgeCount.items())

    fig, ax = plt.subplots()
    #plt.bar(weight, cnt, width=0.80, color='b')

    plt.hist(allWeights,100)
    plt.title("Weight Histogram")
    plt.ylabel("Count")
    plt.xlabel("Weights")
    ax.set_xticks([d  for d in weight])
    ax.set_xticklabels(weight)

    for index, label in enumerate(ax.xaxis.get_ticklabels()):
        if weight[index] not in [16777216,9000000,5000000,2000000,70000]:
            label.set_visible(False)

    plt.show()

def powerList(my_list, exponent):
    return [ x**exponent for x in my_list ]

def degreeDistribution(G):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    # print "Degree sequence", degree_sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    pk=[]
    for x in cnt:
        pk.append(x/G.order())

    plt.plot(np.array(deg), np.array(pk),  '.', color='red', label = 'Empyrical degree distribution')
    plt.plot(np.array(deg), np.array(powerList(deg,-2.1387317708757214)), '.', color='green', label='Scale-free approximation') #1.4425779698352683, -2.2495735

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
    plt.show()

def shortestPaths(G):
    shortPaths=nx.algorithms.shortest_paths.generic.shortest_path(G)
    print(shortPaths)

#draw graph in inset
def drawGraph(G):
     #plt.axes([0.4, 0.4, 0.5, 0.5])
     Gcc = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]
     pos = nx.spring_layout(G)
     #pos = nx.drawing.nx_pydot.graphviz_layout(G)
     fig, ax = plt.subplots()
     plt.axis('off')

     degrees = G.degree()
     nodes = G.nodes()
     n_color = np.asarray([degrees[n] for n in nodes])

     #nx.draw(G, nodelist=d.keys(), node_size=[v * 100 for v in d.values()])


     sc=nx.draw_networkx_nodes(G, pos, nodelist=nodes, ax=ax, node_size=n_color, node_color=n_color, cmap='viridis')
     nx.draw_networkx_edges(G, pos, alpha=0.2)

     sc.set_norm(mcolors.LogNorm())
     fig.colorbar(sc)

if __name__ == "__main__":
    main()

def defineGraph(data) -> object:
    global allMoneyStaked
    G = nx.Graph()
    for x in range(len(data)):
        G.add_edge(data[x]['node2_pub'], data[x]['node1_pub'], capacity=data[x]['capacity'])
        allWeights.append(int(data[x]['capacity']))
        allMoneyStaked+=int(data[x]['capacity'])
    return G


##https://graph.lndexplorer.com/api/graph
def readFile() -> object:
    with open('graph20190418.json') as f:
        data = json.load(f)
    return data['edges']

def smetric(G) -> object:
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    Gmax = li_smax_graph(list(degree_sequence))
    return smetricNonNormalized(G)/smetricNonNormalized(Gmax)

def smetricNonNormalized(G) -> object:
    return float(sum([G.degree(u) * G.degree(v) for (u, v) in G.edges()]))



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
                raise networkx.NetworkXError("THIS SHOULD NOT HAPPEN! - not graphable")
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
