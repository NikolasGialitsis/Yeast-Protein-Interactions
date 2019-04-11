#!/usr/bin/env python
import pandas as pd
import copy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import networkx as nx
import sys
from scipy.optimize import curve_fit

def powerlaw(x, alpha, beta, x0):
	if beta < 0 :
		beta = -beta
	return (x + x0)**alpha * np.exp(-beta *x)



Propagation_matrix = []
colored_nodes = []


G = nx.Graph()


def show_degree_centrality(G):


	degrees = G.degree() #Dict with Node ID, Degree
	nodes = G.nodes()
	edges = G.edges()
	n_color = np.asarray([degrees[n] for n in nodes])
	
	print "\n\nmin degree ="+ str(min(n_color))
	print "max degree ="+ str(max(n_color))
	print str([x == min(n_color)  for x in n_color].count(1))+ " nodes have minimum degree centralities"
	print str([x == max(n_color)  for x in n_color].count(1))+ " nodes nodes have maximum degree centralities"	
	maxitem = 0
	for x in n_color:
		if x == max(n_color):
			break
		maxitem+=1
	start = maxitem
	print 'node with maximum degrees = ' + str(maxitem) + " with num of degrees = "+ str(max(n_color))
	
	fixed_positions = {1: (maxitem,max(n_color))}#dict with two of the positions set
	fixed_nodes = fixed_positions.keys()
	pos = nx.spring_layout(G,pos=fixed_positions, fixed = fixed_nodes)

	plt.figure("Protein interactions in Yeast")
	plt.subplot(211)
	nx.draw(G,pos = pos ,node_color=n_color, nodelist = G.nodes(),cmap='coolwarm',with_labels=False,node_size= n_color,width = 0.1)
	
	Affected = []
	
	for i in xrange(0,max(nodes)+1):
		Propagation_matrix.append(0)
		Affected.append(0)

	propagate_color(start,edges)
	print "Propagation from point " + str(start) + "  affected "+str(Propagation_matrix.count(1)) + " nodes"
	Affected[maxitem] = Propagation_matrix.count(1)
	prop_array = np.asarray([p for p in Propagation_matrix])
	n_color_2 =  np.asarray([degrees[n] for n in colored_nodes])
	plt.subplot(212)
	G2 = G.copy()
	for e in edges:
		if (e[0] not in colored_nodes ) and (e[1] not in colored_nodes):
			G2.remove_edge(e[0],e[1])
	nx.draw(G2,pos = pos,node_color=n_color_2, nodelist = colored_nodes,cmap='coolwarm',with_labels=False,node_size= n_color_2,width = 0.1)
	plt.show(block=True)




def show_color_propagation(G):
	global Propagation_matrix

	Affected = []
	nodes = G.nodes()
	edges = G.edges()
	if Propagation_matrix == []:
		for i in xrange(0,max(nodes)+1):
			Propagation_matrix.append(0)
			Affected.append(0)

	for start in G.nodes():
		for i in xrange(0,max(nodes)+1):
			Propagation_matrix[i] = 0
		propagate_color(start,edges)
		print "Propagation from point " + str(start) + "  affected "+str(Propagation_matrix.count(1)) + " nodes"
		Affected[start] = Propagation_matrix.count(1)

	plt.xlabel("Starting Node")
	plt.ylabel("Nodes affected")
	plt.xlim(-1,2000)
	plt.ylim(-1,2000)
	plt.title("Number of nodes that can be affected from different starting points ")
	plt.plot(Affected,'ro')	


def show_degree_distribution(G):
	nodes = G.nodes()
	degrees = G.degree()
	n_color = np.asarray([degrees[n] for n in nodes])
	plt.figure("Power Law - degrees frequency")
	plt.xlabel("Nodes")
	plt.ylabel("Degrees")
	x ,y  = np.unique(n_color, return_counts=True)
	plt.scatter(x,y,label = "degree count")
	parameters, fit  = curve_fit(powerlaw, x, y,maxfev = 3000)
	plt.plot(x, powerlaw(x, *parameters), label = "power_law",color = "r")
	plt.legend()
	plt.show()


def show_graph_info(G):
	nodes = G.nodes()
	edges = G.edges()
	#print "Degree centrality = "+str(nx.algorithms.degree_centrality(G))
	print "Graph density = " + str(nx.density(G))
	print "Number of edges : " + str(len(edges))
	print "Number of nodes: " + str(len(nodes))

def show_closeness_centrality(G):
	nodes = G.nodes()	
	closeness_centralities = []
	degrees = G.degree()
	#print nodes
	for n in nodes:
		closeness_centralities.append(nx.closeness_centrality(G, u=n))

	plt.figure("closeness centralities")
	plt.xlabel("nodes")
	plt.ylabel("degree centrality")

	print str([x == min(closeness_centralities)  for x in closeness_centralities].count(1))+ " nodes have minimum closeness centralities"
	print str([x == max(closeness_centralities)  for x in closeness_centralities].count(1))+ " nodes have maximum closeness centralities"

	max_centrality_node = 0
	for x in closeness_centralities:
		if x == max(closeness_centralities):
			break
		max_centrality_node+=1

	print "max centrality node = "+str(max_centrality_node)+ " with value : "+str(max(closeness_centralities)) + " and degree = " \
		+str(degrees[max_centrality_node])
	plt.scatter(max_centrality_node,max(closeness_centralities),color = "red",label = "best node",s = 20)
	plt.scatter(nodes,closeness_centralities,color = "blue" , label = "centralities",s = 2)
	plt.legend()
	plt.show(block=True)

	if len(closeness_centralities) != len(nodes):
		print 'error size'
		exit(-1)
	G3 = G.copy()
	edges = G.edges()
	for e in edges:
		if closeness_centralities[e[0]-1] == 0.0:
			G3.remove_edge(e[0],e[1])
	sizes = [x * 50 for x in closeness_centralities]
	plt.figure("Closeness centralities Graph")
	nx.draw(G3,node_color= closeness_centralities ,nodelist = nodes,cmap='coolwarm',with_labels=False, \
		node_size= sizes, width = 0.1)
	plt.show(block=True)



def show_betweeness_centrality(G):
	nodes = G.nodes()
	degrees = G.degree()
	temp = nx.betweenness_centrality(G)
	inbetween_centralities = []
	for x in temp:
		inbetween_centralities.append(temp.get(x))

	plt.figure("inbetween centralities")
	plt.xlabel("nodes")
	plt.ylabel("degree centrality")

	print str([x == min(inbetween_centralities)  for x in inbetween_centralities].count(1))+ " nodes have minimum closeness centralities"
	print str([x == max(inbetween_centralities)  for x in inbetween_centralities].count(1))+ " nodes have maximum closeness centralities"

	max_centrality_node = 0
	for x in inbetween_centralities:
		if x == max(inbetween_centralities):
			break
		max_centrality_node+=1

	print "max centrality node = "+str(max_centrality_node)+ " with value : "+str(max(inbetween_centralities)) + " and degree = " \
		+str(degrees[max_centrality_node])
	plt.scatter(max_centrality_node,max(inbetween_centralities),color = "red",label = "best node",s = 20)
	plt.scatter(nodes,inbetween_centralities,color = "blue" , label = "centralities",s = 2)
	plt.legend()
	plt.show(block=True)

	if len(inbetween_centralities) != len(nodes):
		print 'error size'
		exit(-1)

	G4 = G.copy()
	edges = G.edges()
	for e in edges:
		if inbetween_centralities[e[0]-1] == 0.0:
			G4.remove_edge(e[0],e[1])

	sizes = [x * 100 for x in inbetween_centralities]
	plt.figure("inbetweeness centralities Graph")
	nx.draw(G4,node_color= inbetween_centralities ,nodelist = nodes,cmap='coolwarm',with_labels=False, \
		node_size= sizes, width = 0.1)
	plt.show(block=True)



def show_pagerank(G):

	nodes = G.nodes()
	ranks = nx.pagerank(G)
	#print ranks
	plt.figure("PageRank")

	x = [r for r in ranks ]
	y = [ranks[r] for r in ranks ]
	print "maximum page rank  = " + str(max(y))

	best_node = 0
	for k in x :
		if ranks[k] == max(y):
			print "best page rank node =  " + str(k)
			best_node = k
			break
	
	plt.scatter(x,y,color = "blue", s = 1)
	plt.scatter(best_node,max(y),color = "red" , s = 5)
	plt.show()
	

	plt.figure("PageRank")
	yarray = np.asarray([k*5000 for k in y])
	xarray = np.asarray([k for k in x])

	nx.draw(G,node_color= yarray ,cmap='coolwarm',with_labels=False,node_size= yarray, width = 0.1)
	plt.show()


def initialize_network(G,dataset):
	
	length = len(dataset[(list(dataset)[0])].tolist())
	edges =[]
	for x in xrange(0,length):
		source = (dataset[(list(dataset)[0])].tolist())[x]
		target = (dataset[(list(dataset)[1])].tolist())[x]
		edges.append((source,target))

	G.add_edges_from(edges)



def propagate_color(node,edges):
	global Propagation_matrix
	global colored_nodes
	#print "propagate to  " +str(node)
	Propagation_matrix[node] = 1
	colored_nodes.append(node)
	neighbors = []
	for e in edges: 
		if e[0] ==node:
			neighbors.append(e[1])
		elif e[1] == node:
			neighbors.append(e[0])

	for e in neighbors:
		#print "neighbor = "+ str(e)
		if Propagation_matrix[e] == 1 :
			continue
		else:
 			propagate_color(e,edges)


 		


if __name__ == '__main__':  

	csv_file = ""
	if len(sys.argv) == 1 or (len(sys.argv) == 2 and str(sys.argv[0]) == "python"):
		csv_file = "./out.protein_protein"
	else:
		if str(sys.argv[0]) == "python":
			csv_file = str(sys.argv[2])
		else:
			csv_file = str(sys.argv[1])
	dataset=pd.read_csv(csv_file,delimiter= " ",header= 1)
	if len(list(dataset)) > 2 :
		min_index = 2
		if list(dataset)[0]== '%':
			print("drop %")
			dataset.drop(['%'],axis=1,inplace = True)
			min_index = 2
		trange = len(list(dataset))
		for x in xrange(min_index,trange):
			col = list(dataset)[min_index] 
			print("drop "+col)
			dataset.drop(col,axis=1, inplace = True)

	initialize_network(G,dataset)
	'''
	show_degree_distribution(G)
	show_graph_info(G)
	show_pagerank(G)

	show_betweeness_centrality(G)
	show_closeness_centrality(G)
	show_color_propagation(G)
	'''
	show_degree_centrality(G)
	