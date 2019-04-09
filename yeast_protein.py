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




Propagation_matrix = []
colored_nodes = []

def show_graph_with_labels(adjacency_matrix,size = 1):
	length = len(dataset[(list(dataset)[0])].tolist())
	edges =[]
	print dataset
	G = nx.Graph()
	for x in xrange(0,length):
		source = (dataset[(list(dataset)[0])].tolist())[x]
		target = (dataset[(list(dataset)[1])].tolist())[x]
		edges.append((source,target))

	G.add_edges_from(edges)
	degrees = G.degree() #Dict with Node ID, Degree
	nodes = G.nodes()
	n_color = np.asarray([degrees[n] for n in nodes])
	
	print "\n\nmin degree ="+ str(min(n_color))
	print "max degree ="+ str(max(n_color))
	print str([x == min(n_color)  for x in n_color].count(1))+ " nodes have minimum number of degrees"
	print str([x == max(n_color)  for x in n_color].count(1))+ " nodes nodes have maximum number of degrees"

	plt.figure("Protein interactions in Yeast")
	plt.subplot(211)
	nx.draw(G,node_color=n_color, nodelist = G.nodes(),cmap='coolwarm',with_labels=False,node_size= n_color,width = 0.1)

	maxitem = 0
	for x in n_color:
		if x == max(n_color):
			maxitem = x
			break

	print maxitem
	start = maxitem


	Affected = []

	for i in xrange(0,max(nodes)+1):
		Propagation_matrix.append(0)
		Affected.append(0)

	propagate_color(start,edges)
	print "Propagation from point " + str(start) + "  affected "+str(Propagation_matrix.count(1))
	Affected[maxitem] = Propagation_matrix.count(1)
	prop_array = np.asarray([p for p in Propagation_matrix])
	n_color =  np.asarray([degrees[n] for n in colored_nodes])
	plt.subplot(212)
	for e in edges:
		if (e[0] not in colored_nodes ) and (e[1] not in colored_nodes):
			G.remove_edge(e[0],e[1])
	nx.draw(G,node_color=n_color, nodelist = colored_nodes,cmap='coolwarm',with_labels=False,node_size= n_color,width = 0.1)
	plt.show()





	for start in nodes:
		for i in xrange(0,max(nodes)+1):
			Propagation_matrix[i] = 0
		propagate_color(start,edges)
		print "Propagation from point " + str(start) + "  affected "+str(Propagation_matrix.count(1)) + " nodes"
		Affected[start] = Propagation_matrix.count(1)

	plt.xlabel("Starting Node")
	plt.ylabel("Nodes affected")
	plt.xlim(-1,2000)
	plt.ylim(-1,2000)
	plt.plot(Affected,'ro')	
	plt.title("Number of nodes that can be affected from different starting points ")
	plt.show()

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
	print(list(dataset))
	if len(list(dataset)) > 2 :
		min_index = 2
		if list(dataset)[0]== '%':
			print("drop %")
			dataset.drop(['%'],axis=1,inplace = True)
			min_index = 2
		trange = len(list(dataset))
		#print list(dataset)
		for x in xrange(min_index,trange):
			col = list(dataset)[min_index] 
			print("drop "+col)
			dataset.drop(col,axis=1, inplace = True)
		#print(list(dataset))
	
	#print dataset

	show_graph_with_labels(dataset,5)


	
	
	#print(length)
	
#	CSV = []

#	for x in xrange(0,length):
#		source = (dataset[(list(dataset)[0])].tolist())[x]
#		target = (dataset[(list(dataset)[1])].tolist())[x]
#		print source,target
#		CSV[source] = target
	
#	print CSV
	
	