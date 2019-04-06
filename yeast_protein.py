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
	plt.title("Protein interactions in yeast")
	plt.xlabel("protein1")
	plt.ylabel("protein2")
	#nx.draw(G,color_map = mapper_list, node_size=[v*v*size +1 for v in M],with_labels=False, width = 0.1 ,arrowsize=0.1)
	degrees = G.degree() #Dict with Node ID, Degree
	nodes = G.nodes()
	#print nodes 
	#exit(0)
	n_color = np.asarray([degrees[n] for n in nodes])
	
	print "\n\nmin degree ="+ str(min(n_color))
	print "max degree ="+ str(max(n_color))
	print str([x == min(n_color)  for x in n_color].count(1))+ " nodes have minimum number of degrees"
	print str([x == max(n_color)  for x in n_color].count(1))+ " nodes nodes have maximum number of degrees"


	maxitem = 0
	for x in n_color:
		if x == max(n_color):
			maxitem = x
			break

	print maxitem


	#for i in xrange(0,len(nodes)):
	#	Propagation_matrix.append(0)
	#propagate_color(maxitem,edges)
	#print "Propagation from max affected "+str(Propagation_matrix.count(1))


	start = maxitem
	for i in xrange(0,len(nodes)):
		Propagation_matrix.append(0)
	propagate_color(start,edges)
	print "Propagation from point " + str(start) + "  affected "+str(Propagation_matrix.count(1))

	prop_array = np.asarray([p for p in Propagation_matrix])
	
	#nx.draw(G,node_color=prop_array, cmap='coolwarm',with_labels=False,node_size= prop_array,width = 0.1)
	nx.draw(G,node_color=n_color, cmap='coolwarm',with_labels=False,node_size= n_color,width = 0.1)
	plt.show()


def propagate_color(node,edges):
	global Propagation_matrix
	print "propagate to" +str(node)
	Propagation_matrix[node] = 1
	neighbors = []
	for e in edges:
		if e[0] ==node:
			neighbors.append(e[1])
		elif e[1] == node:
			neighbors.append(e[0])

	for e in neighbors:
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
	
	