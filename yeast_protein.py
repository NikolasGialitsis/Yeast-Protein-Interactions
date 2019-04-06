#!/usr/bin/env python
import pandas as pd
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import networkx as nx
import sys

def show_graph_with_labels(adjacency_matrix,size = 1):
	length = len(dataset[(list(dataset)[0])].tolist())
	edges =[]
	print dataset
	G = nx.DiGraph()
	for x in xrange(0,length):
		source = (dataset[(list(dataset)[0])].tolist())[x]
		print x
		target = (dataset[(list(dataset)[1])].tolist())[x]
		#print source,target
		edges.append((source,target))

	color_map = ['orange']
	G.add_edges_from(edges)
	d = G.in_degree()

	M = []
	for item in d:
		M.append(item)
	#plt.figure(figsize=(50,50))
	plt.title("Protein interactions in yeast")
	plt.xlabel("protein1")
	plt.ylabel("protein2")
	nx.draw_circular(G,node_color = color_map,node_size=[v[1]*500*size for v in M],with_labels=True, width = 1,arrowsize=1)
	plt.show()


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
		print list(dataset)
		for x in xrange(min_index,trange):
			col = list(dataset)[min_index] 
			print("drop "+col)
			dataset.drop(col,axis=1, inplace = True)
		#print(list(dataset))
	
	#print dataset

	show_graph_with_labels(dataset,1)


	
	
	#print(length)
	
#	CSV = []

#	for x in xrange(0,length):
#		source = (dataset[(list(dataset)[0])].tolist())[x]
#		target = (dataset[(list(dataset)[1])].tolist())[x]
#		print source,target
#		CSV[source] = target
	
#	print CSV
	
	