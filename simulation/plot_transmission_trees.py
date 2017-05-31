#!/usr/bin/env python
import networkx as nx
import pygraphviz as pgv
from optparse import OptionParser

if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-f', '--filepath', action = 'store', type='string', dest = 'pathDir', default = '.' )
	
	opts, args = parser.parse_args()
	
	#plot true tree
	file = opts.pathDir+'/input/siLog.txt'
	f = open(file, 'r')
	header = f.next() #patientid	t_inf	source	source_type	t_sample
	log = [l.strip().split() for l in f]
	cases = [int(l[0]) for l in log]
	
	G = pgv.AGraph(strict=False,directed=True)
	G.add_nodes_from(cases)
	
	for l in log:
		if l[3]!='0':
			if l[3]=='1':
				#ward link
				G.add_edge(int(l[2]), int(l[0]), color='red', label=l[1])
			elif l[3]=='2':
				 #i.e. hospital link
				G.add_edge(int(l[2]), int(l[0]), color='blue', label=l[1])
			else:
				#spore link
				G.add_edge(int(l[2]), int(l[0]), color='green', label=l[1])

	G.edge_attr['len']=1.3
	G.layout()
	G.draw(opts.pathDir+'/input/true_transmisson_tree.pdf')