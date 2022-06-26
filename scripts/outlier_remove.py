#! /usr/bin/env python
import kmeans1d
import sys

# tab seperated no header
# first columns is sequence naem
# second is cluster id
# third is marker occupancy

clusters = dict()
with open(sys.argv[1]) as f:
	for line in f.readlines():
		l_elem = line.strip().split("\t")
		name = l_elem[0]
		clust = int(l_elem[1])
		numcopy = int(l_elem[2])
		if clust == -1:
			continue
		if clust in clusters:
			clusters[clust]+= [(numcopy,name)]
		else:
			clusters[clust]= [(numcopy,name)]

clusters_sorted = {k:sorted(v) for k,v in clusters.items()}

k=2
for key, v in clusters_sorted.items():	
	while True:
		clusters, centroids = kmeans1d.cluster([i[0] for i in v], k)
		if clusters[0] > 0:
			break
		rat = (centroids[1]-centroids[0])/centroids[1]
		if rat > 0.5:
			print(str(key) + "\t" + v[0][1])
			v = v[1:]
		else:
			break


