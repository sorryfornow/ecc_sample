# ecc_sample
Eccentricity calculate for temporal graph.

Eccentricity definition: the ecc(v) = max(dist(u,v)) for all vertices u in graph G.

For a temporal graph G with timestamp in [1,max], this code calculate the query ecc(v)_[i,j] for the subgraph of any given interval [i,j] with i,j in range[1,max], where v is a arbitrary vertex in G.

The code IFECC is replicated the study:
On Scalable Computation of Graph Eccentricities SIGMOD'22

The .PY code in tempG_generation is writen for pre-processing
