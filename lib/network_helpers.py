import random
import numpy as np
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from fastcluster import linkage

def create_representative_set(network: nx.Graph) -> (list, dict): 
    random.seed(42)
    """
    This function takes a NX network object, takes a list of all connected components, 
    and then chooses the first element as a representative for that connected component.
    """
    S = list(nx.connected_components(network))
    representatives = []
    element_map = {}
    for subgraph in S:
        representative = list(subgraph)[0]
        representatives.append(representative)

        # Get the remainder of elements within S1: 
        for element in subgraph:
            # if element == representative: 
            #     continue

            element_map[element] = representative

    return representatives, element_map


def extract_representatives(element_map, representative_list): 
	unique_map_values, unique_map_counts = np.unique([ *element_map.values()], return_counts=True)
	count_dict = {}
	for i in range(len(unique_map_values)): 
		count_dict[unique_map_values[i]] = unique_map_counts[i]
		
	for rep in representative_list: 
		if rep not in count_dict.keys():
			count_dict[rep] = 1
	
	return count_dict

def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))
    
def compute_serial_matrix(dist_mat,method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = np.zeros((N,N))
    a,b = np.triu_indices(N,k=1)
    seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b,a] = seriated_dist[a,b]
    
    return seriated_dist, res_order, res_linkage
