import networkx as nx
import numpy as np
import pandas as pd
import random

def create_dominant_match_df(df, genotype_columns): 
	matching_df = df[genotype_columns].copy()
	for column in genotype_columns: 
		matching_df[column] = df['dominant'] == df[column]
		matching_df[column] = matching_df[column].apply(lambda x: 1 if x else 0) 

	return matching_df


def calculate_matching_clusters(matching_df):
	"""
	This creates a matrix in which we see which sequences match which other sequences: 
	e.g. 
			  a b c d			  a b c d
			  1 1 1 1			a 0 0 1 0
			  0 1 0 1			b 0 0 0 1
			  0 1 0 1			c 1 0 0 0
			  0 0 0 0			d 0 1 0 0	
	"""
	matching_distance_mat = np.zeros((len(matching_df.columns), len(matching_df.columns)))

	for column_idx_i in range(len(matching_df.columns)): 
		col_i = matching_df.columns[column_idx_i]
		for column_idx_j in range(column_idx_i+1, len(matching_df.columns)): 
			col_j = matching_df.columns[column_idx_j]
			match = all(matching_df[col_i] == matching_df[col_j])
			matching_distance_mat[column_idx_i,column_idx_j] = match
			
	adj = pd.DataFrame(matching_distance_mat, columns=matching_df.columns)

	adj.index = adj.columns

	return adj
