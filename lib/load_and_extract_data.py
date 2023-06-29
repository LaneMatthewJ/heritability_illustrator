import pandas as pd
import networkx as nx
from sklearn.metrics import DistanceMetric
from scipy.spatial.distance import pdist, squareform
import numpy as np
import imgkit

from lib.data_loaders import read_hapmap
from lib.df_helpers import get_snps_within_frame, get_genotype_columns
from lib.heterozygote_to_homozygote import extract_heterozygous_columns, replace_heterozygotes_with_dom
from lib.haploid_squashing import squash_dip_to_hap
from lib.dominant_matching import create_dominant_match_df, calculate_matching_clusters
from lib.network_helpers import create_representative_set, extract_representatives, compute_serial_matrix
from lib.class_analysis import extract_names, convert_rep_class_map_to_df, calculate_pct_match_gene_group_to_classes


def get_non_dominant_strains(dataframe, max_by = 'count'):
	"""
	This grabs the non-dominant straits defined either by "max counts" which simply 
	plots the strains that are not in the largets grouping (almost always the maximum)
	but, the availability to plot the non-dominant strains by using the match df also
	exists by supplying `max_by` = `row`
	"""
	
	if max_by == 'count':
		counts = dataframe['counts'].apply(int)
		max_count = counts.max()
		values_below_max_count = list(dataframe[ counts < max_count ].index)
		return values_below_max_count

	if max_by == 'row': 
		snp_columns = list(filter(lambda x: x != 'counts', dataframe.columns))
		sums = dataframe[snp_columns].sum(axis=1)

		values_below_max_count = list(dataframe[sums < sums.max()].index)
		return values_below_max_count

def load_and_extract_data(filename, gene_name, start, stop, create_styled_table=True):
	df = read_hapmap(filename)
	
	framed_df = get_snps_within_frame(df, start, stop)
	genotype_columns = get_genotype_columns(df)

	homozygous_df = replace_heterozygotes_with_dom(framed_df, genotype_columns)
	haploid_df = squash_dip_to_hap(homozygous_df, genotype_columns)
	dominant_match = create_dominant_match_df(haploid_df, genotype_columns)
	matching_cluster_adj = calculate_matching_clusters(dominant_match)

	g = nx.from_pandas_adjacency(matching_cluster_adj)

	repset_list, element_map = create_representative_set(g)

	count_dict = extract_representatives(element_map, repset_list)

	dist = DistanceMetric.get_metric('hamming')
	distances = dist.pairwise(matching_cluster_adj[repset_list].T)

	ordered_dist_mat, res_order, res_linkage = compute_serial_matrix(distances,method='ward')

	ordered_repset = np.array(repset_list)[res_order]

	dominant_match.index = framed_df['pos']

	transposed = dominant_match[ordered_repset].T
	transposed['counts'] = [ f'{count_dict[x]}' for x in ordered_repset ] 

	styled_table = None
	outfile_name = None
	if create_styled_table: 
		outfile_name =  f'./styled_table_{gene_name}.png'
		styled_table = transposed.style.background_gradient(cmap='coolwarm')
		html = styled_table.render()
		imgkit.from_string(html, outfile_name)

	return styled_table, transposed, element_map, outfile_name

def analyze_output(filename, outfile_name, element_map): 
	class_encodings = pd.read_excel(filename) 
	class_encodings = extract_names(class_encodings)
	rep_df = convert_rep_class_map_to_df(element_map)
	df = calculate_pct_match_gene_group_to_classes(rep_df, class_encodings)

	df.to_csv(f'./{outfile_name}.tsv', sep='\t')

	return df, class_encodings