import pandas as pd

def extract_names(
	class_encodings, 
	column_name = 'Unnamed: 0', 
	name_split_fn = lambda x: x.split('_')[3].split('.')[0].upper()): 
	
	names_series  = class_encodings[column_name]
	updated_names = names_series.apply(name_split_fn)
	class_encodings['name'] = updated_names

	return class_encodings

def convert_rep_class_map_to_df(rep_dict): 
	return pd.DataFrame({
		'class': rep_dict.values(),
		'elements': rep_dict.keys()
	})


def calculate_pct_match_gene_group_to_classes(rep_df, class_encodings, non_sample_cols= ['Unnamed: 0', 'name']): 
	"""
	This function takes in a repset of genes (i.e. which genotypes align on this particular frame)
	and a user defined class_encoding list. 

	Each user defined classes: 
		Each of the major groups of matching snp profiles is tested to see how many within 
		that snp profile group matches to the class. 
	"""

	class_columns = list(
		filter( 
			lambda x: x not in non_sample_cols, 
			class_encodings.columns)
		)

	df_list = {
		'class_name': [], 
		'group_key':[], 
		'n_matched':[], 
		'n_total':[], 
		'matched_pct':[],
		'matched': []
	}

	for class_name in class_columns: 
		matched = class_encodings[class_name]
		samples_in_class = class_encodings['name'].loc[matched.dropna().index]


		for unique_val in rep_df['class'].unique(): 
			samples_in_class_matching_unique_group = rep_df[ rep_df['class'] == unique_val ]['elements']
			samples_in_unique_repset = samples_in_class_matching_unique_group.isin(samples_in_class)
			n_samples_in_unique_repset = samples_in_unique_repset.sum()
			total_samples_in_whole_set =  rep_df[ rep_df['class'] == unique_val ]['elements'].isin(class_encodings['name']).sum()

			matched = '' if n_samples_in_unique_repset == 0 else '|'.join(list(samples_in_class_matching_unique_group.values))

			df_list['class_name'].append(class_name)
			df_list['group_key'].append(unique_val)
			df_list['n_matched'].append(n_samples_in_unique_repset)
			df_list['n_total'].append(total_samples_in_whole_set)
			df_list['matched_pct'].append(n_samples_in_unique_repset / total_samples_in_whole_set)
			df_list['matched'].append(matched)

	total_vals_df =  pd.DataFrame(df_list)
	return total_vals_df[ total_vals_df['n_matched'] > 0 ]
