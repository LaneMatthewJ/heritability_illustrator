import pandas as pd

def squash_dip_to_hap(df: pd.DataFrame, genotype_columns: list): 
	"""
	Changes diploid data to haploid data
	"""
	for column in genotype_columns: 
		df[column] = df[column].apply( lambda x: x[1])

	return df



