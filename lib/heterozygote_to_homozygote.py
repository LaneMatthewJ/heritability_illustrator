import pandas as pd

"""
## Heterozygotes
Only a few heterozygotes: using dominant nucleotide 
(i.e. the first nucleotide form the alleles column 
"M/N" - "M" is the dominant)
"""


def extract_heterozygous_columns(df: pd.DataFrame, genotype_columns: list): 
	"""
	This function finds the columns that are heterozygous
	"""

	heterozygous_df = df[genotype_columns].copy()
	for genotype in genotype_columns:
		heterozygous_df[genotype] = df[genotype].apply( lambda x: x[0] != x[1] )     
		
	columns_of_heterozygosity = heterozygous_df.columns[ heterozygous_df.sum() == 1 ] 
	
	return columns_of_heterozygosity


def replace_heterozygotes_with_dom(df: pd.DataFrame, genotype_columns: list): 
	"""
	Replace heterozygotes w/ dominant references. 
	"""
	columns_of_heterozygosity = extract_heterozygous_columns(df, genotype_columns)

	for genotype in columns_of_heterozygosity:
		df[genotype] = df[genotype].apply( lambda x: x[1]+x[1] if x[0] != x[1] else x )

	return df


