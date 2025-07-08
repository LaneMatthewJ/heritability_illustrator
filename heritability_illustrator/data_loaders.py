import pandas as pd


def add_dominant_nondominant_to_df(df: pd.DataFrame) -> pd.DataFrame: 
	""" 
	This function extracts dominant and nondominant snps and stores them
	in their own respective columns
	"""
	df['dominant'] = df['alleles'].apply(lambda x: x.split('/')[0])
	df['recessive'] = df['alleles'].apply(lambda x: x.split('/')[1])

	return df 

def read_hapmap(filepath: str) -> pd.DataFrame: 
	"""
	This function reads in a hapmap file and applies dominant and extracts
	dominant and non-domant snps and pasts them into the file. 
	"""

	df = pd.read_csv(filepath, sep='\t') 

	df = add_dominant_nondominant_to_df(df)
	return df
