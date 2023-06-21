import pandas as pd



def get_snps_within_frame(df: pd.DataFrame, start: int, stop: int):
    """
    assumes df has a position column
    """
    ge_mask = df['pos'] > start
    lt_mask = df['pos'] < stop

    sliced = df[ ge_mask & lt_mask ].copy()
    
    return sliced


def get_genotype_columns(df: pd.DataFrame, colstart: int = 11, colstop:int = -2):
	"""
	This function asumes genotypes start at the 11th column. 
	This function additionally assumes that the dominant and non dominant 
	snps have been added as the final 2 columns

	Otherwise, the columns will need to be supplied.
	"""

	genotype_columns = df.columns[colstart:colstop]


	return genotype_columns