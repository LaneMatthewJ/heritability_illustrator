import pandas as pd
import numpy as np


def convert_vcf_to_hapmaplike(input_file, output_filename): 
    """
    Reads in the input file and converts it to a hapmap-esque 
    file

    
    """
    df_orig = pd.read_csv(input_file, sep='\t', skiprows=42)
    df = df_orig.copy()
    
    data_cols = df.columns[9:]

    # VCF files contain combinations of 0, 1, 2, 3, etc. 
    # 0 refers to the reference genome. 
    # The slash assumes diploid, but as these are haploid creatures, 
    # we take only the first value. 
    
    for column in data_cols: 
        df[column] = df_orig[column].apply(lambda x: 'AA' if x[0] == '0' else 'CC')

    df['rs#'] = df[["#CHROM", "POS"]].apply(lambda x: f"{x['#CHROM']}:{x['POS']}", axis=1)
    df['alleles'] = "A/C" # a=0, c=1
    df['chrom'] = 0
    df['pos'] = df['POS']
    df['strand'] = 0
    df['assembly#'] = np.NaN
    df['center'] = np.NaN
    df['protLSID'] = np.NaN
    df['assayLSID'] = np.NaN
    df['panelLSID'] = np.NaN
    df['QCcode'] = np.NaN
    
    desired_cols = ['rs#', 'alleles', 'chrom', 'pos', 'strand', 'assembly#', 'center',
       'protLSID', 'assayLSID', 'panelLSID', 'QCcode', *data_cols]

    
    df[desired_cols].to_csv(output_filename, sep='\t', index=None)