from .load_and_extract_data import load_and_extract_data, analyze_output, get_non_dominant_strains
from .plot_images import plot_representative_set
from .create_gene_region_profile import generate_profile
import time
import os
import shutil
import pandas as pd

def replace_0s_and_1s(row, desired_columns): 
    ref = row['REF']
    alt = row['ALT']

    for col in desired_columns: 
        row[col] = row[col].replace('0', ref)
        row[col] = row[col].replace('1', alt)
    return row

def calculate_start_stops(x, gene_start, gene_stop):
    """
    This function is used to calculate the "starts" and "ends" for circos plots.
    Instead of showing the reference one color and the alternate another, the 
    reference color is left empty to not crowd the space.
    """
    
    start_stop_dict = {
        'indv': [],
        'start': [],
        'end': [], 
        'insertion_size': []
    }
    
    start = -1 
    stop = -1 
    previous_matching = -1
    previous_index = -1
    
    for i in range(len(x)):
        index = x.index[i]
        matching = x[index]
        if i == 0 and matching == 0:
            start = index# gene_start
            previous_index = start
            previous_matching = matching
            continue
        
        if matching == 0 and start == -1: 
            start = index
            previous_matching = matching
            continue
            
        if matching == 1 and previous_matching == 0 and stop == -1: 
            stop = index
            start_stop_dict['start'].append(start)
            start_stop_dict['end'].append(stop)
            start_stop_dict['insertion_size'].append(int(stop)-int(start))
            
            start = -1
            stop = -1
            previous_matching = 1
            continue
        
    if start != -1 and previous_matching == 0:  
        stop = index            
        start_stop_dict['start'].append(start)
        start_stop_dict['end'].append(stop)
        
        
        start_stop_dict['insertion_size'].append(int(stop)-int(start))

    if len(start_stop_dict['start']) > 0: 
        start_stop_dict['indv'] = x.name
        
    return pd.DataFrame(start_stop_dict)
# data.apply(lambda x: calculate_start_stops(x), axis=1)

def add_circular_insertion_to_create_start_stop(df, gene_start, gene_stop, percentage):
    """
    If there is an insertion near the end of the gene region and
    one near the beginning (within the specified percentage of the total region),
    this function adds a new "bridging" insertion that circulates around zero.
    
    Parameters:
        df (pd.DataFrame): DataFrame with columns ['indv', 'start', 'end', 'insertion_size'].
        gene_start (int): The start coordinate of the gene region.
        gene_stop (int): The stop coordinate of the gene region.
        percentage (float): A value between 0 and 1 that defines the threshold 
                            distance (as a fraction of the gene length) to consider
                            an insertion as near the boundary.
                            
    Returns:
        pd.DataFrame: Updated DataFrame with the bridging insertion added if conditions are met.
    """
    # Calculate the total gene region length and threshold distance.
    total_length = gene_stop - gene_start
    threshold = percentage * total_length

    # If the DataFrame is empty, return it unchanged.
    if df.empty:
        return df

    # Assume the DataFrame is already sorted by the start coordinate.
    first_start = df['start'].iloc[0]
    last_end = df['end'].iloc[-1]

    # Check if the last insertion is close to the gene_stop and the first is close to the gene_start.
    if gene_start == 0 and ((gene_stop - last_end) <= threshold) and ((first_start - gene_start) <= threshold):
        # For a circular region, the bridging insertion starts at the insertion from the end
        # and ends at the beginning insertion's end adjusted by the total gene length.
        bridging_start = df['end'].iloc[-1]
        bridging_end = df['start'].iloc[0]

        # Create a new row for the bridging insertion.
        new_rows = pd.DataFrame([{
            'indv': df['indv'].iloc[0],  # Assuming all rows refer to the same individual.
            'start': bridging_start,
            'end': gene_stop,
            'insertion_size': gene_stop - bridging_start
        }, 
        {
            'indv': df['indv'].iloc[0],  # Assuming all rows refer to the same individual.
            'start': 0,
            'end': bridging_end,
            'insertion_size': bridging_end
        }])

        # Append the new row to the DataFrame.
        df = pd.concat([df, new_rows])

        # Optionally, sort the DataFrame if needed.
        df.sort_values(by='start', inplace=True)

    return df


def calculate_insertions_with_overlap(x,  gene_start, gene_stop, pct): 
    out = calculate_start_stops(x,  gene_start, gene_stop)
    out['start'] = out['start'].apply(int)
    out['end'] = out['end'].apply(int)
    out = add_circular_insertion_to_create_start_stop(out, gene_start, gene_stop, pct)
    return out 
    
def create_start_stop_df(matchings_by_group_file_path, gene_start, gene_stop, output_file, pct):
    start_stop_data = pd.read_csv(
        matchings_by_group_file_path, 
        sep='\t', 
        index_col=0)
    
    no_counts = start_stop_data[start_stop_data.columns[:-1]]
    out = no_counts.apply(
        lambda x: calculate_insertions_with_overlap(x,  gene_start, gene_stop, pct), 
        axis=1
    )
    output_df = pd.concat(out.values)

    
    output_df['start'] = output_df['start'].apply(int)
    output_df['end'] = output_df['end'].apply(int)
    sorted_start_stop = output_df.sort_values(by='start')
    
    sorted_start_stop.to_csv(f"{output_file}_total_start_stop_df.tsv", index=False)

    grouped = output_df[['indv','insertion_size']].groupby('indv').sum()
    mean_size = grouped['insertion_size'].mean()
    left_group = grouped[grouped['insertion_size'] < mean_size] 
    right_group = grouped[grouped['insertion_size'] > mean_size]

    left_df = output_df[output_df['indv'].isin(left_group.index)]
    left_df.to_csv(f'{output_file}_total_start_stop_df_MINOR.tsv', index=False)
    right_df = output_df[output_df['indv'].isin(right_group.index)]
    right_df.to_csv(f'{output_file}_total_start_stop_df_MAJOR.tsv', index=False)


def read_vcf_and_create_profile(file_path, output_name): 
    df = pd.read_csv(file_path, skiprows=42, sep='\t') 
    # 1. extract the columns 
    desired_columns = df.columns[9:]
    for col in desired_columns: 
        df[col] = df[col].apply(lambda x: x.split(":")[0])
        df[col] = df[col].apply(lambda x: x.replace('/',''))
    
    df2 = df.apply(lambda x: replace_0s_and_1s(x, desired_columns), axis=1)
    
    df2.tail() # 51997	
    df2['alleles'] = df['REF'] + '/' + df['ALT']
    df2.columns = list(map(lambda x: x.lower(), df2.columns) )
    df2.to_csv(output_name, sep='\t', index=False) 
    return df2


def generate_profile_from_vcf(vcf_filepath, start_pos, stop_pos, class_encoding_filename, input_image_directory):
    base_name = os.path.basename(vcf_filepath)
    no_ext = base_name.split(".vcf")[0]
    output_name = f"{no_ext}_profile.tsv"
    file_name = 'v'
    read_vcf_and_create_profile( vcf_filepath, output_name )
    
    gene_name = os.path.basename(no_ext)
    
    output = generate_profile(
        output_name,
        gene_name,
        start_pos,
        stop_pos, 
        create_styled_table=False,
        class_encodings_filename=class_encoding_filename, 
        input_image_directory = input_image_directory, 
        create_image_plots=False
    )
    return output
