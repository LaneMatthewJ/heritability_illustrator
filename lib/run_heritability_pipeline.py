from lib.load_and_extract_data import load_and_extract_data, analyze_output, get_non_dominant_strains
from lib.plot_images import plot_representative_set
from lib.create_gene_region_profile import generate_profile
import time
import os
import shutil
import pandas as pd




def run_pipeline(): 



file_name = '../ctherm/test_ctherm_DSM1313_fromvcf.tsv'
gene_name = 'ctherm_DSM1313_test_fromvcf'
start_pos = 0
stop_pos =  3_570_000
class_encoding_filename = None# '../image_data/sorted/class_encodings.xlsx'
input_image_directory = '../swarm_image_data/' # '../image_data/images_background_removed'

    
start = time.time()
generate_profile(
    file_name,
    gene_name,
    start_pos,
    stop_pos, 
    create_styled_table=False,
    class_encodings_filename=class_encoding_filename, 
    input_image_directory = input_image_directory, 
    create_image_plots=False

)
stop = time.time() - start

print(stop) 

ctherm_data = pd.read_csv('./ctherm_DSM1313_test_fromvcf_0_to_3570000_visualizations/ctherm_DSM1313_test_fromvcf_0_3570000_snp_matchings_by_group.tsv', sep='\t', index_col=0)
gene_start = 1
gene_stop = 3_840_000
# out = snp_data.apply(
#     lambda x: calculate_start_stops(
#         x, 
#         actual_gene_start, 
#         actual_gene_stop), 
#     axis=1)
no_counts = ctherm_data[ctherm_data.columns[:-1]]
out = no_counts.apply(
    lambda x: calculate_start_stops(
        x, 
        gene_start, 
        gene_stop), 
    axis=1)

output_df = pd.concat(out.values)
output_df['start'] = output_df['start'].apply(int)
output_df['end'] = output_df['end'].apply(int)
output_df.sort_values(by='start').to_csv("ctherm_DSM1313_vcf.tsv", index=False)



grouped = output_df[['indv','insertion_size']].groupby('indv').sum()
mean_size = grouped['insertion_size'].mean()
left_group = grouped[grouped['insertion_size'] < mean_size] 
right_group = grouped[grouped['insertion_size'] > mean_size]

left_df = output_df[output_df['indv'].isin(left_group.index)]
left_df.to_csv('ctherm_DSM1313_vcf_minor.tsv', index=False)
right_df = output_df[output_df['indv'].isin(right_group.index)]
right_df.to_csv('ctherm_DFM1313_vcf_major.tsv', index=False)
