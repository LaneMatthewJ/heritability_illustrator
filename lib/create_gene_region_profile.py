import os
import shutil
import pandas as pd
from lib.load_and_extract_data import load_and_extract_data, analyze_output, get_non_dominant_strains
from lib.plot_images import plot_representative_set


def generate_profile(file_name: str, gene_name: str, start: int, stop: int, create_styled_table=False, class_encodings_filename=None, input_image_directory='../image_data/images_background_removed'):
	"""
	This function takes a particular gene of interest, calculates its region from the hapmap files, and then 
	creates summary groupings and visualizes the output SNP regions. 
	"""
	outfile_list = []
	print("Generating Data")
	# Create groupints / visualization (if specified)
	visualization, dataframe, element_map, outfile_name = load_and_extract_data(
		filename = file_name, 
		gene_name = gene_name, 
		start = start,
		stop = stop, 
		create_styled_table=create_styled_table)

	# Append file name to list
	if outfile_name is not None: 
		outfile_list.append(outfile_name)

	print("Saving Group Statistics")
	# Save counts / snp groupings to file. 
	dataframe_filename = f"{gene_name}_{start}_{stop}_snp_matchings_by_group.tsv"
	dataframe.to_csv(dataframe_filename, sep='\t', index=True)

	outfile_list.append(dataframe_filename)

	# Create element 
	element_map_df = pd.DataFrame({"Gene": element_map.keys(), "Group Value": element_map.values()})
	element_map_filename = f"{gene_name}_{start}_{stop}_groupings.tsv"
	element_map_df.to_csv(element_map_filename, sep='\t', index=False)

	outfile_list.append(element_map_filename)

	print("Creating Non Dominant Plots")
	non_dominant_strains = get_non_dominant_strains(dataframe, max_by='row')
	for non_dominant_representative in non_dominant_strains: 

		
		figs_per_row = 5
		
		n_counts = int(dataframe.loc[non_dominant_representative]['counts'])
		
		if n_counts < figs_per_row:
			figs_per_row = n_counts
		
		outfile = plot_representative_set(
			element_map,
			gene_name,
			non_dominant_representative,
			n_figures_per_row = figs_per_row,
			input_image_directory = input_image_directory,
			file_list_path='../image_data/sorted/file_list.xlsx')
		
		outfile_list.append(outfile)


	## Analyze Output: 
	if class_encodings_filename is not None: 
		class_base = class_encodings_filename.split('.')[0].split('/')[-1]
		analysis_file_out = f"{gene_name}_{gene_name}_{start}_to_{stop}_{class_base}.tsv"
		
		analyze_output(class_encodings_filename, analysis_file_out, element_map)

		outfile_list.append(analysis_file_out)

	gene_dir = f"{gene_name}_{start}_to_{stop}_visualizations"
    
	os.mkdir(gene_dir)

	for file in outfile_list: 
		shutil.move(file, f"{gene_dir}/{file}")
