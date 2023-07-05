from lib.class_analysis import extract_names
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.gridspec as gridspec
import numpy as np
import os 
import pandas as pd
from PIL import Image
from PIL import ImageChops
import pydash
import seaborn as sns


def create_representative_set_dataframe(element_map, class_encodings = None): 
	df = pd.DataFrame({'Representative': element_map.values(), 'Sample': list(element_map.keys())})
	
	if class_encodings == None:
		df['matches']  = df['Sample'].map(lambda x: '')
		return df

	df['matches'] = df['Sample'].apply(lambda x: (class_encodings[ class_encodings['name'] == x]).dropna(axis=1))
	df['matches'] = df['matches'].apply(lambda x: list(filter(lambda y: y != "Unnamed: 0" and y!= 'name', x)))                                             
	df['filename'] = df['Sample'].apply(lambda x: 'NA' if len(class_encodings[class_encodings['name'] == x]['Unnamed: 0']) == 0 else
															class_encodings[class_encodings['name'] == x]['Unnamed: 0'].values[0])

	return df

def get_class_per_representative(representative, desired_class, df): 
	matches_class_mask = df[df['Representative'] == representative ]['matches'].apply(lambda x: desired_class in x)
	return df[ df['Representative'] == representative ]['matches'][matches_class_mask]


def trim_image(image): 
	vals = np.array(image)
	ge_0 = np.argwhere(vals > 0)
	y_min = ge_0[:, 0].min()
	y_max = ge_0[:, 0].max()
	x_min = ge_0[:, 1].min()
	x_max = ge_0[:, 1].max()
	return Image.fromarray(vals[y_min:y_max, x_min:x_max, :])

def create_figure(n_rows, n_cols, df, gene_name, title, match, representative, input_image_directory='../image_data/images_background_removed/', verbose=False): 
	
	image_files = os.listdir(input_image_directory)

	actual_image_list = list(filter( lambda x: x in image_files, df['filename'].values))
	actual_image_count = len(actual_image_list)

	if verbose: 	
		print("Desired amount of images: \t\t\t\t", len(df))
		print('Available Matching Images in Directory: \t\t', len(actual_image_list))
		print("TOTAL ROWS:", n_rows, "\tTOTAL COLS:", n_cols)

	plt.ioff()

	fig2 = plt.figure(figsize=(int(n_cols*5), int(n_rows*5)), facecolor='black')
	spec2 = gridspec.GridSpec(ncols=n_cols, nrows=n_rows, figure=fig2)

	image_count = 0
	row_count = 0
	
	row = 0
	images_per_row = 0
	
	
	for i in range(actual_image_count): 
		try: 
			file = actual_image_list[i]
			
			matches = df[df['filename'] == file]['matches'].values[0]
			
			
			sub_title = f"{file}\n{matches}"

			col = image_count % n_cols
			
			if images_per_row == n_cols: 
				row = row + 1
				images_per_row = 0
			if verbose: 
				print(file, '\t', 'row', row, ':\t', 'col:', col)
			im = Image.open(f'{input_image_directory}/{file}')

			im = trim_image(im) 
			# ax.imshow(im, aspect='auto')

			f2_ax1 = fig2.add_subplot(spec2[row, col])
			f2_ax1.imshow(im)
			f2_ax1.axes.get_xaxis().set_visible(False)
			f2_ax1.axes.get_yaxis().set_visible(False)
			f2_ax1.set_title(sub_title, color='white')
			images_per_row = images_per_row + 1 
			
			
			image_count = image_count + 1
		except Exception as ex: 
			print("Exception: ", ex)
			# print()
		

	size = 32 if n_rows == 1 else 64
	plt.suptitle(title, color='white', size=size)
	filename = f'{gene_name}_{representative}_{match}_superplot.png'

	plt.savefig(filename)

	plt.ion()

	return filename

def save_count_plot(df, title, gene_name, match, representative): 
	# print(df, title, gene_name, match, representative)
	# print(pydash.flatten_deep(df[ df['Representative'] == representative ]['matches'].values))
	plt.figure(figsize=(10,10))
	sns.countplot(x=pydash.flatten_deep(df[ df['Representative'] == representative ]['matches'].values) )
	plt.xticks(rotation=90)
	plt.suptitle(title)
	plt.savefig(f'{gene_name}_{representative}_{match}_class_counts.png')
	
def print_images_from_set(representative, representative_set, class_encodings, gene_name, class_match, n_figures_per_row=5, input_image_directory='', verbose=False):
	dataframe = create_representative_set_dataframe(representative_set, class_encodings)
	df_slice = dataframe[dataframe['Representative'] == representative]

	if len(df_slice) == 0: 
		print("No reprsentatives match")
		return

	no_matches_class_mask = df_slice['matches'].apply(lambda x: class_match not in x)
	matches_class_mask = df_slice['matches'].apply(lambda x: class_match in x)

	if verbose: 
		print("Matches: ")
		print(df_slice[matches_class_mask])
		print(len(df_slice[matches_class_mask]))

	if len(df_slice[matches_class_mask]) != 0:
		n_rows = int(np.ceil(len(df_slice[matches_class_mask]) / n_figures_per_row))
		n_columns = n_figures_per_row
		
		save_count_plot(df_slice[matches_class_mask], f"{gene_name} {class_match} Matches", gene_name, "match", representative)
		create_figure(n_rows, n_columns, df_slice[matches_class_mask], gene_name, f"{gene_name} {class_match} Matches to rep {representative}", "match", representative, input_image_directory=input_image_directory, verbose=verbose)

	if verbose: 
		print("NO Matches: ")
		print(df_slice[no_matches_class_mask])
		print(len(df_slice[no_matches_class_mask]))

	if len(df_slice[no_matches_class_mask]) != 0: 
		n_rows = int(len(df_slice[no_matches_class_mask]) / n_figures_per_row)
		n_columns = n_figures_per_row

		save_count_plot(df_slice[no_matches_class_mask], f"{gene_name} {class_match} No Matches to rep {representative}", gene_name, "no_match", representative)

		create_figure(
			n_rows,
			n_columns,
			df_slice[no_matches_class_mask],
			gene_name,
			f"{gene_name} {class_match} No Matches to rep {representative}",
			"no_match",
			representative,
			input_image_directory=input_image_directory,
			verbose=verbose
		)
		
def plot_representative_set(representative_set, gene_name, representative, n_figures_per_row=5, input_image_directory='../image_data', file_list_path = None, verbose=False): 
	dataframe = create_representative_set_dataframe(representative_set)
	file_list_df = pd.read_excel(file_list_path, header=None)
	named = extract_names(file_list_df, 0)
	
	dataframe['filename'] = dataframe['Sample'].apply(lambda x: np.nan if len(named[ named['name'] == x]) == 0 else named[ named['name'] == x].values[0][0])

	df_slice = dataframe[dataframe['Representative'] == representative]
	df_slice = df_slice[ df_slice['filename'].notna() ]

	n_rows = int(np.ceil(len(df_slice) / n_figures_per_row))
	n_rows = n_rows if n_rows > 0 else 1
	n_columns = n_figures_per_row

	# save_count_plot(
	# 	df_slice,
	# 	f"{gene_name} grouped to {representative}",
	# 	gene_name,
	# 	"grouped_set",
	# 	representative
	# )

	outfile_name = create_figure(
		n_rows,
		n_columns,
		df_slice,
		gene_name,
		f"{gene_name} grouped to {representative}",
		"grouped_set",
		representative,
		input_image_directory=input_image_directory,
		verbose=verbose
	)

	return outfile_name