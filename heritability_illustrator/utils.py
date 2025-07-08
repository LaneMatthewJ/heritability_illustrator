import numpy as np

def extract_file_from_directory_list(sample_name, file_list):
    files = list(filter(lambda x: sample_name in x, file_list))

    if len(files) > 2: 
        print(files)

    return np.NaN if len(files) == 0 else files
    