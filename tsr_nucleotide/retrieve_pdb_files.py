import urllib.request
import multiprocessing
from joblib import Parallel, delayed
import pandas as pd
import os

# Function to download a single PDB file
def download_pdb(fname, out_dir):
    print(f"Downloading {fname}...")
    try:
        urllib.request.urlretrieve(f'http://files.rcsb.org/download/{fname.upper()}.pdb', os.path.join(out_dir, f'{fname.upper()}.pdb'))
        print(f"Downloaded: {fname}")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"Skipping {fname} - Protein file not found.")
        else:
            raise

# Function to handle the parallel downloading of PDB files
def download_pdb_files(protein_list, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores, verbose=50)(delayed(download_pdb)(fname, out_dir) for fname in protein_list)
    print("All downloads completed.")

# Function to read PDB IDs from a file
def read_pdb_from_file(filename):
    df = pd.read_csv(filename)
    return df['protein'].tolist()

# Main function that takes a single PDB ID, list of IDs, or filename with IDs
def retrieve_pdb_files(input_data, out_dir='Dataset/'):
    if isinstance(input_data, str):
        if os.path.isfile(input_data):
            # If input_data is a filename, read the list from the file
            print(f"Reading PDB IDs from file: {input_data}")
            protein_list = read_pdb_from_file(input_data)
        else:
            # If input_data is a single PDB ID
            print(f"Single PDB ID provided: {input_data}")
            protein_list = [input_data]
    elif isinstance(input_data, list):
        # If input_data is a list of PDB IDs
        print(f"List of PDB IDs provided: {input_data}")
        protein_list = input_data
    else:
        raise ValueError("Input data must be a single PDB ID, a list of IDs, or a file path.")

    # Download the PDB files
    download_pdb_files(protein_list, out_dir)