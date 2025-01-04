import os
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import urllib.request


# Function to download a single PDB file
def download_pdb(fname, out_dir, timeout=10):
    pdb_url = f'https://files.rcsb.org/download/{fname.upper()}.pdb'
    save_path = os.path.join(out_dir, f'{fname.upper()}.pdb')
    
    # Check if the file already exists
    if os.path.exists(save_path):
        print(f"Skipping {fname} - File already exists.")
        return
    
    print(f"Downloading {fname}...")
    
    try:
        urllib.request.urlretrieve(
            f'http://files.rcsb.org/download/{fname.upper()}.pdb',
            os.path.join(out_dir, f'{fname.upper()}.pdb')
        )
        print(f"Downloaded: {fname}")
    
    except urllib.error.HTTPError as http_err:
        if http_err.code == 404:
            print(f"Skipping {fname} - Protein file not found (404).")
        else:
            print(f"HTTP error occurred while downloading {fname}: {http_err}")
    except urllib.error.URLError as url_err:
        print(f"Connection error occurred while downloading {fname}: {url_err}")
    except Exception as e:
        print(f"An unexpected error occurred while downloading {fname}: {e}")


# Function to handle the parallel downloading of PDB files
def download_pdb_files(protein_list, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores, verbose=0)(
        delayed(download_pdb)(fname, out_dir) for fname in protein_list
    )
    print("All downloads completed.")

# Function to read PDB IDs from a file
def read_pdb_from_file(filename):
    try:
        df = pd.read_csv(filename)
        if 'protein' not in df.columns:
            raise ValueError("CSV file must contain a 'protein' column.")
        return df['protein'].tolist()
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return []

# Main function that takes a single PDB ID, list of IDs, or filename with IDs
def PDB_DL(input_data, out_dir='Dataset/'):
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


'''
Example usage:
input_files = ["1GTA", "1gtb", "1lbe"]
PDB_DL(input_files)

OR

csv_file = "sample_details.csv"
PDB_DL(csv_file)
'''