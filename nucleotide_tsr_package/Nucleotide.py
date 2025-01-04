from .analyzer import AminoAcidAnalyzer
import os
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing

def NucleotideTSR(data_dir, input_files, chain=None, output_option='both', output_subdir='nucleotide_results', mirror_image=False):
        dtheta = 29
        dLen = 18
        numOfLabels = 112
        analyzer = AminoAcidAnalyzer(dtheta, dLen, numOfLabels)    
        os.makedirs(os.path.join(data_dir, output_subdir), exist_ok=True)
        analyzer.readDrugLexicalCsv("drug_atom_lexical_txt.csv")
        chain_dict = {}
        # Handle single file input
        if isinstance(input_files, str):
            if input_files.endswith('.csv'):
                # Read CSV
                df = pd.read_csv(input_files)
                chain_dict = dict(zip(df['protein'].str.upper(), df['chain']))
                input_files = df['protein'].str.upper().tolist()
            else:
                input_files = [input_files]
        if chain:
            if isinstance(chain, list):
                chain_dict = {f.upper(): c for f, c in zip(input_files, chain)}
            elif isinstance(chain, str):
                chain_dict = {f.upper(): chain for f in input_files}
        numOfCores = multiprocessing.cpu_count()

        def generate_keys_and_triplets(data_dir, file_name, chain, output_subdir, output_option, mirror_image):
            analyzer.readSeqAndIdentityChain(data_dir, file_name, chain)
            for seq_value, chain_identity in analyzer.seqchainIdentity.items():
                analyzer.calcuTheteAndKey(data_dir, file_name, chain, seq_value, chain_identity, output_subdir, output_option, mirror_image)

        # Using Parallel to process files concurrently
        Parallel(n_jobs=numOfCores, verbose=50)(
            delayed(generate_keys_and_triplets)(data_dir, file_name.upper(), chain_dict.get(file_name.upper(), chain), output_subdir, output_option, mirror_image)
            for file_name in input_files
        )
       

    # data_dir = 'Dataset/' 
    # input_files = ["2R93"]
    # chain = ["R"]
    # output_option = "both"
    # NucleotideTSR(data_dir, input_files, chain=chain, output_option=output_option, mirror_image=True)