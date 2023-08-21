from enzymemap import make_initial
import enzymemap
import json
import os

data_dir = os.path.join(os.path.dirname(enzymemap.__file__), 'data')  # Define the data directory

df, compound_to_smiles = make_initial(file_loc = os.path.join(data_dir, 'brenda_2023_1.txt'),  # Load BRENDA data
                                       file_loc_inchi = os.path.join(data_dir, 'brenda_ligands.csv'),  # Load InChI data
                                       file_loc_chebi = os.path.join(data_dir, 'brenda_ligands_chebi.csv'),  # Load ChEBI data
                                       manual_corrections=True)  # Apply manual corrections if necessary

df.to_csv(os.path.join(data_dir, "raw_reactions.csv"),index=False)  # Save the raw reactions data to a CSV file

with open(os.path.join(data_dir, "compound_to_smiles.json"), 'w') as f:  # Open a new JSON file to store compound to SMILES mapping
    json.dump(compound_to_smiles,f)  # Dump the compound to SMILES mapping into the JSON file

with open(os.path.join(data_dir, "ec_nums.csv"),'w') as f:  # Open a new CSV file to store EC numbers
    for ec_num in sorted(list(set(df['EC_NUM']))):  # Iterate over the unique EC numbers in the dataframe
        print(ec_num,file=f)  # Print each EC number to the CSV file
