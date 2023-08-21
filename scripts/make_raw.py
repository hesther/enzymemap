from enzymemap import make_initial
import enzymemap
import json
import os

# Define the data directory
data_dir = os.path.join(os.path.dirname(enzymemap.__file__), 'data')  

# go through and parse the BRENDA flatfile, map compounds to smiles, balance reactions, and correct reactions
df, compound_to_smiles = make_initial(file_loc = os.path.join(data_dir, 'brenda_2023_1.txt'), 
                                       file_loc_inchi = os.path.join(data_dir, 'brenda_ligands.csv'),
                                       file_loc_chebi = os.path.join(data_dir, 'brenda_ligands_chebi.csv'),
                                       manual_corrections=True)

# Save the raw reactions data to a CSV file
df.to_csv(os.path.join(data_dir, "raw_reactions.csv"),index=False)  

# Open a new JSON file to store compound to SMILES mapping
with open(os.path.join(data_dir, "compound_to_smiles.json"), 'w') as f:  
    json.dump(compound_to_smiles,f)

# Open a new CSV file to store EC numbers
with open(os.path.join(data_dir, "ec_nums.csv"),'w') as f:  
    for ec_num in sorted(list(set(df['EC_NUM']))):  # Iterate over the unique EC numbers in the dataframe
        print(ec_num, file=f)  # Print each EC number to the CSV file
