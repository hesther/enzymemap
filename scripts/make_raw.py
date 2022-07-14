from enzymemap import make_initial
import pandas as pd
import json

df, compound_to_smiles = make_initial(file_loc = '../data/brenda_2023_1.txt',
                                       file_loc_inchi = '../data/brenda_ligands.csv',
                                       file_loc_chebi = '../data/brenda_ligands_chebi.csv',
                                       manual_corrections=True)
df.to_csv("../data/raw_reactions.csv",index=False)
with open("../data/compound_to_smiles.json", 'w') as f:
    json.dump(compound_to_smiles,f)

with open("ec_nums.csv",'w') as f:
    for ec_num in sorted(list(set(df['EC_NUM']))):
        print(ec_num,file=f)
