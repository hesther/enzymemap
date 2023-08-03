import pandas as pd
import ast
from rdkit import Chem

#enzymemap brenda
df = pd.read_csv("../../../data/processed_reactions.csv")
rxns_enzymemap = list(set(list(df['unmapped'].values)))

#raw kegg before mapping with enzymemap workflow (i.e. balanced and cofactors added):
df=pd.read_csv("../data/raw_reactions.csv")
rxns_kegg_raw = [item for x in df['BALANCED_RXNS'] for item in ast.literal_eval(x)]
print(len(rxns_kegg_raw))


#processed kegg after mapping with enzymemap workflow (i.e. only those that were resolved, plus reversed and suggested):
df=pd.read_csv("../data/processed_reactions.csv")
rxns_kegg_processed = list(df['unmapped'].values)


ctr1=0
ctr2=0
for rxn in rxns_enzymemap:
    if rxn in rxns_kegg_raw:
        ctr1+=1
    if rxn in rxns_kegg_processed:
        ctr2+=1
print('rxns also in raw kegg',ctr1)
print('rxns also in processed kegg',ctr2)
