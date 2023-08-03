import pandas as pd
import ast
from rdkit import Chem

def unmap_mol(smi):
    m=Chem.MolFromSmiles(smi)
    [a.SetAtomMapNum(0) for a in m.GetAtoms()]
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(m)))
def unmap_rxn(rxn):
    r,_,p=rxn.split(">")
    return unmap_mol(r)+">>"+unmap_mol(p)

#enzymemap brenda
df = pd.read_csv("../../../data/processed_reactions.csv")
rxns_enzymemap = list(set(list(df['unmapped'].values)))


#raw metacyc before mapping with enzymemap workflow (i.e. balanced and cofactors added):
df=pd.read_csv("../data/raw_reactions.csv")
rxns_metacyc_raw = [ast.literal_eval(x)[0] for x in df['BALANCED_RXNS']]


#processed metacyc after mapping with enzymemap workflow (i.e. only those that were resolved, plus reversed and suggested):
df=pd.read_csv("../data/processed_reactions.csv")
rxns_metacyc_processed = list(df['unmapped'].values)


ctr1=0
ctr2=0

for rxn in rxns_enzymemap:
    if rxn in rxns_metacyc_raw:
        ctr1+=1
    if rxn in rxns_metacyc_processed:
        ctr2+=1


print('rxns also in raw metacyc',ctr1)
print('rxns also in processed metacyc',ctr2)
