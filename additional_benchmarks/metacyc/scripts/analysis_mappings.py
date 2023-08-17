import pandas as pd
import enzymemap
from rdkit import Chem
from copy import deepcopy

def same_through_sim(rxn1,rxn2):
    #must be standardized, i.e. reactants the same
    if rxn1.split(">")[0]!=rxn2.split(">")[0]:
        return False
    reac = Chem.MolFromSmiles(rxn1.split(">")[0])
    reac_map = dict([(a.GetIdx(),a.GetAtomMapNum()) for a in reac.GetAtoms()]) 
    prod_smi1 = rxn1.split(">")[-1]
    prod_smi2 = rxn2.split(">")[-1]
    prod1 = Chem.MolFromSmiles(prod_smi1)
    prod2 = Chem.MolFromSmiles(prod_smi2)
    for match in reac.GetSubstructMatches(reac,uniquify=False):
        remap = dict([(reac_map[i],reac_map[j]) for i, j in enumerate(match)])
        prod_remap=deepcopy(prod2)
        [a.SetAtomMapNum(remap.get(a.GetAtomMapNum(),0)) for a in prod_remap.GetAtoms()]
        prod_remap_smi=Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(prod_remap)))
        if prod_remap_smi==prod_smi1:
            return True
    return False


df1=pd.read_csv("../data/processed_reactions.csv")
df2=pd.read_csv("../data/raw_reactions.csv")

for n in range(1,7):
    ctr_true=0
    ctr_false=0
    ctr_better=0
    ctr_worse=0
    for i in df1.index: 
        if df1['ec_num'][i][0]!=str(n):
            continue
        rxn_idx=df1['rxn_idx'][i]
        if df1['source'][i] == 'direct' and df1['steps'][i] != 'single from multi':
            try:
                m_m=enzymemap.helpers_map.standardize_mapping_rxn(df2['ORIG_RXNS'][rxn_idx])
            except:
                ctr_false+=1
                continue
            m_e=df1['mapped'][i]
            if m_m==m_e or same_through_sim(m_m,m_e):
                ctr_true+=1
            else:
                ctr_false+=1
                if sum(enzymemap.helpers_rdkit.bond_edit_stats(m_m).values()) >sum(enzymemap.helpers_rdkit.bond_edit_stats(m_e).values()):
                    ctr_better += 1
                elif sum(enzymemap.helpers_rdkit.bond_edit_stats(m_m).values()) < sum(enzymemap.helpers_rdkit.bond_edit_stats(m_e).values()):
                    ctr_worse += 1
            #print(i,m_m==m_e,sum(enzymemap.helpers_rdkit.bond_edit_stats(m_m).values()),sum(enzymemap.helpers_rdkit.bond_edit_stats(m_e).values()),df1['ec_num'][i])
    print("EC",n)
    print("same",ctr_true,ctr_true/(ctr_true+ctr_false)*100)
    print("diff",ctr_false,ctr_false/(ctr_true+ctr_false)*100)
    print("diff better",ctr_better,ctr_better/(ctr_true+ctr_false)*100)
    print("diff better or same",100-ctr_worse/(ctr_true+ctr_false)*100)
    print("diff worse",ctr_worse,ctr_worse/(ctr_true+ctr_false)*100)
