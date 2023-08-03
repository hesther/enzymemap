import pandas as pd
import json
from enzymemap import helpers_brenda, helpers_resolve_smiles, helpers_rdkit, helpers_map
from enzymemap.helpers_brenda import extract_reaction
import itertools

def process_entry(ec_num,buff):
    reactions = []
    if not '-' in ec_num:
        tmp = pd.DataFrame([x.split("  ") for x in buff.split('\n')[:-1]])[[4,5]]
        for i in tmp.index:
            line = tmp[5][i].replace('<=>','=')
            desc = extract_reaction(line)
            desc['EC_NUM']=ec_num
            desc['ID']=tmp[4][i]
            desc['ORIG_RXN_TEXT']=desc['RXN_TEXT']
            reactions.append(desc)
    return reactions
    
def parse_kegg(file_loc):
    """parse_kegg.

    Parses a KEGG text file

    Args:
        file_loc (str): Location of KEGG download text file

    Return: 
        Pandas dataframe containing EC numbers, lists of substrates and products, the reaction texts and the reversibility tag, as well as a dataframe of compounds
    """
    df = pd.DataFrame(columns=['EC_NUM','SUBSTRATES','PRODUCTS','RXN_TEXT','REVERSIBLE','ORIG_RXN_TEXT','ID'])

    buff = ''
    ec_num = ''
    with open(file_loc) as fp:
        new_line = fp.readline()
        while (new_line):
            if new_line[0]=="D":
                if buff != '':
                    reactions = process_entry(ec_num,buff)
                    df = pd.concat([df,pd.DataFrame(reactions)],ignore_index=True)
                ec_num = new_line.strip().split()[-1]
                buff = ''
            elif new_line[0]=='E':
                buff += new_line

            new_line = fp.readline()
            if not new_line:
                reactions = process_entry(ec_num,buff)
                df = pd.concat([df,pd.DataFrame(reactions)],ignore_index=True)
                
    df = df.loc[df.astype(str).drop_duplicates().index] # Drop duplicates (cast as str to avoid problem with list)
    df = df[~df['EC_NUM'].str.startswith('7')] # Drop transferases
    df = df.reset_index(drop=True) #Reset index

    print("Read in ")
    print(df)

    compounds=set()
    for x in df['SUBSTRATES']:
        for xx in x:
            compounds.add(xx)
    for x in df['PRODUCTS']:
        for xx in x:
            compounds.add(xx)

    compound_df = pd.DataFrame(compounds,columns=['compound'])

    return df, compound_df

def make_initial(file_loc, file_loc_inchi, file_loc_chebi, manual_corrections=False):
    """make_initial.

    Function to read BRENDA database and resolve names to SMILES.

    Args:
        file_loc (str): Location of BRENDA download text file
        file_loc_inchi (str): Location of BRENDA download of inchi keys
        file_loc_chebi (str): Location of BRENDA download of chebi keys
        manual_corrections (bool): Whether to add a set of manual corrections for SMILES entries. Do not use for custom datasets.

    Returns:
        Pandas Dataframes of uncorrected, raw reactions and compound to smiles dictionary
    """

    df, compounds_df = parse_kegg(file_loc)

    compound_df = helpers_resolve_smiles.resolve_all(compounds_df, file_loc_inchi, file_loc_chebi)
    compound_df = helpers_resolve_smiles.standardize_compound_df(compound_df)
    if manual_corrections:
        compound_df = helpers_resolve_smiles.manual_corrections_compounds(compound_df)

    compound_to_smiles={}
    for i in compound_df.index:
        compound_to_smiles[compound_df['compound'][i]] = compound_df['smiles_neutral'][i]
    
    reactions = []
    for i in df.index:
        #Make all possible reactions from different compound_to_smiles results
        reacs = [".".join(x) for x in list(itertools.product(*[compound_to_smiles[x] for x in df['SUBSTRATES'][i]]))]
        prods = [".".join(x) for x in list(itertools.product(*[compound_to_smiles[x] for x in df['PRODUCTS'][i]]))]
        reaction = [">>".join(x) for x in list(itertools.product(*[reacs, prods]))]
        reaction = [helpers_rdkit.put_h_last(r) for r in reaction]
        reactions.append(reaction)
    df['POSSIBLE_RXNS'] = reactions

    reduced_to_oxidized, oxidized_to_reduced = helpers_rdkit.get_strip_list(compound_to_smiles)
    
    # Correct and balance reactions
    cached=dict()
    print("Correcting and balancing")
    balanced_rxn_list = []
    for i in df.index:
        print(i,end='\r')
        if df['RXN_TEXT'][i] not in cached:
            rxns = helpers_rdkit.correct_reaction(df['POSSIBLE_RXNS'][i], df['RXN_TEXT'][i], reduced_to_oxidized, oxidized_to_reduced)
            rxns = [helpers_rdkit.put_h_last(r) for r in rxns]
            cached[df['RXN_TEXT'][i]]=rxns
        else:
            rxns = cached[df['RXN_TEXT'][i]]
        balanced_rxn_list.append(rxns)
    df['BALANCED_RXNS'] = balanced_rxn_list

    print("not resolvable:", sum([len(x) == 0 for x in df['POSSIBLE_RXNS']]),"out of",len(df),"entries")
    print("resolvable but unbalanced:", sum([len(x) == 0 for x in df['BALANCED_RXNS']])-sum([len(x) == 0 for x in df['POSSIBLE_RXNS']]),"out of",len(df),"entries")
    
    return df, compound_to_smiles


df, compound_to_smiles = make_initial(file_loc = '../data/kegg_reactions_2023_07_24.keg',
                                       file_loc_inchi = '../../../data/brenda_ligands.csv',
                                       file_loc_chebi = '../../../data/brenda_ligands_chebi.csv',
                                       manual_corrections=False)
df['PROTEIN_REFS'] = [[]]*len(df)
df['NATURAL']=None
df['ORGANISM']=None
df['PROTEIN_DB']=None

df.to_csv("../data/raw_reactions.csv",index=False)
with open("../data/compound_to_smiles.json", 'w') as f:
    json.dump(compound_to_smiles,f)

with open("ec_nums.csv",'w') as f:
    for ec_num in sorted(list(set(df['EC_NUM']))):
        print(ec_num,file=f)
