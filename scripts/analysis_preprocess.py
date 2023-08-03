import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from joblib import Parallel, delayed
from templatecorr import canonicalize_mol, get_templates, switch_direction
from templatecorr import switch_direction, canonicalize_template, correct_templates, split_data_df
import random
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants
import os

def downsample(d):
    targets=[1000,3000,5000,10000]
    indices = list(set(d['index'].values))
    random.shuffle(indices)

    ds = []
    for i in range(len(indices)-1):
        d_i = d[d['index'].isin(indices[:i])]
        d_i1 = d[d['index'].isin(indices[:i+1])]
        for t in targets:
            if len(d_i)< t and len(d_i1)>=t:
                ds.append((t,d_i1))
    return ds

def remove_unmapped(rxn_smi):
    new_reacs = []
    reacs = [Chem.MolFromSmiles(x) for x in rxn_smi.split(">>")[0].split('.')]
    prod = Chem.MolFromSmiles(rxn_smi.split(">>")[-1]) #Must be single product
    prod_mapnos = [atom.GetAtomMapNum() for atom in prod.GetAtoms()]
    for reac in reacs:
        spectator = True
        for atom in reac.GetAtoms():
            mapno = atom.GetAtomMapNum()
            if mapno not in prod_mapnos:
                atom.SetAtomMapNum(0)
            else:
                spectator = False
        if not spectator:
            new_reacs.append(Chem.MolToSmiles(reac))
    return ".".join(new_reacs)+">>"+rxn_smi.split(">>")[-1]

def make_singleprod(path,column='mapped',system='brenda',n_cpus=48):
    '''
    Makes reactions with a single product.
    '''

    if system=='brendadirectsingle':
        data=pd.read_csv(path)
        data=data[data['source']=='direct']
        data=data[data['steps']=='single']
        rxns = data[column]
    elif system=='brendadirectall':
        data=pd.read_csv(path)
        data=data[data['source']=='direct']
        rxns = data[column]
    elif system=='brendanotsuggestedall':
        data=pd.read_csv(path)
        data=data[~data['source'].str.contains('suggest')]
        rxns = data[column]
    elif system=='brendanotreversedall':
        data=pd.read_csv(path)
        data==data[~data['source'].str.contains('reverse')]
        rxns = data[column]
    elif system=='brendaallsingle':
        data=pd.read_csv(path)
        data=data[data['steps']=='single']
        rxns = data[column]
    else:
        rxns = pd.read_csv(path)[column]

    df = pd.DataFrame(list(zip(rxns)), columns=['rxn_smiles'])
    df = df.drop_duplicates(subset=['rxn_smiles'])
    print("Read in",len(df),"reactions (without duplicates)")
    
    
    #Unmap, count occurence, remove everything that occurs more than 100 times to get rid of e.g. cofactors
    print("Count occurences")
    unmapped_prods = []
    unmapped_reacs = []
    for i,r in enumerate(df['rxn_smiles']):
        print(i,end='\r')
        xs = []
        for x in r.split(">")[0].split("."):
            mol = Chem.MolFromSmiles(x)
            if mol:
                [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
                smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
                xs.append(smi)
            else:
                xs.append(None)
        unmapped_reacs.append(xs)
    
        xs = []
        for x in r.split(">")[-1].split("."):
            mol = Chem.MolFromSmiles(x)
            if mol:
                [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
                smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
                xs.append(smi)
            else:
                xs.append(None)
        unmapped_prods.append(xs)
    
    df['u_r'] = unmapped_reacs
    df['u_p'] = unmapped_prods

    flat_list = [item for sublist in unmapped_reacs+unmapped_prods for item in sublist]
    unique = list(set(flat_list))
    d={}
    for i,x in enumerate(unique):
        d[x] = flat_list.count(x)
    remove = []
    for x in unique:
        if d[x] >= 100:
            remove.append(x)
    ctr=0
    new_rxn_smiles = []
    for i in df.index:
        ctr2=0
        s = df['rxn_smiles'][i].split(">>")[0] + ">>"
        for p,mp in zip(df['u_p'][i],df['rxn_smiles'][i].split(">>")[-1].split(".")):
            if p not in remove:
                ctr2+=1
                s += mp
        if ctr2==1: 
            try:
                s = remove_unmapped(s)
                new_rxn_smiles.append(s)
                ctr+= 1
            except:
                new_rxn_smiles.append(None)
        else:
            new_rxn_smiles.append(None)
    df['rxn_smiles'] = new_rxn_smiles
    df = df[df['rxn_smiles'].notnull()]
    print("retained",len(df),"reactions")
    
    print("Drop duplicates")
    #Drop duplicates (same reactions but different mappings)
    df["reac_smiles"] = Parallel(n_jobs=n_cpus, verbose=1)(delayed(canonicalize_mol)(rxn_smi,0) for rxn_smi in df["rxn_smiles"])
    df["prod_smiles"] = Parallel(n_jobs=n_cpus, verbose=1)(delayed(canonicalize_mol)(rxn_smi,-1) for rxn_smi in df["rxn_smiles"])
    df = df.drop_duplicates(subset=['prod_smiles', 'reac_smiles'])
    print("retained",len(df),"reactions after dropping duplicates")
    df[['rxn_smiles','reac_smiles','prod_smiles']].to_csv("../data/"+system+"_singleprod.csv", index=False)


def remove_unmatched(rxn_smi):
    reac=Chem.MolFromSmiles(rxn_smi.split(">")[0])
    prod=Chem.MolFromSmiles(rxn_smi.split(">")[-1])
    prod_map_list=[a.GetAtomMapNum() for a in prod.GetAtoms()]
    for a in reac.GetAtoms():
        if a.GetAtomMapNum() not in prod_map_list:
            a.SetAtomMapNum(0)
        
    return Chem.MolToSmiles(reac)+">>"+Chem.MolToSmiles(prod)

def mapped_rxn_smi(reac_smi,template,remove_lgs=False):
    rxn = rdchiralReaction(template)
    rct = rdchiralReactants(reac_smi)
    reac_mapped=Chem.MolToSmiles(rct.reactants)
    prods_mapped = rdchiralRun(rxn, rct, keep_mapnums=True, combine_enantiomers=False)
    rxn_smiles_list=[]
    for prod_mapped in prods_mapped:
        prod_unmapped = unmap(prod_mapped)
        rxn_smi=reac_mapped+">>"+prod_mapped
        if remove_lgs:
            rxn_smi=remove_unmatched(rxn_smi)
        rxn_smiles_list.append((rxn_smi,prod_unmapped))
    return rxn_smiles_list

def unmap(smi):
    mol = Chem.MolFromSmiles(smi)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    return smi

def read_and_save(system='brenda',n_cpus=48):
    """
    Processes the single-product reactions to compute fingerprints and templates
    
    :param system: String of system name.
    """

    getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),2,useChirality=True, useFeatures=True)
    template_choices = ["default", "r1", "r0"]
    data = pd.read_csv("../data/"+system+"_singleprod.csv")
    
    #Get templates
    print("Calculate templates")
    data["template_default"] = Parallel(n_jobs=n_cpus, verbose=1)(delayed(get_templates)(data["rxn_smiles"][idx], data["reac_smiles"][idx], False, 1) for idx in data.index)
    data["template_r1"] = Parallel(n_jobs=n_cpus, verbose=1)(delayed(get_templates)(data["rxn_smiles"][idx], data["reac_smiles"][idx], True, 1) for idx in data.index)
    data["template_r0"] = Parallel(n_jobs=n_cpus, verbose=1)(delayed(get_templates)(data["rxn_smiles"][idx], data["reac_smiles"][idx], True, 0) for idx in data.index)

    #Delete all lines without a template:
    data = data.dropna(subset=["template_default","template_r0","template_r1"])
    print(len(data))

    #Only keep necessary columns:
    template_columns=["template_"+choice for choice in template_choices] 
    save_columns=["rxn_smiles", "prod_smiles", "reac_smiles"] + template_columns
    data=data[save_columns]
    if not os.path.exists('data_analysis'):
            os.makedirs('data_analysis')
    data.to_csv("data_analysis/"+system+".csv")

def correct_loop(df,column_name2, template):
    """
    Calls correct_templates function for a set of templates where data[column_name1]==template

    :param df: Pandas dataframe.
    :param column_name2: Name of column with more specific templates
    :param template: Template

    :return: Indices of dataframe, corrected templates
    """
    templates=correct_templates(df[column_name2])
    return df.index, templates
    
def correct_all_templates(data,column_name1,column_name2,n_cpus=48):
    """
    Computes a corrected set of templates for templates of different specificity levels 

    :param data: Pandas dataframe
    :param column_name1: Name of column with more general templates
    :param column_name2: Name of column with more specific templates
    :return: List of new templates in order of templates in dataframe
    """
    unique_templates=sorted(list(set(data[column_name1].values)))
    large_unique_templates=sorted(list(set(data[column_name2].values)))
    data["new_t"] = None
    print("...Unique templates in column",column_name1,":",len(unique_templates))
    print("...Unique templates in column",column_name2,":",len(large_unique_templates))
    print("...Correcting templates in column",column_name2)

    results = Parallel(n_jobs=n_cpus, verbose=1)(delayed(correct_loop)(data[data[column_name1]==template].copy(), column_name2, template) for template in unique_templates)
    for result in results:
        idxs, templates = result
        ctr=0
        for idx in idxs:
            data.at[idx,"new_t"]=templates[ctr]
            ctr+=1  
    new_unique_templates=set(data["new_t"].values)
    print("...Unique corrected templates in column",column_name2,":",len(new_unique_templates))
    print("")
    return list(data["new_t"].values)

def count_reacs_per_template(data,column):
    """
    Computes the number of reactions associated with a template
    
    :param data: Pandas dataframe
    :param column: Name of column with templates
    """

    templates=list(data[column].values)
    unique_templates=sorted(list(set(templates)))

    reactions_per_template=[]
    for template in unique_templates:
        reactions_per_template.append(len(data[data[column]==template]))
    counts=np.bincount(reactions_per_template)

    count_1=counts[1]/len(reactions_per_template)
    count_2_5=sum(counts[2:6])/len(reactions_per_template)
    count_6_10=sum(counts[6:11])/len(reactions_per_template)
    count_more=sum(counts[11:]) /len(reactions_per_template)
    print("%7s %7s %7s %7s" % ("1","2-5","6-10",">10"))
    print("%7.3f %7.3f %7.3f %7.3f" % (count_1,count_2_5,count_6_10,count_more))

def number_unique_templates(data,column,n):
    """
    Computes the number of unique (based on string) templates in n reactions
    
    :param data: Pandas dataframe
    :param column: Name of column with templates
    :param n: Number of reactions

    :return: Number of unique templates
    """    
    data_subset = data.sample(n = n)
    templates=list(data_subset[column].values)
    unique_templates=set(templates)

    return len(unique_templates)

def process(system,n_cpus=48):
    """
    Processes the data for either uspto_50k or uspto_460k
    
    :param system: String of system name.
    """
    
    template_choices = ["default", "r1", "r0"]
    
    data=pd.read_csv("data_analysis/"+system+".csv")


    #Correct regular templates
    print("Correct noncanonical templates")
    data["corrected_template_r1"] = correct_all_templates(data,"template_r0","template_r1",n_cpus=n_cpus)
    data["corrected_template_default"] = correct_all_templates(data,"corrected_template_r1","template_default",n_cpus=n_cpus)

    #Only keep necessary columns:
    template_columns=["corrected_template_default"]

    save_columns=["rxn_smiles",
                  "prod_smiles",
                  "reac_smiles"]+template_columns
    data=data[save_columns]

    # Find unique template indices
    print("Preprocess data for ML-fixed algorithm")
    lengths={}
    for column in template_columns:
        unique_templates=sorted(list(set(data[column].values)))
        with open("data_analysis/"+system+"_"+column+"_unique_templates.txt", "w") as f:
            for item in unique_templates:
                f.write("%s\n" % item)
        lengths[column]=len(unique_templates)
        template_ids=[]
        for template in data[column]:
            template_ids.append(unique_templates.index(template))
        data[column+"_id"]=template_ids
    
    #Split data to train/val/test
    print("Splitting dataset")
    data = split_data_df(data)
    data.to_csv("data_analysis/"+system+"_all.csv",index=False)

    datasub = data.loc[data["dataset"] == "train"]
    datasub_val = data.loc [data["dataset"] == "val"]
    datasub_test = data.loc [data["dataset"] == "test"]


    #Save for ML
    for column in template_columns:
        np.save("data_analysis/"+system+"num_classes_"+column+".npy",lengths[column])

    for column in template_columns:
        if "forward" in column:
            datasub[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_train.csv",index=False)
            datasub_val[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_val.csv",index=False)
            datasub_test[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_test.csv",index=False)
        else:
            datasub[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_train.csv",index=False)
            datasub_val[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_val.csv",index=False)
            datasub_test[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_test.csv",index=False)


    #Preprocess for CGR
    print("preprocess for cgr-chemprop")
    X=[]
    index=[]
    y=[]

    for i in data.index:
        print(i,end="\r")
        try:
            outcomes = mapped_rxn_smi(data['reac_smiles'][i],switch_direction(data['corrected_template_r1'][i]),remove_lgs=True)
        except:
            outcomes=[]
        if len(outcomes) <= 1:
            continue
        if len(set([outcome[1] for outcome in outcomes])) <= 1:
            continue
        if not True in [outcome[1] == data['prod_smiles'][i] for outcome in outcomes]:
            print("no result for", i)
            continue

        taken = []
        for outcome in outcomes:
            if outcome[1] not in taken:
                X.append(outcome[0])
                index.append(i)
                if outcome[1] == data['prod_smiles'][i]:
                    y.append(1)
                else:
                    y.append(0)
                taken.append(outcome[1])

    df3 = pd.DataFrame(list(zip(X, y, index)),columns=['smiles','class','index'])
    indices = list(set(df3['index'].values))
    for i in indices:
        if sum(df3[df3['index']==i]['class'].values) != 1:
            print("WARNING, wrong number of correct products",i)
        
    df3.to_csv("data_analysis/"+system+"_cgr_class.csv",index=False)

    df3=pd.read_csv("data_analysis/"+system+"_cgr_class.csv")
    
    for ii in range(5,10):
        indices = list(set(df3['index'].values))
        random.shuffle(indices)
        df3['set'] = None
        for i in df3.index:
            if df3['index'][i] in indices[:int(0.8*len((indices)))]:
                df3['set'][i]='train'
            elif df3['index'][i] in indices[int(0.9*len((indices))):]:
                df3['set'][i]='test'
            else:
                df3['set'][i]='val'

        df3[df3['set']=='train'].to_csv("data_analysis/"+system+"_cgr_class_train_split"+str(ii)+".csv",index=False)
        df3[df3['set']=='val'].to_csv("data_analysis/"+system+"_cgr_class_val_split"+str(ii)+".csv",index=False)
        df3[df3['set']=='test'].to_csv("data_analysis/"+system+"_cgr_class_test_split"+str(ii)+".csv",index=False)
    
        print("Downsample")

        d = pd.read_csv("data_analysis/"+system+"_cgr_class_train_split"+str(ii)+".csv")
        ds = downsample(d)
        for t,d in ds:
            print(system, t)
            d.to_csv("data_analysis/"+system+"_cgr_class_train_split"+str(ii)+"_"+str(t)+".csv",index=False)


def process_same_test(system,test,n_cpus=48):
    """
    Processes the data for either uspto_50k or uspto_460k
    
    :param system: String of system name.
    :param test: String of test file
    """
    
    template_choices = ["default", "r1", "r0"]
    template_columns=["corrected_template_default"]
    
    data=pd.read_csv("data_analysis/"+system+"_all.csv")
    data_test=pd.read_csv("data_analysis/"+test+"_all.csv")
    data_test=data_test.loc[data_test["dataset"] == "test"]
    list_test=list(data_test['rxn_smiles'].values)
    print("External test",len(data_test))


    #Drop rows which are in test
    print(len(data))
    data=data.drop(data[data['rxn_smiles'].isin(list_test)].index)
    print(len(data))
 
    #Split data to train/val/test
    print("Splitting dataset")
    data = split_data_df(data,val_frac=0.2,test_frac=0.0)

    data=pd.concat([data,data_test],ignore_index=True)


    
    # Find unique template indices
    print("Preprocess data for ML-fixed algorithm")
    lengths={}
    for column in template_columns:
        unique_templates=sorted(list(set(data[column].values)))
        with open("data_analysis/"+system+"_"+column+"_unique_templates.txt", "w") as f:
            for item in unique_templates:
                f.write("%s\n" % item)
        lengths[column]=len(unique_templates)
        template_ids=[]
        for template in data[column]:
            template_ids.append(unique_templates.index(template))
        data[column+"_id"]=template_ids
    

    datasub = data.loc[data["dataset"] == "train"]
    datasub_val = data.loc [data["dataset"] == "val"]
    datasub_test = data.loc [data["dataset"] == "test"]


    system+='_sametest'
    #Save for ML
    for column in template_columns:
        np.save("data_analysis/"+system+"num_classes_"+column+".npy",lengths[column])

    for column in template_columns:
        if "forward" in column:
            datasub[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_train.csv",index=False)
            datasub_val[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_val.csv",index=False)
            datasub_test[["reac_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_test.csv",index=False)
        else:
            datasub[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_train.csv",index=False)
            datasub_val[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_val.csv",index=False)
            datasub_test[["prod_smiles",column+"_id"]].to_csv("data_analysis/"+system+"_"+column+"_test.csv",index=False)

            
if __name__ == '__main__':
    #enzymemap
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brenda',n_cpus=48)
    read_and_save(system='brenda',n_cpus=48)
    process(system='brenda',n_cpus=48)
    
    #enzymemap direct single                                                                                                                                                                                                              
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brendadirectsingle',n_cpus=48)
    read_and_save(system='brendadirectsingle',n_cpus=48)
    process(system='brendadirectsingle',n_cpus=48)
    process_same_test(system='brendadirectsingle',test='brenda',n_cpus=48)

    #enzymemap direct all                                                                                                                                                                                                              
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brendadirectall',n_cpus=48)
    read_and_save(system='brendadirectall',n_cpus=48)
    process(system='brendadirectall',n_cpus=48)
    process_same_test(system='brendadirectall',test='brenda',n_cpus=48)

    #enzymemap not suggested all                                                                                                                                                                                                              
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brendanotsuggestedall',n_cpus=48)
    read_and_save(system='brendanotsuggestedall',n_cpus=48)
    process(system='brendanotsuggestedall',n_cpus=48)
    process_same_test(system='brendanotsuggestedall',test='brenda',n_cpus=48)

    #enzymemap not reversed all                                                                                                                                                                                                              
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brendanotreversedall',n_cpus=48)
    read_and_save(system='brendanotreversedall',n_cpus=48)
    process(system='brendanotreversedall',n_cpus=48)
    process_same_test(system='brendanotreversedall',test='brenda',n_cpus=48)

    #enzymemap all single                                                                                                                                                                                                           
    make_singleprod(path="../data/processed_reactions.csv",column='mapped',system='brendaallsingle',n_cpus=48)
    read_and_save(system='brendaallsingle',n_cpus=48)
    process(system='brendaallsingle',n_cpus=48)
    process_same_test(system='brendaallsingle',test='brenda',n_cpus=48)

    #rhea
    make_singleprod(path="../data/RHEA_reactions.csv",column='rxn_smiles',system='rhea',n_cpus=48)
    read_and_save(system='rhea',n_cpus=48)
    process(system='rhea',n_cpus=48)
    process_same_test(system='rhea',test='brenda',n_cpus=48)

    #metamdb
    make_singleprod(path="../data/MetAMDB_reactions.csv",column='rxn_smiles',system='metamdb',n_cpus=48)
    read_and_save(system='metamdb',n_cpus=48)
    process(system='metamdb',n_cpus=48)
    process_same_test(system='metamdb',test='brenda',n_cpus=48)
