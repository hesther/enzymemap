import pandas as pd
import json
from enzymemap import helpers_brenda, helpers_resolve_smiles, helpers_rdkit, helpers_map
from enzymemap.helpers_brenda import extract_reaction
import re
import itertools
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import collections

def is_ec(s):
    if s:
        if re.match(r"\d+\.\d+\.\d+\.\d+",s):
            return True
    return False

def process_entry(buff):
    buff = buff.strip().replace("\n/"," ").replace("\n^COEFFICIENT","#COEFFICIENT").split("\n")
    d={}
    for line in buff:
        k = line.split(" - ")[0]
        v = " - ".join(line.split(" - ")[1:])
        if k in d.keys():
            d[k].append(v)
        else:
            d[k]=[v]
    return d

def parse(file_loc):
    compounds={}
    buff = ''
    with open(file_loc, errors='replace') as fp:
        new_line = fp.readline()
        while (new_line):
            # Get the next line if we have a blank
            if buff == "":
                if new_line[0]!='#':
                    buff += new_line
            else:
                # If we reach a new BRENDA header
                if new_line == '//\n':
                    d = process_entry(buff)
                    compounds[d['UNIQUE-ID'][0]]=d.get('SMILES',[None])[0]
                    buff = ''
                elif len(buff) > 0 and new_line[0]!='#':
                    buff += new_line
            new_line = fp.readline()
    return compounds

def parse_reac(file_loc):
    ctr=0
    reactions={}
    buff = ''
    with open(file_loc, errors='replace') as fp:
        new_line = fp.readline()
        while (new_line):
            # Get the next line if we have a blank
            if buff == "":
                if new_line[0]!='#':
                    buff += new_line
            else:
                # If we reach a new BRENDA header
                if new_line == '//\n':
                    ctr+=1
                    d = process_entry(buff)
                    if not 'LEFT' in d.keys() or not 'RIGHT' in d.keys():
                        pass
                    else:
                        c = get_coeffs(d['LEFT'],d['RIGHT'])
                        reactions[d['UNIQUE-ID'][0]]=(make_list(d['LEFT'],c[0]),make_list(d['RIGHT'],c[1]),d.get('REACTION-DIRECTION',['UNKNOWN'])[0])
                    buff = ''
                elif len(buff) > 0 and new_line[0]!='#':
                    buff += new_line
            new_line = fp.readline()
    return reactions

def make_list(test,c):
    l=[]
    for x,factor in zip(test,c):
        c = x.split('#COEFFICIENT - ')[0]
        for i in range(factor):
            l.append(c)
    return l

def try_coeffs(coeffs,exchange):
    coeffs = ([str(x) for x in coeffs[0]], [str(x) for x in coeffs[1]])
    for k in exchange.keys():
        coeffs = ([x.replace(k,str(exchange[k]))  for x in coeffs[0]],
                 [x.replace(k,str(exchange[k]))  for x in coeffs[1]])
    coeffs = ([str(eval(x)) for x in coeffs[0]],[str(eval(x)) for x in coeffs[1]])
    if "0.5" in coeffs:
        factor = 2
    else:
        factor = 1
    coeffs = ([int(float(x)*factor) for x in coeffs[0]], [int(float(x)*factor) for x in coeffs[1]])
    return coeffs

def get_coeffs(l,r):
    initial_coeffs = (['1' if '#COEFFICIENT - ' not in x else x.split('#COEFFICIENT - ')[1].replace("<i><b>n</i></b>","n") for x in l],
     ['1' if '#COEFFICIENT - ' not in x else x.split('#COEFFICIENT - ')[1].replace("<i><b>n</i></b>","n") for x in r] )

    exchange={'n':1,'m':1,'q':0,'l':1,'k':1,'j':1,"N":1,"x":1,"y":1,"z":1}
    
    coeffs = try_coeffs(initial_coeffs,exchange)
    while False in [x>0 for side in coeffs for x in side] and exchange['n']<5:
        coeffs = initial_coeffs
        exchange['n']+=1
        coeffs = try_coeffs(coeffs,exchange)
    
    if False in [x>0 for side in coeffs for x in side] and exchange['n']<5:
        sys.exit()
    return coeffs

def get_mols(rxn):
    reacs = [Chem.MolFromSmiles(x) for x in rxn.split(">")[0].split(".")]
    prods = [Chem.MolFromSmiles(x) for x in rxn.split(">")[-1].split(".")]
    return (reacs, prods)
def is_same(m1,m2):
    return m1.HasSubstructMatch(m2) and m2.HasSubstructMatch(m1)
def match(mols,mols_ref):
    matched_ref = []
    matched=[]
    matches=[]
    for i,m1 in enumerate(mols):
        for j,m2 in enumerate(mols_ref):
            if not m2 or j in matched_ref:
                continue
            if is_same(m1,m2):
                matches.append((i,j))
                matched_ref.append(j)
                matched.append(i)
                break
    if len(matched) != len(mols):
        ok=False
    else:
        ok=True
    return ok, matched_ref, matches

def same_formula(r,p):
    mr = CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(r)))
    mp = CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(p)))
    return mr.split('+')[0].split('-')[0]==mp.split('+')[0].split('-')[0]

def same_formula_noh(r,p):
    mr = re.sub(r"H\d*","",CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(r))))
    mp = re.sub(r"H\d*","",CalcMolFormula(Chem.AddHs(Chem.MolFromSmiles(p))))
    return mr.split('+')[0].split('-')[0]==mp.split('+')[0].split('-')[0]

def match_and_correct(mapped,unmapped_ref):
    mols = get_mols(mapped)
    mols_ref = ([Chem.MolFromSmiles(compounds.get(x,"UNK")) if compounds.get(x,"UNK") else None for x in unmapped_ref[0]],
     [Chem.MolFromSmiles(compounds.get(x,"UNK")) if compounds.get(x,"UNK") else None  for x in unmapped_ref[1]])

    #Reactants
    missing_reacs=[]
    #Try native direction:
    ok, matched_ref, matches = match(mols[0],mols_ref[0])
    true_reacs = mapped.split('>')[0]

    reverse=False
    if not ok:
        #try reverse direction:
        ok, matched_ref, matches = match(mols[1],mols_ref[0])
        true_reacs = mapped.split('>')[-1]
        reverse=True
        if not ok:
            print("error!")
    for j in range(len(mols_ref[0])):
        if j not in matched_ref:
            x = reactions[name][0][j]
            if x=='PROTON':
                true_reacs+='.[H+]'
            else:
                missing_reacs.append(x)

    #Products
    missing_prods=[]
    if reverse:
        ok, matched_ref, matches = match(mols[0],mols_ref[1])
        true_prods = mapped.split('>')[0]
    else:
        ok, matched_ref, matches = match(mols[1],mols_ref[1])
        true_prods = mapped.split('>')[-1]
    if not ok:
        print("error!")
    for j in range(len(mols_ref[1])):
        if j not in matched_ref:
            x = reactions[name][1][j]
            if x=='PROTON':
                true_prods+='.[H+]'
            else:
                missing_prods.append(x)
    return missing_reacs, missing_prods, reverse, true_reacs, true_prods

def add_mapnum(smi,add):
    if add==0:
        return smi
    m=Chem.MolFromSmiles(smi)
    [a.SetAtomMapNum(a.GetAtomMapNum()+add) for a in m.GetAtoms() if a.GetAtomMapNum()!=0]
    return Chem.MolToSmiles(m)

def make_rxn_with_cof(missing_reacs, missing_prods, true_reacs, true_prods,cofactors):
    add=0
    for r in missing_reacs:
        if r in cofactors.keys():
            if add>0:
                print('r',r)
            true_reacs+='.'+add_mapnum(cofactors[r],add)
            add+=100
    add=0
    for p in missing_prods:
        if p in cofactors.keys():
            if add>0:
                print('p',p)
            true_prods+='.'+add_mapnum(cofactors[p],add)
            add+=100
    return true_reacs, true_prods

def update(corrected_rxns,reacs,prods,rev,ec):
    if "LEFT-TO-RIGHT" in rev:
        corrected_rxns.append({'ORIG_RXN_TEXT':name,
                               'ORIG_RXNS':remap_duplicates(reacs+">>"+prods),
                               'REVERSIBLE':'ir',
                               'EC_NUM':ec})
    elif rev == "UNKNOWN":
        corrected_rxns.append({'ORIG_RXN_TEXT':name,
                               'ORIG_RXNS':remap_duplicates(reacs+">>"+prods),
                               'REVERSIBLE':'?',
                               'EC_NUM':ec})
    elif "RIGHT-TO-LEFT" in rev:
        corrected_rxns.append({'ORIG_RXN_TEXT':name,
                               'ORIG_RXNS':remap_duplicates(prods+">>"+reacs),
                               'REVERSIBLE':'ir',
                               'EC_NUM':ec})
    elif rev == "REVERSIBLE":
        corrected_rxns.append({'ORIG_RXN_TEXT':name,
                               'ORIG_RXNS':remap_duplicates(reacs+">>"+prods),
                               'REVERSIBLE':'r',
                               'EC_NUM':ec})
    else:
        raise Exception("Unknown reversibility type")
    return corrected_rxns

def remap_duplicates(rxn):
    r=Chem.MolFromSmiles(rxn.split(">")[0])
    dups = [item for item, count in collections.Counter([a.GetAtomMapNum() for a in r.GetAtoms() if a.GetAtomMapNum()!=0]).items() if count > 1]
    if len(dups)==0:
        return rxn
    changes_r={}
    for d in dups:
        changes_r[d]=[]
        for a in r.GetAtoms():
            if a.GetAtomMapNum()==d:
                changes_r[d].append([a.GetIdx(),a.GetSymbol()+str(a.GetFormalCharge()),False])
    p=Chem.MolFromSmiles(rxn.split(">")[-1])
    changes_p={}
    for d in dups:
        changes_p[d]=[]
        for a in p.GetAtoms():
            if a.GetAtomMapNum()==d:
                changes_p[d].append([a.GetIdx(),a.GetSymbol()+str(a.GetFormalCharge()),False])
    num=2001
    for d in dups:
        for i in range(len(changes_r[d])):
            for j in range(len(changes_p[d])):
                if not changes_r[d][i][2] and not changes_p[d][j][2] and changes_r[d][i][1]==changes_p[d][j][1]:
                    r.GetAtomWithIdx(changes_r[d][i][0]).SetAtomMapNum(num)
                    p.GetAtomWithIdx(changes_p[d][j][0]).SetAtomMapNum(num)
                    changes_r[d][i][2]=True
                    changes_p[d][j][2]=True
                    num+=1  
        if False in [x[2] for x in changes_r[d]] or False in [x[2] for x in changes_p[d]]:
            #Second round where we simply assign without checking symbols:
            for i in range(len(changes_r[d])):
                for j in range(len(changes_p[d])):
                    if not changes_r[d][i][2] and not changes_p[d][j][2]:
                        r.GetAtomWithIdx(changes_r[d][i][0]).SetAtomMapNum(num)
                        p.GetAtomWithIdx(changes_p[d][j][0]).SetAtomMapNum(num)
                        changes_r[d][i][2]=True
                        changes_p[d][j][2]=True
                        num+=1 
        if False in [x[2] for x in changes_r[d]] or False in [x[2] for x in changes_p[d]]:
            print("Error in remapping",rxn)
            return rxn
    dups = [item for item, count in collections.Counter([a.GetAtomMapNum() for a in r.GetAtoms() if a.GetAtomMapNum()!=0]).items() if count > 1]
    if len(dups)!=0:
        print("Error in remapping",rxn)
        return rxn
    return Chem.MolToSmiles(r)+">>"+Chem.MolToSmiles(p)







links = pd.read_csv('../data/reaction-links.dat', sep='\t', header=None, skiprows=3, usecols=[0,1])
links = links.dropna()
name_to_ec={}
for i in links.index:
    name_to_ec[links[0][i]] = links[1][i].split("EC-")[-1]

df = pd.read_csv('../data/atom-mappings-smiles.dat', sep='\t', header=None)
df.columns = ['ORIGIN','MAPPED_ORIGINAL']
print("Detected",len(df), "reactions")
df['EC_NUM']=[name_to_ec.get(df['ORIGIN'][i],None) for i in df.index]
df['EC_NUM']=[s if is_ec(s) else None for s in df['EC_NUM'] ]
df=df.dropna()
print("Kept",len(df), "reactions after dropping those with partial/invalid ECs")

all_rxns = []
for i in df.index:
    print(i,end='\r')
    rxns=[]
    if not "<i>" in df['MAPPED_ORIGINAL'][i]:
        for rxn in df['MAPPED_ORIGINAL'][i].split(" "):
            try:
                s,_,p = rxn.split(">")
                rxns.append(helpers_rdkit.unmap(s)+">>"+helpers_rdkit.unmap(p))
            except:
                pass
    if len(rxns)==0:
        rxns=None
    all_rxns.append(rxns)

df['BALANCED_RXNS'] = all_rxns
df=df.dropna()
print("Found",len(df),"reactions that are valid with RDKit")

#drop transferases
df = df[~df['EC_NUM'].str.startswith('7')] # Drop transferases
print("Found",len(df),"reactions after dropping transferases")

#metacyc doesn't contain cofactors and there are a number of erroneous reactions ... correct in the following using the reactions.dat and compounds.dat.file:
compounds = parse("../data/compounds.dat")
reactions = parse_reac("../data/reactions.dat")
df_m = df

#manual corrections. Some of these were just not balanced correctly,
#but some were also unbelievably wrong (like merging two completely different parts,
#eg RXN-12625 or TANNASE-RXN)
idx = df_m[df_m['ORIGIN']=='1.10.3.7-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:1][C:8]2(\\[CH:4]=[C:11]([C:15]3(\\[C:16](=[O:21])[C@:14]1([C:13](/[O:23][CH3:2])=[CH:7]\\[C:9](=[O:18])/[CH:6]=[C:10]1/[C:17]([O:24][CH3:3])=[O:22])[O:20][C:12](/[CH:5]=2)=3))\\[OH:19]).[CH3:101][C:108]2(\\[CH:104]=[C:111]([C:115]3(\\[C:116](=[O:121])[C@:114]1([C:113](/[O:123][CH3:102])=[CH:107]\\[C:109](=[O:118])/[CH:106]=[C:110]1/[C:117]([O:124][CH3:103])=[O:122])[O:120][C:112](/[CH:105]=2)=3))\\[OH:119]).[OH2:25].[OH2:26]>>[CH3:1][C:8]2(\\[CH:4]=[C:11]([C:15](\\[C:16](=[O:21])[C:14]1(/[C:10](/[C:17]([O:24][CH3:3])=[O:22])=[CH:6]\\[C:9](\\[OH:18])=[CH:7]/[C:13](/[O:23][CH3:2])=1))=[C:12](\\[CH:5]=2)\\[O-:20])\\[OH:19]).[CH3:101][C:108]2(\\[CH:104]=[C:111]([C:115](\\[C:116](=[O:121])[C:114]1(/[C:110](/[C:117]([O:124][CH3:103])=[O:122])=[CH:106]\\[C:109](\\[OH:118])=[CH:107]/[C:113](/[O:123][CH3:102])=1))=[C:112](\\[CH:105]=2)\\[O-:120])\\[OH:119]).[O:25]=[O:26]'

idx = df_m[df_m['ORIGIN']=='1.10.3.8-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:1][C:8]3(\\[CH:5]=[C:12]1([C:15](\\[C:16](=[O:21])[C@:14]2([O:20]1)([C:13](/[O:23][CH3:2])=[CH:7]\\[C:9](=[O:18])/[CH:6]=[C:10]2/[C:17]([O:24][CH3:3])=[O:22]))=[C:11](\\[CH:4]=3)\\[OH:19])).[CH3:101][C:108]3(\\[CH:105]=[C:112]1([C:115](\\[C:116](=[O:121])[C@:114]2([O:120]1)([C:113](/[O:123][CH3:102])=[CH:107]\\[C:109](=[O:118])/[CH:106]=[C:110]2/[C:117]([O:124][CH3:103])=[O:122]))=[C:111](\\[CH:104]=3)\\[OH:119])).[OH2:25].[OH2:26]>>[CH3:1][C:8]2(\\[CH:4]=[C:11]([C:15](\\[C:16](=[O:21])[C:14]1(/[C:10](/[C:17]([O:24][CH3:3])=[O:22])=[CH:6]\\[C:9](\\[OH:18])=[CH:7]/[C:13](/[O:23][CH3:2])=1))=[C:12](\\[CH:5]=2)\\[O-:20])\\[OH:19]).[CH3:101][C:108]2(\\[CH:104]=[C:111]([C:115](\\[C:116](=[O:121])[C:114]1(/[C:110](/[C:117]([O:124][CH3:103])=[O:122])=[CH:106]\\[C:109](\\[OH:118])=[CH:107]/[C:113](/[O:123][CH3:102])=1))=[C:112](\\[CH:105]=2)\\[O-:120])\\[OH:119]).[O:25]=[O:26]'

idx = df_m[df_m['ORIGIN']=='1.21.3.1-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:3][C:9]1([CH3:4])([S:22][C@@H:8]2([C@@H:11]([C:14](=[O:21])[N:19]([C@@H:13]([C:16](=[O:2])[O-:2])1)2)[NH:18][C:12](=[O:20])[CH2:7][CH2:5][CH2:6][C@H:10]([NH3+:17])[C:15](=[O:1])[O-:1])).[OH2:23].[OH2:24]>>[CH3:3][CH:9]([CH3:4])[C@@H:13]([NH:19][C:14](=[O:21])[C@H:11]([CH2:8][SH:22])[NH:18][C:12](=[O:20])[CH2:7][CH2:5][CH2:6][C@H:10]([NH3+:17])[C:15](=[O:1])[O-:1])[C:16](=[O:2])[O-:2].[O:23]=[O:24]'

idx = df_m[df_m['ORIGIN']=='BENZOIN-ALDOLASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH:9]2(\\[CH:11]=[CH:13]/[C:15](/[CH:14]([C:6]([C:7]1(\\[CH:5]=[CH:3]/[CH:1]=[CH:2]\\[CH:4]=1))=[O:8])[OH:16])=[CH:12]\\[CH:10]=2)>>[CH:14](=[O:16])[C:15]1(\\[CH:13]=[CH:11]/[CH:9]=[CH:10]\\[CH:12]=1).[CH:6](=[O:8])[C:7]1(\\[CH:5]=[CH:3]/[CH:1]=[CH:2]\\[CH:4]=1)'

idx = df_m[df_m['ORIGIN']=='BILIRUBIN-OXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:3]=[CH:9][C:22]1(\\[C:29](\\[NH:39][C:34]([C:21](/[CH3:8])=1)=[O:40])=[CH:16]\\[C:27]4(\\[NH:37][C:31](/[CH:17]=[C:30]3([C:24](/[CH2:11][CH2:13][C:32](=[O:1])[O-:1])=[C:19]([CH3:6])/[C:26](\\[CH:15]=[C:28]2([C:18](/[CH3:5])=[C:23]([CH:10]=[CH2:4])/[C:35](=[O:41])[NH:38]2))=[N:36]\\3))=[C:25]([CH2:12][CH2:14][C:33](=[O:2])[O-:2])/[C:20](\\[CH3:7])=4)).[CH2:103]=[CH:109][C:122]1(\\[C:129](\\[NH:139][C:134]([C:121](/[CH3:108])=1)=[O:140])=[CH:116]\\[C:127]4(\\[NH:137][C:131](/[CH:117]=[C:130]3([C:124](/[CH2:111][CH2:113][C:132](=[O:101])[O-:101])=[C:119]([CH3:106])/[C:126](\\[CH:115]=[C:128]2([C:118](/[CH3:105])=[C:123]([CH:110]=[CH2:104])/[C:135](=[O:141])[NH:138]2))=[N:136]\\3))=[C:125]([CH2:112][CH2:114][C:133](=[O:102])[O-:102])/[C:120](\\[CH3:107])=4)).[OH2:42].[OH2:43]>>[CH2:3]=[CH:9][C:22]1(\\[C:29](\\[NH:39][C:34]([C:21](/[CH3:8])=1)=[O:40])=[CH:16]\\[C:27]4(\\[NH:37][C:31](/[CH2:17][C:30]3(\\[NH:36][C:26](/[CH:15]=[C:28]2([C:18](/[CH3:5])=[C:23]([CH:10]=[CH2:4])/[C:35](=[O:41])[NH:38]2))=[C:19]([CH3:6])/[C:24](\\[CH2:11][CH2:13][C:32]([O-:1])=[O:1])=3))=[C:25]([CH2:12][CH2:14][C:33]([O-:2])=[O:2])/[C:20](\\[CH3:7])=4)).[CH2:103]=[CH:109][C:122]1(\\[C:129](\\[NH:139][C:134]([C:121](/[CH3:108])=1)=[O:140])=[CH:116]\\[C:127]4(\\[NH:137][C:131](/[CH2:117][C:130]3(\\[NH:136][C:126](/[CH:115]=[C:128]2([C:118](/[CH3:105])=[C:123]([CH:110]=[CH2:104])/[C:135](=[O:141])[NH:138]2))=[C:119]([CH3:106])/[C:124](\\[CH2:111][CH2:113][C:132]([O-:101])=[O:101])=3))=[C:125]([CH2:112][CH2:114][C:133]([O-:102])=[O:102])/[C:120](\\[CH3:107])=4)).[O:42]=[O:43]'

idx = df_m[df_m['ORIGIN']=='CATAL-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[O:2]=[O:1].[OH2:3].[OH2:4]>>[OH:3][OH:2].[OH:4][OH:1]'

idx = df_m[df_m['ORIGIN']=='CATECHOL-OXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH:1]1(\\[CH:2]=[CH:4]/[C:6]([C:5](\\[CH:3]=1)=[O:7])=[O:8]).[CH:101]1(\\[CH:102]=[CH:104]/[C:106]([C:105](\\[CH:103]=1)=[O:107])=[O:108]).[OH2:9].[OH2:10]>>[CH:1]1(\\[CH:2]=[CH:4]/[C:6](/[OH:8])=[C:5](\\[CH:3]=1)/[OH:7]).[CH:101]1(\\[CH:102]=[CH:104]/[C:106](/[OH:108])=[C:105](\\[CH:103]=1)/[OH:107]).[O:9]=[O:10]'

idx = df_m[df_m['ORIGIN']=='DIPEPTIDYL-DIPEPTIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:3][C@H:7]([NH3+:15])[C:11]([NH:14][CH2:5][C:9]([O-:20])=[O:16])=[O:18].[CH3:2][C@H:6]([NH3+:12])[C:10]([NH:13][CH2:4][C:8]([O-:19])=[O:1])=[O:17]>>[CH3:2][C@H:6]([NH3+:12])[C:10]([NH:13][CH2:4][C:8]([NH:15][C@H:7]([C:11]([NH:14][CH2:5][C:9]([O-:20])=[O:16])=[O:18])[CH3:3])=[O:1])=[O:17].[OH2:19]'

idx = df_m[df_m['ORIGIN']=='FATTY-ACID-PEROXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][CH2:3][CH2:4][CH2:5][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH2:12][CH2:13][CH2:14][CH2:15][CH:16]=[O:21].[C:17](=[O:1])=[O:1].[OH2:19].[OH2:20].[OH2:22]>>[CH3:2][CH2:3][CH2:4][CH2:5][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH2:12][CH2:13][CH2:14][CH2:15][CH2:16][C:17]([O-:1])=[O:1].[OH:20][OH:21].[OH:19][OH:22]'

idx = df_m[df_m['ORIGIN']=='HYDROXYBUTYRATE-DIMER-HYDROLASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:3][C@H:7]([CH2:5][C:8]([O-:1])=[O:1])[OH:12].[CH3:2][C@H:6]([CH2:4][C:9]([O-:13])=[O:11])[OH:10]>>[CH3:2][C@H:6]([CH2:4][C:9]([O:12][C@H:7]([CH3:3])[CH2:5][C:8]([O-:1])=[O:1])=[O:11])[OH:10].[OH2:13]'

idx = df_m[df_m['ORIGIN']=='IODIDE-PEROXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[OH2:1].[OH2:2].[I:3][I:4]>>[OH:1][OH:2].[I-:4].[I-:3]'

idx = df_m[df_m['ORIGIN']=='LIGNOSTILBENE-ALPHA-BETA-DIOXYGENASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][O:20][C:16]1(\\[C:14](\\[OH:18])=[CH:8]/[CH:6]=[C:12](/[CH:10]=1)/[CH:4]=[O:22]).[CH3:1][O:19][C:15]1(\\[C:13](\\[OH:17])=[CH:7]/[CH:5]=[C:11](/[CH:9]=1)/[CH:3]=[O:21])>>[CH3:1][O:19][C:15]2(\\[C:13](\\[OH:17])=[CH:7]/[CH:5]=[C:11](/[CH:3]=[CH:4]/[C:12]1(/[CH:10]=[C:16]([C:14](/[OH:18])=[CH:8]\\[CH:6]=1)/[O:20][CH3:2]))/[CH:9]=2).[O:21]=[O:22]'

idx = df_m[df_m['ORIGIN']=='LIMONIN-D-RING-LACTONASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:3][C:23]7([O:34][C@H:18]1([CH2:11][C:19]([O:28][CH2:13][C@:26]61([C@@H:16]([CH2:10][C:17](=[O:29])[C@:25]4([CH3:5])([C@H:15]([CH2:6][CH2:8][C@@:24]2([CH3:4])([C@@H:20]([O:30][C:22](=[O:32])[C@H:21]3([O:35][C@@:27]234))[C:14]5(\\[CH:7]=[CH:9]/[O:33]/[CH:12]=5)))6))7))=[O:1]))([CH3:2]).[OH2:1].[OH2:31]>>[CH3:2][C:23]5([O:34][C@@H:18]([CH2:11][C:19](=[O:1])[O-:1])[C@@:26]4([CH2:13][OH:28])([C@@H:16]([CH2:10][C:17](=[O:29])[C@:25]3([CH3:5])([C@H:15]([CH2:6][CH2:8][C@@:24]([CH3:4])([C@@H:20]([OH:30])[C:14]1(\\[CH:7]=[CH:9]/[O:33]/[CH:12]=1))[C@@:27]2([O:35][C@H:21]([C:22](=[O:32])[O-:31])2)3)4))5))([CH3:3])'

idx = df_m[df_m['ORIGIN']=='MANGANESE-PEROXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[OH2:1].[OH2:2].[Mn+3:3].[Mn+3:4]>>[OH:1][OH:2].[Mn+2:3].[Mn+2:4]'

idx = df_m[df_m['ORIGIN']=='OAHTHAUERA-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:5][C:32]([CH3:6])([C@@H:27]([OH:44])[C:30](=[O:45])[NH:35][CH2:11][CH2:10][C:21](=[O:41])[NH:34][CH2:12][CH2:13][S:54][C:23](=[O:42])[C:14]1(\\[C:22]([CH2:9][CH2:7][CH2:8]/[CH:19]=1)=[O:1]))[CH2:16][O:47][P:53](=[O:4])([O:50][P:52](=[O:3])([O:46][CH2:15][C@@H:20]2([C@@H:26]([O:49][P:51]([O-:2])(=[O:2])[O-:2])[C@@H:25]([OH:43])[C@@H:31]([O:48]2)[N:39]3([C:29]4(\\[N:37]=[CH:17]/[N:36]=[C:28]([C:24](\\[N:38]=[CH:18]/3)=4)/[NH2:33]))))[O-:3])[O-:4].[OH2:1].[OH2:40]>>[CH3:5][C:32]([CH3:6])([C@@H:27]([OH:44])[C:30](=[O:45])[NH:35][CH2:11][CH2:10][C:21](=[O:41])[NH:34][CH2:12][CH2:13][S:54][C:23]([CH2:14][CH:19]([CH2:8][CH2:7][CH2:9][C:22]([O-:1])=[O:1])[OH:40])=[O:42])[CH2:16][O:47][P:53](=[O:4])([O:50][P:52](=[O:3])([O:46][CH2:15][C@@H:20]1([C@@H:26]([O:49][P:51]([O-:2])(=[O:2])[O-:2])[C@@H:25]([OH:43])[C@@H:31]([O:48]1)[N:39]2([C:29]3(\\[N:37]=[CH:17]/[N:36]=[C:28]([C:24](\\[N:38]=[CH:18]/2)=3)/[NH2:33]))))[O-:3])[O-:4]'

idx = df_m[df_m['ORIGIN']=='ORSELLINATE-DEPSIDE-HYDROLASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][C:8]1(\\[C:15](\\[C:17](=[O:21])[O-:23])=[C:12](/[CH:6]=[C:10](\\[CH:4]=1)/[OH:18])\\[OH:19]).[CH3:3][C:9]1(\\[C:14](\\[C:16](=[O:1])[O-:1])=[C:13](/[CH:7]=[C:11](\\[CH:5]=1)/[OH:22])\\[OH:20])>>[CH3:2][C:8]2(/[CH:4]=[C:10](/[CH:6]=[C:12]([C:15](/[C:17](=[O:21])[O:22][C:11]1(\\[CH:7]=[C:13]([C:14](/[C:16](=[O:1])[O-:1])=[C:9]([CH3:3])\\[CH:5]=1)\\[OH:20]))=2)/[OH:19])/[OH:18]).[OH2:23]'

idx = df_m[df_m['ORIGIN']=='PROPIOIN-SYNTHASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][CH2:4][CH:6]=[O:8].[CH3:1][CH2:3][CH:5]=[O:7]>>[CH3:1][CH2:3][CH:5]([C:6](=[O:8])[CH2:4][CH3:2])[OH:7]'

idx = df_m[df_m['ORIGIN']=='RXN-12265'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:4][CH2:8][CH2:9][CH2:12][C@@H:23]([CH3:6])[C@@H:35]([OH:48])[C@@H:30]([OH:47])[CH2:16][C@@H:22]([CH3:5])[CH2:15][C@H:27]([OH:39])[CH2:13][CH2:10][CH2:11][CH2:14][C@@H:28]([OH:40])[CH2:21][C@H:29]([OH:41])[C@@H:24]([NH3+:38])[CH3:7].[CH2:20]([C:34](=[O:44])[O-:2])[CH:26]([C:37](=[O:45])[O-:46])[CH2:18][C:32](=[O:43])[O-:42].[CH2:19]([C:33](=[O:102])[O-:101])[CH:25]([C:36](=[O:3])[O-:3])[CH2:17][C:31](=[O:1])[O-:1]>>[CH3:4][CH2:8][CH2:9][CH2:12][C@@H:23]([CH3:6])[C@@H:35]([O:48][C:34]([CH2:20][C@H:26]([C:37]([O-:46])=[O:45])[CH2:18][C:32](=[O:43])[O-:42])=[O:44])[C@@H:30]([O:47][C:33](=[O:102])[CH2:19][C@H:25]([C:36]([O-:3])=[O:3])[CH2:17][C:31]([O-:1])=[O:1])[CH2:16][C@@H:22]([CH3:5])[CH2:15][C@H:27]([OH:39])[CH2:13][CH2:10][CH2:11][CH2:14][C@@H:28]([OH:40])[CH2:21][C@H:29]([OH:41])[C@@H:24]([NH3+:38])[CH3:7].[OH2:2].[OH2:101]'

idx = df_m[df_m['ORIGIN']=='RXN-12316'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[N:2](=[O:1])(=[O:8])[O-:9].[NH:6]=[O:7].[NH:3]=[O:5].[OH2:4]>>[N:2]([O-:8])=[O:1].[N:3]([O-:4])=[O:5].[N:6]([O-:9])=[O:7]'

idx = df_m[df_m['ORIGIN']=='RXN-12625'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][C:6](=[O:22])[NH:18][C@H:10]1([CH:16]([OH:30])[O:28][C@H:7]([CH2:3][OH:19])[C@@H:11]([OH:23])[C@H:12]([OH:24])1).[CH3:1][C:5](=[O:21])[NH:17][C@H:9]1([CH:15]([OH:26])[O:27][C@H:8]([CH2:4][OH:20])[C@@H:14]([OH:29])[C@H:13]([OH:25])1)>>[CH3:1][C:5](=[O:21])[NH:17][C@H:9]1([CH:15]([OH:26])[O:27][C@H:8]([CH2:4][OH:20])[C@H:14]([C@H:13]([OH:25])1)[O:29][C@@H:16]2([O:28][C@@H:7]([C@@H:11]([OH:23])[C@H:12]([OH:24])[C@@H:10]([NH:18][C:6]([CH3:2])=[O:22])2)[CH2:3][OH:19])).[OH2:30]'

idx = df_m[df_m['ORIGIN']=='RXN-3962'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:11]([CH2:7][CH2:3][CH2:5][CH2:9][C:13](=[O:16])[O-:17])[NH3+:15].[CH2:10]([CH2:6][CH2:2][CH2:4][CH2:8][C:12](=[O:1])[O-:18])[NH3+:14]>>[CH2:11]([CH2:7][CH2:3][CH2:5][CH2:9][C:13](=[O:16])[O-:17])[NH:15][C:12](=[O:1])[CH2:8][CH2:4][CH2:2][CH2:6][CH2:10][NH3+:14].[OH2:18]'

idx = df_m[df_m['ORIGIN']=='RXN-6701'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:23]([S:39][S:19][CH2:5][C@@H:7]([C:10]([NH:14][CH2:4][C:8]([NH2:13])=[O:16])=[O:18])[NH:15][C:9](=[O:17])[CH2:3][CH2:2][C@H:6]([NH3+:12])[C:11]([O-:1])=[O:1])[C@@H:25]([C:28]([NH:32][CH2:22][C:26]([NH2:31])=[O:34])=[O:36])[NH:33][C:27](=[O:35])[CH2:21][CH2:20][C@H:24]([NH3+:30])[C:29]([O-:37])=[O:38].[OH2:40].[OH2:41]>>[CH2:23]([SH:39])[C@@H:25]([C:28]([NH:32][CH2:22][C:26]([NH2:31])=[O:34])=[O:36])[NH:33][C:27](=[O:35])[CH2:21][CH2:20][C@H:24]([NH3+:30])[C:29]([O-:38])=[O:37].[CH2:5]([SH:19])[C@@H:7]([C:10]([NH:14][CH2:4][C:8]([NH2:13])=[O:16])=[O:18])[NH:15][C:9](=[O:17])[CH2:3][CH2:2][C@H:6]([NH3+:12])[C:11]([O-:1])=[O:1].[OH:40][OH:41]'

idx = df_m[df_m['ORIGIN']=='RXN0-1461'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:12]=[CH:8][C:24]2(\\[C:20](/[CH3:4])=[C:28]1([CH2:16][C:29]5(\\[NH:41][C:34](/[CH2:19][C:35]4(\\[NH:43][C:31](/[CH2:18][C:33]3(\\[NH:42][C:30](/[CH2:17][C:32](\\[NH:40]1)=2)=[C:22]([CH3:6])/[C:25](\\[CH:9]=[CH2:13])=3))=[C:23]([CH3:7])/[C:27](\\[CH2:11][CH2:15][C:39]([O-:3])=[O:3])=4))=[C:26]([CH2:10][CH2:14][C:38](=[O:2])[O-:2])/[C:21](\\[CH3:5])=5))).[C:36](=[O:1])=[O:1].[C:37](=[O:45])=[O:44].[OH2:46].[OH2:47]>>[CH3:4][C:20]2(\\[C:24](/[CH2:8][CH2:12][C:36]([O-:1])=[O:1])=[C:32]1([CH2:17][C:30]5(\\[NH:42][C:33](/[CH2:18][C:31]4(\\[NH:43][C:35](/[CH2:19][C:34]3(\\[NH:41][C:29](/[CH2:16][C:28](\\[NH:40]1)=2)=[C:21]([CH3:5])/[C:26](\\[CH2:10][CH2:14][C:38]([O-:2])=[O:2])=3))=[C:27]([CH2:11][CH2:15][C:39]([O-:3])=[O:3])/[C:23](\\[CH3:7])=4))=[C:25]([CH2:9][CH2:13][C:37](=[O:44])[O-:45])/[C:22](\\[CH3:6])=5))).[O:46]=[O:47]'

idx = df_m[df_m['ORIGIN']=='RXN0-1483'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[Fe+3:1].[Fe+3:2].[Fe+3:3].[Fe+3:6].[OH2:4].[OH2:5]>>[Fe+2:1].[Fe+2:2].[Fe+2:3].[Fe+2:6].[O:4]=[O:5]'

idx = df_m[df_m['ORIGIN']=='SUPEROX-DISMUT-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[OH:1][OH:2].[O:4]=[O:3]>>[O-:3][OH:1].[O-:2][OH:4]'

idx = df_m[df_m['ORIGIN']=='TANNASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[C:15]([O-:23])(=[O:21])[C:7]1(\\[CH:4]=[C:9]([C:12](\\[OH:19])=[C:8](/[CH:3]=1)\\[OH:16])\\[OH:17]).[C:14]([O-:1])(=[O:1])[C:6]1(\\[CH:5]=[C:11]([C:13](\\[OH:20])=[C:10](/[CH:2]=1)\\[OH:18])\\[OH:22])>>[C:14](=[O:1])([O-:1])[C:6]1(\\[CH:2]=[C:10]([C:13](\\[OH:20])=[C:11](/[CH:5]=1)\\[O:22][C:15]([C:7]2(\\[CH:4]=[C:9]([C:12](\\[OH:19])=[C:8](/[CH:3]=2)\\[OH:16])\\[OH:17]))=[O:21])\\[OH:18]).[OH2:23]'

idx = df_m[df_m['ORIGIN']=='UROGENDECARBOX-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:19][C:31]2(\\[C:27](/[CH2:11][CH2:15][C:43]([O-:4])=[O:4])=[C:35]1([CH2:22][C:37]5(\\[NH:50][C:34](/[CH2:21][C:36]4(\\[NH:48][C:32](/[CH2:20][C:33]3(\\[NH:49][C:38](/[CH2:23][C:39](\\[NH:51]1)=2)=[C:30]([CH3:18])/[C:25](\\[CH2:9][CH2:13][C:41]([O-:2])=[O:2])=3))=[C:24]([CH2:8][CH2:12][C:40]([O-:1])=[O:1])/[C:28](\\[CH3:16])=4))=[C:26]([CH2:10][CH2:14][C:42](=[O:3])[O-:3])/[C:29](\\[CH3:17])=5))).[C:44](=[O:5])=[O:5].[C:47](=[O:53])=[O:52].[C:45](=[O:6])=[O:6].[C:46](=[O:7])=[O:7]>>[C:40](=[O:1])([O-:1])[CH2:12][CH2:8][C:24]2(\\[C:28](/[CH2:16][C:44](=[O:5])[O-:5])=[C:36]1([CH2:21][C:34]5(\\[NH:50][C:37](/[CH2:22][C:35]4(\\[NH:51][C:39](/[CH2:23][C:38]3(\\[NH:49][C:33](/[CH2:20][C:32](\\[NH:48]1)=2)=[C:25]([CH2:9][CH2:13][C:41](=[O:2])[O-:2])/[C:30](\\[CH2:18][C:46](=[O:7])[O-:7])=3))=[C:31]([CH2:19][C:47](=[O:53])[O-:52])/[C:27](\\[CH2:11][CH2:15][C:43](=[O:4])[O-:4])=4))=[C:29]([CH2:17][C:45]([O-:6])=[O:6])/[C:26](\\[CH2:10][CH2:14][C:42](=[O:3])[O-:3])=5)))'

idx = df_m[df_m['ORIGIN']=='H2-METHYLENE-THMPT-DEHYDRO-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:4][C@H:16]4([C@@H:25]3([C@@H:17]([CH3:5])[N+:39](/[C:19]2(\\[CH:9]=[CH:7]/[C:18](\\[CH2:12][C@H:20]([OH:41])[C@H:27]([OH:43])[C@H:21]([OH:42])[CH2:13][O:47][C@@H:33]1([C@H:29]([OH:45])[C@@H:28]([C@@H:23]([CH2:14][O:48][P:51]([O-:3])(=[O:3])[O:50][C@H:22]([C:32](=[O:2])[O-:2])[CH2:10][CH2:11][C:24](=[O:1])[O-:1])[O:49]1)[OH:44]))=[CH:6]\\[CH:8]=2))=[CH:15]\\[N:40]3[C:26]5(\\[C:31](=[O:46])[NH:38][C:34](/[NH2:35])=[N:37]/[C:30](\\[NH:36]4)=5))).[H][H]>>[CH3:4][C@H:16]4([C@@H:25]3([C@@H:17]([CH3:5])[N:39]([C:19]2(\\[CH:9]=[CH:7]/[C:18](\\[CH2:12][C@H:20]([OH:41])[C@H:27]([OH:43])[C@H:21]([OH:42])[CH2:13][O:47][C@@H:33]1([C@H:29]([OH:45])[C@@H:28]([C@@H:23]([CH2:14][O:48][P:51]([O-:3])(=[O:3])[O:50][C@H:22]([C:32](=[O:2])[O-:2])[CH2:10][CH2:11][C:24](=[O:1])[O-:1])[O:49]1)[OH:44]))=[CH:6]\\[CH:8]=2))[CH2:15][N:40]3[C:26]5(\\[C:31](=[O:46])[NH:38][C:34](/[NH2:35])=[N:37]/[C:30](\\[NH:36]4)=5)))'

idx = df_m[df_m['ORIGIN']=='PROTOPORGENOXI-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:3]=[CH:9][C:23]4(\\[C:19](/[CH3:5])=[C:27]3(\\[CH:15]=[C:28]5([C:21](/[CH3:7])=[C:25]([CH2:11][CH2:13][C:35]([O-:1])=[O:1])/[C:33](\\[CH:18]=[C:34]1([C:26](/[CH2:12][CH2:14][C:36]([O-:2])=[O:2])=[C:22]([CH3:8])/[C:30](/[NH:40]1)=[CH:17]/[C:32]2(/[C:24](\\[CH:10]=[CH2:4])=[C:20]([CH3:6])\\[C:29](\\[N:38]=2)=[CH:16]\\[C:31](\\[NH:37]3)=4)))=[N:39]\\5))).[OH:42][OH:41].[OH:43][OH:44].[OH:45][OH:46]>>[CH2:3]=[CH:9][C:23]2(\\[C:19](/[CH3:5])=[C:27]1([CH2:15][C:28]5(\\[NH:39][C:33](/[CH2:18][C:34]4(\\[NH:40][C:30](/[CH2:17][C:32]3(\\[NH:38][C:29](/[CH2:16][C:31](\\[NH:37]1)=2)=[C:20]([CH3:6])/[C:24](\\[CH:10]=[CH2:4])=3))=[C:22]([CH3:8])/[C:26](\\[CH2:12][CH2:14][C:36]([O-:2])=[O:2])=4))=[C:25]([CH2:11][CH2:13][C:35](=[O:1])[O-:1])/[C:21](\\[CH3:7])=5))).[O:41]=[O:42].[O:43]=[O:44].[O:45]=[O:46]'

idx = df_m[df_m['ORIGIN']=='RXN-8699'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH2:8]1([CH2:7][C:14]4(/[CH:10]=[C:19]([C:20](/[O:25][CH3:3])=[CH:11]\\[C:15](\\[C:17]2(\\[N+:22]\\1=[CH:12]/[C:16]3(\\[C:13](/[CH:9]=2)=[CH:5]/[CH:6]=[C:18]([C:21](\\[O:26][CH3:4])=3)/[O:23][CH3:1])))=4)/[O:24][CH3:2])).[OH:27][OH:28].[OH:29][OH:30]>>[CH2:8]1([CH2:7][C:14]4(/[CH:10]=[C:19]([C:20](/[O:25][CH3:3])=[CH:11]\\[C:15](\\[C@@H:17]2([CH2:9][C:13]3(\\[C:16](\\[CH2:12][N:22]12)=[C:21]([C:18](/[O:23][CH3:1])=[CH:6]\\[CH:5]=3)\\[O:26][CH3:4])))=4)/[O:24][CH3:2])).[O:27]=[O:28].[O:29]=[O:30]'

idx = df_m[df_m['ORIGIN']=='TETRAHYDROBERBERINE-OXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:1][O:22][C:17]1(/[CH:4]=[CH:3]\\[C:12]2(\\[C:15](\\[C:20](\\[O:23][CH3:2])=1)=[CH:10]/[N+:21]3(\\[CH2:6][CH2:5][C:13]4(\\[C:14](\\[C:16](/[CH:7]=2)=3)=[CH:9]/[C:19]5(/[O:25][CH2:11][O:24][C:18](/[CH:8]=4)=5))))).[OH:27][OH:26].[OH:29][OH:28]>>[CH2:6]1([CH2:5][C:13]4(\\[C:14](\\[C@@H:16]2([CH2:7][C:12]3(\\[CH:3]=[CH:4]/[C:17](/[O:22][CH3:1])=[C:20]([C:15](/[CH2:10][N:21]12)=3)/[O:23][CH3:2])))=[CH:9]/[C:19]5(/[O:25][CH2:11][O:24][C:18](/[CH:8]=4)=5))).[O:26]=[O:27].[O:28]=[O:29]'

idx = df_m[df_m['ORIGIN']=='THIAMIN-OXIDASE-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:2][C:9]2(\\[N+:17](/[CH2:7][C:11]1(\\[CH:6]=[N:15]/[C:10](\\[CH3:3])=[N:16]\\[C:13](/[NH2:14])=1))=[CH:8]/[S:18][C:12](/[CH2:4][C:5](=[O:1])[O-:1])=2).[OH:19][OH:20].[OH:21][OH:22]>>[CH3:2][C:9]2(\\[N+:17](/[CH2:7][C:11]1(\\[C:13](\\[NH2:14])=[N:16]/[C:10](\\[CH3:3])=[N:15]\\[CH:6]=1))=[CH:8]/[S:18][C:12](/[CH2:4][CH2:5][OH:1])=2).[OH2:1].[O:19]=[O:20].[O:21]=[O:22]'

idx = df_m[df_m['ORIGIN']=='H2-METHYLENE-THMPT-DEHYDRO-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:4][C@H:16]4([C@@H:25]3([C@@H:17]([CH3:5])[N+:39](/[C:19]2(\\[CH:9]=[CH:7]/[C:18](\\[CH2:12][C@H:20]([OH:41])[C@H:27]([OH:43])[C@H:21]([OH:42])[CH2:13][O:47][C@@H:33]1([C@H:29]([OH:45])[C@@H:28]([C@@H:23]([CH2:14][O:48][P:51]([O-:3])(=[O:3])[O:50][C@H:22]([C:32](=[O:2])[O-:2])[CH2:10][CH2:11][C:24](=[O:1])[O-:1])[O:49]1)[OH:44]))=[CH:6]\\[CH:8]=2))=[CH:15]\\[N:40]3[C:26]5(\\[C:31](=[O:46])[NH:38][C:34](/[NH2:35])=[N:37]/[C:30](\\[NH:36]4)=5))).[H][H]>>[CH3:4][C@H:16]4([C@@H:25]3([C@@H:17]([CH3:5])[N:39]([C:19]2(\\[CH:9]=[CH:7]/[C:18](\\[CH2:12][C@H:20]([OH:41])[C@H:27]([OH:43])[C@H:21]([OH:42])[CH2:13][O:47][C@@H:33]1([C@H:29]([OH:45])[C@@H:28]([C@@H:23]([CH2:14][O:48][P:51]([O-:3])(=[O:3])[O:50][C@H:22]([C:32](=[O:2])[O-:2])[CH2:10][CH2:11][C:24](=[O:1])[O-:1])[O:49]1)[OH:44]))=[CH:6]\\[CH:8]=2))[CH2:15][N:40]3[C:26]5(\\[C:31](=[O:46])[NH:38][C:34](/[NH2:35])=[N:37]/[C:30](\\[NH:36]4)=5)))'

idx = df_m[df_m['ORIGIN']=='H2SOXIREDPYRO-RXN'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[S:1].[H][H]>>[SH2:1]'

idx = df_m[df_m['ORIGIN']=='RXN-18521'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:4][C@@H:13]([O:19][P:20]([O-:2])(=[O:2])[O-:2])[C@H:15]([NH:17][C:14]([CH2:9][CH2:7][CH2:5][CH2:6][CH2:8][CH2:10][SH:21])=[O:18])[C:16](=[O:1])[O-:1].[CH2:12]([CH2:11][SH:22])[S:23](=[O:3])(=[O:3])[O-:3]>>[CH3:4][C@H:13]([C@H:15]([NH:17][C:14]([CH2:9][CH2:7][CH2:5][CH2:6][CH2:8][CH2:10][S:21][S:22][CH2:11][CH2:12][S:23](=[O:3])(=[O:3])[O-:3])=[O:18])[C:16](=[O:1])[O-:1])[O:19][P:20]([O-:2])(=[O:2])[O-:2].[H][H].[H][H]'

idx = df_m[df_m['ORIGIN']=='RXN-7733'].index[0]
df_m['MAPPED_ORIGINAL'][idx]='[CH3:1][C:28]([CH3:2])=[CH:13][CH2:9][CH2:14][C:29](\\[CH3:3])=[CH:15]/[CH2:10][CH2:16][C:30](\\[CH3:4])=[CH:17]/[CH2:11][CH2:18][C:31](\\[CH3:5])=[CH:19]/[CH2:12][CH2:20][CH:32]([CH2:25][CH2:26][O:40][C:33]2(\\[CH:23]=[CH:24]/[C:36]3(/[N:38]=[C:34]1(\\[CH:21]=[CH:7]/[CH:8]=[CH:22]\\[C:35]/1=[N:39]\\[C:37](/[CH:27]=2)=3))))[CH3:6].[H][H]>>[CH3:1][C:28]([CH3:2])=[CH:13][CH2:9][CH2:14][C:29](\\[CH3:3])=[CH:15]/[CH2:10][CH2:16][C:30](\\[CH3:4])=[CH:17]/[CH2:11][CH2:18][C:31](\\[CH3:5])=[CH:19]/[CH2:12][CH2:20][CH:32]([CH3:6])[CH2:25][CH2:26][O:40][C:33]1(\\[CH:23]=[CH:24]/[C:36]2(/[NH:38][C:34]3(\\[C:35](\\[NH:39][C:37](/[CH:27]=1)=2)=[CH:22]/[CH:8]=[CH:7]\\[CH:21]=3)))'


#Missing compounds (manual obtained from atom-mapped file)
compounds['Oxidized-flavodoxins']='Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-])c2cc1C'
compounds['Reduced-flavodoxins']='Cc1cc2c(cc1C)N(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-])c1[nH]c(=O)[nH]c(=O)c1N2'
compounds['Oxidized-Azurins']='[Cu+2]'
compounds['Reduced-Azurins']='[Cu+]'
compounds['HYDROGEN-MOLECULE']='[H][H]'
compounds['Flavodoxins-Semiquinones']='Cc1cc2c(cc1C)N(C[C@H](O)[C@H](O)[C@H](O)COP(=O)([O-])[O-])c1[nH]c(=O)[nH]c(=O)c1N2'
compounds['R-2-Hydroxyisocaproyl-CoA-Dehydratases']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe]2[SH]3[Fe+]41'
compounds['Inactive-R-2OH-Isocaproyl-CoA-Dehydratas']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe+]2[SH]3[Fe+]41'
compounds['Lact-Coa']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe]2[SH]3[Fe+]41'
compounds['Inactive-lactyl-Coa-dehydratases']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe+]2[SH]3[Fe+]41'
compounds['Inactive-R-Phenyllactoyl-CoA-Dehydratase']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe+]2[SH]3[Fe+]41'
compounds['R-Phenyllactoyl-CoA-Dehydratases']='[SH]12[Fe]3[SH]4[Fe]1[SH]1[Fe]2[SH]3[Fe+]41'

cofactors={'NADH': '[NH2:1004][C:1005](=[O:1006])[C:1007]1=[CH:1008][N:1009]([C@@H:1010]2[O:1011][C@H:1012]([CH2:1013][O:1014][P:1015](=[O:1016])([OH:1017])[O:1018][P:1019](=[O:1020])([OH:1021])[O:1022][CH2:1023][C@H:1024]3[O:1025][C@@H:1026]([n:1027]4[cH:1028][n:1029][c:1030]5[c:1031]([NH2:1032])[n:1033][cH:1034][n:1035][c:1036]45)[C@H:1037]([OH:1038])[C@@H:1039]3[OH:1040])[C@@H:1041]([OH:1042])[C@H:1043]2[OH:1044])[CH:1045]=[CH:1046][CH2:1047]1',
           'NAD': '[NH2:1004][C:1005](=[O:1006])[c:1007]1[cH:1008][n+:1009]([C@@H:1010]2[O:1011][C@H:1012]([CH2:1013][O:1014][P:1015](=[O:1016])([OH:1017])[O:1018][P:1019](=[O:1020])([OH:1021])[O:1022][CH2:1023][C@H:1024]3[O:1025][C@@H:1026]([n:1027]4[cH:1028][n:1029][c:1030]5[c:1031]([NH2:1032])[n:1033][cH:1034][n:1035][c:1036]45)[C@H:1037]([OH:1038])[C@@H:1039]3[OH:1040])[C@@H:1041]([OH:1042])[C@H:1043]2[OH:1044])[cH:1045][cH:1046][cH:1047]1',
           'NADPH': '[NH2:1010][C:1011](=[O:1012])[C:1013]1=[CH:1014][N:1015]([C@@H:1016]2[O:1017][C@H:1018]([CH2:1019][O:1020][P:1021](=[O:1022])([OH:1023])[O:1024][P:1025](=[O:1026])([OH:1027])[O:1028][CH2:1029][C@H:1030]3[O:1031][C@@H:1032]([n:1033]4[cH:1034][n:1035][c:1036]5[c:1037]([NH2:1038])[n:1039][cH:1040][n:1041][c:1042]45)[C@H:1043]([O:1044][P:1045](=[O:1046])([OH:1047])[OH:1048])[C@@H:1049]3[OH:1050])[C@@H:1051]([OH:1052])[C@H:1053]2[OH:1054])[CH:1055]=[CH:1056][CH2:1057]1',
           'NADP': '[NH2:1010][C:1011](=[O:1012])[c:1013]1[cH:1014][n+:1015]([C@@H:1016]2[O:1017][C@H:1018]([CH2:1019][O:1020][P:1021](=[O:1022])([OH:1023])[O:1024][P:1025](=[O:1026])([OH:1027])[O:1028][CH2:1029][C@H:1030]3[O:1031][C@@H:1032]([n:1033]4[cH:1034][n:1035][c:1036]5[c:1037]([NH2:1038])[n:1039][cH:1040][n:1041][c:1042]45)[C@H:1043]([O:1044][P:1045](=[O:1046])([OH:1047])[OH:1048])[C@@H:1049]3[OH:1050])[C@@H:1051]([OH:1052])[C@H:1053]2[OH:1054])[cH:1055][cH:1056][cH:1057]1',
           'FAD': '[CH3:1001][c:1002]1[cH:1003][c:1004]2[n:1005][c:1006]3[c:1007](=[O:1008])[nH:1009][c:1010](=[O:1011])[n:1012][c:1013]-3[n:1014]([CH2:1015][C@H:1016]([OH:1017])[C@H:1018]([OH:1019])[C@H:1020]([OH:1021])[CH2:1022][O:1023][P:1024](=[O:1025])([OH:1026])[O:1027][P:1028](=[O:1029])([OH:1030])[O:1031][CH2:1032][C@H:1033]3[O:1034][C@@H:1035]([n:1036]4[cH:1037][n:1038][c:1039]5[c:1040]([NH2:1041])[n:1042][cH:1043][n:1044][c:1045]45)[C@H:1046]([OH:1047])[C@@H:1048]3[OH:1049])[c:1050]2[cH:1051][c:1052]1[CH3:1053]',
           'FADH2': '[CH3:1001][c:1002]1[cH:1003][c:1004]2[c:1050]([cH:1051][c:1052]1[CH3:1053])[N:1014]([CH2:1015][C@H:1016]([OH:1017])[C@H:1018]([OH:1019])[C@H:1020]([OH:1021])[CH2:1022][O:1023][P:1024](=[O:1025])([OH:1026])[O:1027][P:1028](=[O:1029])([OH:1030])[O:1031][CH2:1032][C@H:1033]1[O:1034][C@@H:1035]([n:1036]3[cH:1037][n:1038][c:1039]4[c:1040]([NH2:1041])[n:1042][cH:1043][n:1044][c:1045]34)[C@H:1046]([OH:1047])[C@@H:1048]1[OH:1049])[c:1013]1[c:1006]([c:1007](=[O:1008])[nH:1009][c:1010](=[O:1011])[nH:1012]1)[NH:1005]2'}


corrected_rxns=[]
for name in df_m['ORIGIN']:
    print(name)
    if df_m[df_m['ORIGIN']==name]['BALANCED_RXNS'].values[0]=='[]':
        continue
    mapped=df_m[df_m['ORIGIN']==name]['MAPPED_ORIGINAL'].values[0]
    ec=df_m[df_m['ORIGIN']==name]['EC_NUM'].values[0]
    unmapped_ref=(reactions[name][0],reactions[name][1])
    rev = reactions[name][2]
    missing_reacs, missing_prods, reverse, true_reacs, true_prods = match_and_correct(mapped,unmapped_ref)
#    if sorted(missing_reacs) == sorted(missing_prods):
#        continue
    
    reacs, prods = make_rxn_with_cof(missing_reacs, missing_prods, true_reacs, true_prods,cofactors)

    if same_formula(reacs,prods):
        corrected_rxns = update(corrected_rxns,reacs,prods,rev,ec)
        
    else:
        if 'NAD-P-OR-NOP' in missing_reacs and 'NADH-P-OR-NOP' in missing_prods:
            #make nop
            reacs, prods = make_rxn_with_cof([x if x!='NAD-P-OR-NOP' else 'NAD' for x in missing_reacs],
                                             [x if x!='NADH-P-OR-NOP' else 'NADH' for x in missing_prods],
                                             true_reacs, true_prods,cofactors)
            if same_formula(reacs,prods):
                corrected_rxns = update(corrected_rxns,reacs,prods,rev,ec)

            #make p 
            reacs, prods = make_rxn_with_cof([x if x!='NAD-P-OR-NOP' else 'NADP' for x in missing_reacs],
                                             [x if x!='NADH-P-OR-NOP' else 'NADPH' for x in missing_prods],
                                             true_reacs, true_prods,cofactors)
            if same_formula(reacs,prods):
                corrected_rxns = update(corrected_rxns,reacs,prods,rev,ec)

        elif 'NADH-P-OR-NOP' in missing_reacs and 'NAD-P-OR-NOP' in missing_prods:
            #make nop
            reacs, prods = make_rxn_with_cof([x if x!='NADH-P-OR-NOP' else 'NADH' for x in missing_reacs],
                                             [x if x!='NAD-P-OR-NOP' else 'NAD' for x in missing_prods],
                                             true_reacs, true_prods,cofactors)
            if same_formula(reacs,prods):
                corrected_rxns = update(corrected_rxns,reacs,prods,rev,ec)

            #make p 
            reacs, prods = make_rxn_with_cof([x if x!='NADH-P-OR-NOP' else 'NADPH' for x in missing_reacs],
                                             [x if x!='NAD-P-OR-NOP' else 'NADP' for x in missing_prods],
                                             true_reacs, true_prods,cofactors)
            if same_formula(reacs,prods):
                corrected_rxns = update(corrected_rxns,reacs,prods,rev,ec)

print("Corrected and added cofactors, obtained",len(corrected_rxns),"reactions")
df = pd.DataFrame(corrected_rxns)
print(df)

df['BALANCED_RXNS'] = [[helpers_rdkit.unmap(rxn.split(">")[0])+">>"+helpers_rdkit.unmap(rxn.split(">")[-1])] for rxn in df['ORIG_RXNS']]
df['POSSIBLE_RXNS'] = df['BALANCED_RXNS']
df['SUBSTRATES'] = [[]]*len(df)
df['PRODUCTS'] = [[]]*len(df)
df['PROTEIN_REFS'] = [[]]*len(df)
df['NATURAL']=None
df['ORGANISM']=None
df['PROTEIN_DB']=None


with open("ec_nums.csv",'w') as f:
    for ec_num in sorted(list(set(df['EC_NUM']))):
        print(ec_num,file=f)


df.to_csv("../data/raw_reactions.csv",index=False)


ec_nums =  sorted(list(set(df['EC_NUM'])))
count = []
for ec_num in ec_nums:
    df2 = df[df['EC_NUM']==ec_num]
    c = sum([df2['BALANCED_RXNS'][i]!=[] for i in df2.index])
    count.append(c)
count_df = pd.DataFrame(list(zip(ec_nums,count)),columns=['ec','count'])
print(count_df)
count_df.to_csv("counts.csv",index=False)
