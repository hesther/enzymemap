from .helpers_rdkit import chiral, achiral, unmap, get_balance
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import itertools
import timeout_decorator
import time
from copy import deepcopy
from rdchiral.template_extractor import extract_from_reaction
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants
import pandas as pd

class Reaction():
    """
    helper class for reaction rules
    """
    def __init__(self, reaction_smarts):
        self.rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        mapNoToReactant = {}
        for i,reactant in enumerate(self.rxn.GetReactants()):
            for atom in reactant.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    mapNoToReactant[atom.GetIntProp('molAtomMapNumber')] = (i,atom.GetIdx())     
        self.mapNoToReactant = mapNoToReactant

class Reactants_old():
    def __init__(self,smi,remap=True):
        self.reactants = [Chem.MolFromSmiles(x) for x in smi.split('.')]
        if remap:
            n=0
            for reac in self.reactants:
                [a.SetAtomMapNum(a.GetIdx()+1001+n) for a in reac.GetAtoms()]
                n+=reac.GetNumAtoms() 
        else:
            for reac in self.reactants:
                [a.SetAtomMapNum(a.GetAtomMapNum()+1000) for a in reac.GetAtoms()]
                
def delete_dupl(reac, prod):
    """delete_dupl.

    Delete molecules that occur on the reactant and product site

    Args:
        reac (str): SMILES string of reactants
        prod (str): SMILES string of products

    Returns:
        str: Reaction SMILES
    """
    reacs = reac.split('.')
    prods = prod.split('.')
    
    new_reacs = '.'.join([r for r in reacs if r not in prods])
    new_prods = '.'.join([p for p in prods if p not in reacs])
    
    return new_reacs + ">>" + new_prods

def balance_rxn_protons(rxn_smi,diff_h,diff_h2):
    """balance_rxn_protons.

    Add protons back to reaction string (needed to be previously balanced!).

    Args:
        rxn_smi (str): Reaction SMILES 
        diff_h (int): difference in protons between reactants and products
        diff_h2 (int): difference in molecular hydrogen between reactants and products

    Returns:
        str: Reaction SMILES
    """
    reac, _, prod = rxn_smi.split(">")
    if diff_h > 0:
        reac+=diff_h*'.[H+]'
    elif diff_h < 0:
        prod+=abs(diff_h)*'.[H+]'  
    if diff_h2 > 0:
        reac+=diff_h2*'.[H][H]'
    elif diff_h2 < 0:
        prod+=abs(diff_h2)*'.[H][H]' 
    return reac + ">>" + prod

def reactionFromPermutations(rxn, reactants):
    """reactionFromPermutations.

    Create reactions from all possible permutations of reactants

    Args:
        rxn: RDchiral reaction object
        reactants: RDchiral reactant object

    Returns:
        List(str), List(str): List of mapped and unmapped reaction SMILES
    """    
    rxns = {}
    unique = []

    for perm in itertools.permutations(reactants.reactants):
        tmp = createReactionInstance(rxn,perm)
        for entry in tmp:
            rxn_smi=AllChem.ReactionToSmiles(entry)
            Chem.rdChemReactions.RemoveMappingNumbersFromReactions(entry)
            rxn_smi_unmapped = AllChem.ReactionToSmiles(entry)
            if rxn_smi_unmapped in rxns.keys():
                if rxn_smi not in rxns[rxn_smi_unmapped]:
                    rxns[rxn_smi_unmapped].append(rxn_smi)
            else:
                rxns[rxn_smi_unmapped] = [rxn_smi]
                unique.append(rxn_smi)
    return rxns, unique

def renumber(tres):
    """renumber.

    Renumbers atom map numbers of reaction

    Args:
        tres: RDKit reaction object

    Returns:
        RDKit reaction object
    """      
    roles = [tres.GetReactants(), tres.GetProducts()]
    mapnums = sorted(set([a.GetAtomMapNum() for role in roles for mol in role for a in mol.GetAtoms()]))
    old_to_new = {k: v+1 for v, k in enumerate(mapnums)}
    [a.SetAtomMapNum(old_to_new[a.GetAtomMapNum()]) for role in roles for mol in role for a in mol.GetAtoms()]
    return tres

def createReactionInstance(rxn,reactants):   
    """createReactionInstance.

    For one specific permutation of reactants, create reactions

    Args:
        rxn: RDChiral reaction
        reactants: List of reactants

    Returns:
        List of reactions
    """
    res = []
    try:
        ps = rxn.rxn.RunReactants(reactants)
    except Exception:
        return res
    
    for pset in ps:
        lreacts = [Chem.Mol(x) for x in reactants]
        tres = AllChem.ChemicalReaction()
        for p in pset:
            for atom in p.GetAtoms():
                if atom.HasProp('old_mapno'):
                    mno = atom.GetIntProp('old_mapno')
                    atom.SetIntProp('molAtomMapNumber',mno)
                    ridx,aidx = rxn.mapNoToReactant[mno] 
                    raidx = int(atom.GetProp("react_atom_idx"))
                    lreacts[ridx].GetAtomWithIdx(raidx).SetIntProp('molAtomMapNumber',mno)
            tres.AddProductTemplate(p)
        for reactant in lreacts:
            tres.AddReactantTemplate(reactant)
        tres = renumber(tres)
        res.append(tres)
    return res

def has_correct(rxn_smi):
    """has_correct.

    Check whether the reaction can be reproduced by extracting an RDChiral template and applying it

    Args:
        rxn_smi (str): reaction SMILES

    Returns:
        bool: Whether the reaction was reproducable
    """
    reaction = rxn_smi.replace('.[H+]','')
    try:
        rxn = rdchiralReaction(get_template(reaction))
        srct = reaction.split(">")[0]
        prod = unmap(reaction.split(">")[-1])
        rct = rdchiralReactants(srct)
        outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False)
        return prod in outcomes
    except:
        return False

def get_template(rxn_smi, radius=1, nsg=True):
    """get_template.

    Extract RDChiral template

    Args:
        rxn_smi (str): reaction SMILES
        radius (int): Radius of rule
        nsg (int): Whether to exclude special groups

    Returns:
        RDChiral template
    """
    rxn_split = rxn_smi.split(">")
    reaction={"_id":0,"reactants":rxn_split[2],"spectator":rxn_split[1],"products":rxn_split[0]}
    try:
        template = extract_from_reaction(reaction, radius=radius, no_special_groups=nsg)["reaction_smarts"]
        template = "(" + template.replace(">>", ")>>")
    except:
        template = None
    return template

def propose_new(rxn_smi):
    """propose_new.

    Tries to correct stereochemistry in a reaction

    Args:
        rxn_smi (str): reaction SMILES

    Returns:
        New reaction SMILES (or None if not possible)
    """
    getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),2,useChirality=False, useFeatures=True)
    getfp_chi = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),2,useChirality=True, useFeatures=True)
    similarity_metric = DataStructs.BulkTanimotoSimilarity

    reaction = rxn_smi
    srct = reaction.split(">")[0]
    prod = unmap(reaction.split(">")[-1])
    p_fp = getfp(prod)
    p_fp_chi = getfp_chi(prod)
    rct = rdchiralReactants(srct)
    
    try:
        rxn = rdchiralReaction(get_template(reaction))
        outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False,keep_mapnums=True)
    except:
        try:
            rxn = rdchiralReaction(get_template(reaction, nsg=False))
            outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False,keep_mapnums=True)
        except:
#            print("could not run reaction", rxn_smi)
            outcomes = []

    #Need to include a loop here to filter out outcomes that are the same without mapnums!
    known=[]
    new_outcome=None
    max_sim_chi=0
    for outcome in outcomes:
        unmapped = unmap(outcome)
        if achiral(unmapped) != achiral(prod):
            continue
        if unmapped not in known:
            known.append(unmapped)
            p_ref_fp = [getfp(outcome)]
            p_ref_fp_chi = [getfp_chi(outcome)]
            sim = similarity_metric(p_fp, p_ref_fp)[0]
            if abs(sim-1)<0.001:
                sim_chi = similarity_metric(p_fp_chi, p_ref_fp_chi)[0]
                if sim_chi>max_sim_chi:
                    new_outcome = outcome
    if new_outcome != None:
        return Chem.MolToSmiles(rct.reactants)+">>"+new_outcome
    else:
        return None
    
@timeout_decorator.timeout(30)
def get_mapped_reacs(reacs, prods, rules):
    """get_mapped_reacs.

    Map reaction via single rule application

    Args:
        reacs (str): SMILES of reactants
        prods (str): SMILES of products
        rules (pandas.DataFrame): Pandas Dataframe of reaction rules

    Returns:
        str, str, int: Mapped reaction SMILES, SMARTS rule, rule id
    """
    #Remove protons since they don't get treated correctly in the reaction rules from Broadbelt et al
    diff_h = reacs.split(".").count('[H+]') - prods.split(".").count('[H+]')
    diff_h2 = reacs.split(".").count('[H][H]') - prods.split(".").count('[H][H]')
    reacs= ".".join([x for x in reacs.split(".") if x!= '[H+]' and x != '[H][H]'])
    prods= ".".join([x for x in prods.split(".") if x!= '[H+]' and x != '[H][H]'])
    
    r_a = achiral(reacs)
    p_a = achiral(prods)
    r_c = chiral(reacs)
    p_c = chiral(prods)
    num_r = len(r_a.split('.'))
    num_p = len(p_a.split('.'))
    mr_c = Chem.MolFromSmiles(r_c)
    mp_c = Chem.MolFromSmiles(p_c)
    
    reactants = Reactants_old(r_a)

    results=[]
    applied_rules = []
    rules_id = []
    for i in rules.index:
        if rules['Num_R'][i] != num_r or rules['Num_P'][i] != num_p:
            #Wrong number of products, don't even need to try template.        
            continue
        rxn = rules['Rxn'][i]    
        d, lst = reactionFromPermutations(rxn, reactants)
        for key in d.keys():
            _, p_u = key.split('>>')
            if p_u == p_a:
                
                for mapped_rxn in d[key][:1]: #Use only first
                    r_m, p_m = mapped_rxn.split('>>')
                    mr_a = Chem.MolFromSmiles(r_m)
                    match=mr_a.GetSubstructMatch(mr_c)
                    [a.SetAtomMapNum(mr_a.GetAtomWithIdx(match[a.GetIdx()]).GetAtomMapNum()) for a in mr_c.GetAtoms()]
                    
                    mp_a = Chem.MolFromSmiles(p_m)
                    match=mp_a.GetSubstructMatch(mp_c)
                    [a.SetAtomMapNum(mp_a.GetAtomWithIdx(match[a.GetIdx()]).GetAtomMapNum()) for a in mp_c.GetAtoms()]
                    
                    rxn_smi = Chem.MolToSmiles(mr_c) + '>>' + Chem.MolToSmiles(mp_c)
                    
                    if not has_correct(rxn_smi):
                        #Correct stereochem
                        rxn_smi = propose_new(rxn_smi)
                        if not rxn_smi:
                            continue
                    rxn_smi = balance_rxn_protons(rxn_smi,diff_h,diff_h2)

                    return rxn_smi, rules['SMARTS'][i], i

    return None, None, None

def initial_map(smi):
    """inital_map.

    Make map numbers for an unmapped compound

    Args:
        smi (str): SMILES string

    Returns:
        str, List: SMILES string, list of individual molecules
    """
    r = [Chem.MolFromSmiles(x) for x in smi.split(".")]
    ctr=1
    for x in r:
        for atom in x.GetAtoms():
            atom.SetAtomMapNum(ctr)
            ctr += 1
    return ".".join([Chem.MolToSmiles(x) for x in r]), r

@timeout_decorator.timeout(20)
def createReactionInstance_multi(rxn,reactants_unshuffled): 
    """createReactionInstance.

    For one specific permutation of reactants, create reactions intended for multi-step rule applciations (which need to track atoms differently)

    Args:
        rxn: RDChiral reaction
        reactants_unshuffled: List of reactants

    Returns:
        List of reactions, list of update in mapnumbers
    """
    res = []
    changed =[]
    for old_reactants in itertools.permutations(reactants_unshuffled):
        if len(old_reactants)>rxn.rxn.GetNumReactantTemplates():
            reactants = old_reactants[:rxn.rxn.GetNumReactantTemplates()]
            unused_reactants = old_reactants[rxn.rxn.GetNumReactantTemplates():]
        else:
            reactants = old_reactants
            unused_reactants = []
        
        try:
            ps = rxn.rxn.RunReactants(reactants)
        except Exception:
            return res, changed
    
        for pset in ps:
            lreacts = [Chem.Mol(x) for x in reactants]
            tres = AllChem.ChemicalReaction()
            updated_mnos = []
            for p in pset:
                for atom in p.GetAtoms():
                    if atom.HasProp('old_mapno'):
                        mno = atom.GetIntProp('old_mapno')
                        ridx,aidx = rxn.mapNoToReactant[mno] 
                        raidx = int(atom.GetProp("react_atom_idx"))
                        mno2 = lreacts[ridx].GetAtomWithIdx(raidx).GetAtomMapNum()
                        atom.SetIntProp('molAtomMapNumber',mno2)
                        updated_mnos.append(mno2)
                tres.AddProductTemplate(p)
            for reactant in lreacts:
                tres.AddReactantTemplate(reactant)
            for r in unused_reactants:
                tres.AddReactantTemplate(r)
                tres.AddProductTemplate(r)
            if is_valid_molecules(tres):   
                res.append(tres)
                changed.append(updated_mnos)
    return res, changed

def is_valid_molecules(tres):
    """is_valid_molecules.

    Check whether molecules in reaction string are all valid in RDKit

    Args:
        tres: RDKit reaction

    Returns:
        bool: Boolean whether molecules in reaction are valid
    """
    rxn_smi = AllChem.ReactionToSmiles(tres)
    for smi in rxn_smi.split(">>"):
        if not Chem.MolFromSmiles(smi):
            return False
        if not Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smi))):
            return False
    return True

def make_mol_from_mapped(smi):
    """make_mol_from_mapped.

    This is necessary because rdkit discerns betweet "default implicit hydrogens", e.g from smi='C' and
    "specified implicit hydrogens, e.g. from smi='[CH4:1]'. Since we sometimes need to make molecules from
    mapped smiles strings, we need to make the hydrogens "default implicit" for the reaction functionality
    to work correctly. This is a rather unelegant workaround, but it works:

    Args:
        smi (str): SMILES string

    Returns:
        List(rdkit.Chem.mol): List of RDKit molecules
    """

    mols = []
    for s in smi.split("."):
        mol = Chem.MolFromSmiles(s)
        mapnos = [atom.GetAtomMapNum() for atom in mol.GetAtoms()]
        [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
        mol2 = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        match = mol2.GetSubstructMatch(mol)
        [mol2.GetAtomWithIdx(m).SetAtomMapNum(mapnos[i]) for i, m in enumerate(match)]
        mols.append(mol2)
    return mols


def get_multistep(rxn, reac_smi, true_prod_smi,MAX_STEPS=2):
    """get_multistep.

    Get atom mapping via multiple mapping steps

    Args:
        rxn: RDChiral reaction
        reac_smi (str): SMILES string of reactants
        true_prod_smi (str): SMILES string of products
        MAX_STEPS (int): maximum number of consecutive steps

    Returns:
        List, List, List: List of mapped reactions, corresponding rules and changed atoms
    """
    reac_smi, reac_mol = initial_map(reac_smi)
    mapped = [(reac_smi, reac_mol, [[]])]
    unfinished = 1
    steps =0
    start=time.time()
    end = time.time()
    
    while unfinished > 0 and steps < MAX_STEPS and end-start<60:
        new_mapped = []
        for m in mapped:
            if end-start>60: #Also break within loop if it takes too long
                break
            reac_smi, reac_mol, c_old = m
            
            if reac_mol is not None:
                results, changed = createReactionInstance_multi(rxn, reac_mol)
                try:
                    results, changed = createReactionInstance_multi(rxn, reac_mol)
                except:
                    results = []
                end = time.time()
                if len(results) > 0:
                    rxns_unmapped = []
                    for entry, c in zip(results, changed):
                        entry_copy = deepcopy(entry)
                        corrected = AllChem.ReactionToSmiles(entry) 
                        Chem.rdChemReactions.RemoveMappingNumbersFromReactions(entry_copy)
                        rxn_smi_unmapped = AllChem.ReactionToSmiles(entry_copy)
                        if rxn_smi_unmapped not in rxns_unmapped:
                            new_mapped.append((reac_smi+">>"+corrected.split(">>")[-1],
                                               make_mol_from_mapped(corrected.split(">>")[-1]),
                                               c_old+[c_old[-1]+c]))
                            rxns_unmapped.append(rxn_smi_unmapped)
                else:
                    new_mapped.append((reac_smi, None, c_old))
            else:
                new_mapped.append((reac_smi, None, c_old))

        mapped = new_mapped
        unfinished = sum([m[1] is not None for m in mapped])
        steps += 1
        end = time.time()

    results = [m[0] for m in mapped if ">>" in m[0]]
    changed = [m[2] for m in mapped if ">>" in m[0]]

    individual = {}
    overall = {}
    true_prod_smi = Chem.MolToSmiles(Chem.MolFromSmiles(true_prod_smi))
    for result,c in zip(results,changed):
        unmapped = [unmap(x) for x in result.split(">>")]
        mapped = result.split(">>")

        for i, (m, um) in enumerate(zip(mapped, unmapped)):
            if i == 0:
                continue
            if um == true_prod_smi:
                unmapped_rxn = unmapped[0]+">>"+um 
                mapped_rxn = mapped[0]+">>"+m
                if unmapped_rxn not in overall.keys():
                    #Found reaction that produces product
                    #Add overall reaction
                    overall[unmapped_rxn] = mapped_rxn
                if i>1:    
                    #Add individual reactions
                    for j in range(i):
                        mapped_rxn = delete_dupl(mapped[j], mapped[j+1])
                        #Need to unmap again because ordering might change when molecules are deleted
                        unmapped_rxn = unmap(mapped_rxn.split(">>")[0]) + ">>" + unmap(mapped_rxn.split(">>")[-1])
                    
                        if unmapped_rxn not in individual.keys():
                            individual[unmapped_rxn] = {'rxn':mapped_rxn, 'changed':(c[j],c[j+1])}

    return list(overall.values()), [individual[k]['rxn'] for k in individual.keys()], [individual[k]['changed'] for k in individual.keys()]

def correct_stereochem(changed_atoms_nos, mol1, mol3, mol2):
    """correct_stereochem.

    Correct stereochemistry for two-step reactions for intermediats

    Args:
        changed_atoms_nos: Numbers of changed atoms in the reactions
        mol1: RDKit molecule (before first step)
        mol2: RDKit molecule (after first step)
        mol3: RDKit molecule (after second step)

    Returns:
        Corrected molecule of intermediates
    """
    d_mol1 = {}
    for atom in mol1.GetAtoms():
        d_mol1[atom.GetAtomMapNum()] = atom.GetChiralTag()
    
    d_mol3 = {}
    for atom in mol3.GetAtoms():
        d_mol3[atom.GetAtomMapNum()] = atom.GetChiralTag()

    for atom in mol2.GetAtoms():
        mapno = atom.GetAtomMapNum()
        if mapno in changed_atoms_nos:
            atom.SetChiralTag(d_mol3[mapno])
        else:
            atom.SetChiralTag(d_mol1[mapno])
    return mol2

@timeout_decorator.timeout(60)
def get_mapped_reacs_multi(reacs, prods, rules, MAX_STEPS=2): 
    """get_mapped_reacs_multi.

    Allows for multistep reactions (up to MAX_STEPS) and extracts the overall reactions
    as well as the corresponding single steps. For single steps, stereochemistry is obtained by copying
    information, so might not be 100% accurate for some edge cases (but shouldn't introduce more
    noise than there already is in BRENDA concerning stereochemistry).

    Args:
        reacs (str): SMILES of reactants
        prods (str): SMILES of products
        rules (pandas.DataFrame): Pandas Dataframe of reaction rules
        MAX_STEPS (int): maximum number of steps for mapping

    Returns:
        str, str, int, List(str): Mapped reaction SMILES, SMARTS rule, rule id, List of single-step reactions
    """


    
    #Remove protons since they don't get treated correctly in the reaction rules from Broadbelt et al
    diff_h = reacs.split(".").count('[H+]') - prods.split(".").count('[H+]')
    diff_h2 = reacs.split(".").count('[H][H]') - prods.split(".").count('[H][H]')
    reacs= ".".join([x for x in reacs.split(".") if x!= '[H+]' and x != '[H][H]'])
    prods= ".".join([x for x in prods.split(".") if x!= '[H+]' and x != '[H][H]'])
        
    r_a = achiral(reacs)
    p_a = achiral(prods)
    r_c = chiral(reacs)
    p_c = chiral(prods)
    mr_c = Chem.MolFromSmiles(r_c)
    mp_c = Chem.MolFromSmiles(p_c)

    for i in rules.index:
        
        rxn = rules['Rxn'][i]  
        overall, individual, changed_individual = get_multistep(rxn, r_a, p_a, MAX_STEPS)
        
        if len(overall) == 0:
            continue
        elif len(overall) >1:
            print("Cauton, more than one mapping possibility detected", reacs, prods)
        mapped_rxn = overall[0]

        r_m, p_m = mapped_rxn.split('>>')
        mr_a = Chem.MolFromSmiles(r_m)
        match=mr_a.GetSubstructMatch(mr_c)
        [a.SetAtomMapNum(mr_a.GetAtomWithIdx(match[a.GetIdx()]).GetAtomMapNum()) for a in mr_c.GetAtoms()]
                    
        mp_a = Chem.MolFromSmiles(p_m)
        match=mp_a.GetSubstructMatch(mp_c)
        [a.SetAtomMapNum(mp_a.GetAtomWithIdx(match[a.GetIdx()]).GetAtomMapNum()) for a in mp_c.GetAtoms()]
                    
        rxn_smi = Chem.MolToSmiles(mr_c) + '>>' + Chem.MolToSmiles(mp_c)

        if not has_correct(rxn_smi):
            #Correct stereochem
            rxn_smi = propose_new(rxn_smi)
            if not rxn_smi:
                continue
        rxn_smi = balance_rxn_protons(rxn_smi,diff_h,diff_h2)
                        
        results_single=[]
        for ii, c in zip(individual, changed_individual):  
            mir_c = correct_stereochem(c[0], mr_c, mp_c, Chem.MolFromSmiles(ii.split(">>")[0]))
            mip_c = correct_stereochem(c[1], mr_c, mp_c, Chem.MolFromSmiles(ii.split(">>")[-1]))
            
            rxn_smi_i = Chem.MolToSmiles(mir_c) + '>>' + Chem.MolToSmiles(mip_c)
            #Cutting the number of protons and H2 in half only makes sense for two-step FIXME
            rxn_smi_i = balance_rxn_protons(rxn_smi_i,int(diff_h/2),int(diff_h2/2))
            results_single.append(rxn_smi_i)

        return rxn_smi, rules['SMARTS'][i], i, results_single
            
    return None, None, None, None

def map(rxns, rules, single=True):
    """map.

    Map reactions via single or multi-step rule application

    Args:
        rxns (List(str)): List of reaction SMILES
        rules (pandas.DataFrame): Pandas Dataframe of reaction rules
        single (bool): Whether to only allow single-step reactions

    Returns:
        List(str), List(str), List(int), List(List(str)): Mapped reaction SMILES, SMARTS rule, rule id, List of single-step reactions
    """
    
    res=[]
    rls=[]
    iis=[]
    indis=[]
    cached=dict()
    for rxn in rxns:
        if rxn not in cached.keys():
            reacs, _, prods = rxn.split(">")

            if single:
                #Try isomerase first:
                result = map_isomerase(rxn)
                if result:
                    rule=None
                    idx=-1
                else:
                    try:
                        result, rule, idx = get_mapped_reacs(reacs, prods, rules)
                    except:
                        result = None
                        rule = None
                        idx = None
                indi = None
            else:
                try:
                    result, rule, idx, indi = get_mapped_reacs_multi(reacs, prods, rules)
                except:
                    result = None
                    rule = None
                    idx = None
                    indi = None

            cached[rxn]=(result, rule, idx, indi)
        else:
            result, rule, idx, indi = cached[rxn]
            
        if not result:
            continue

        res.append(result)
        rls.append(rule)
        iis.append(idx)
        indis.append(indi)
    
    return res, rls, iis, indis       


def make_templates_for_suggestions(reactions_in_ec):
    """make_templates_for_suggestions.

    Make very general RDChiral templates for suggesting alternate products for unbalanced or unmapped reactions

    Args:
        reactions_in_ec (List(str)): List of reaction SMILES

    Returns:
        List(rdchiralReaction), dict, dict, List(str): List of reactions, dictionary with hydrogen counts, dictionary with reactants, sorted templates
    """

    templates = [get_template(x, radius=0, nsg=True) for x in reactions_in_ec]
    temp2reac = {}
    temp2h={}
    for x,r in zip(templates,reactions_in_ec):
        if x is None:
            continue
        if x not in temp2reac.keys():
            temp2reac[x]=[r]
            temp2h[x]=r.split(">>")[-1].count('.[H+]')
        else:
            temp2reac[x].append(r)
            if r.split(">>")[-1].count('.[H+]') != temp2h[x]:
                print('error')
    template_s = sorted(list(set([x for x in templates if x is not None])))
    templates = [rdchiralReaction(x) for x in template_s]

    return templates, temp2h, temp2reac, template_s

def compute_sim(p, r, reactions):
    """compute_sim.

    Computes fingerprint similarities to all reactions

    Args:
        p (str): SMILES of product
        r (str): SMILES of reactant
        reactions (List(str)): List of reaction SMILES

    Returns:
        int: maximum similarity to reference reactions
    """
    getfp = lambda smi: AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smi),2,useChirality=True, useFeatures=True)
    similarity_metric = DataStructs.BulkTanimotoSimilarity

    p_fp = getfp(p)
    p_ref_fp = [getfp(x.split(">>")[-1]) for x in reactions]
    p_sims = similarity_metric(p_fp, p_ref_fp)
    
    r_fp = getfp(r)
    r_ref_fp = [getfp(x.split(">>")[0]) for x in reactions]
    r_sims = similarity_metric(r_fp, r_ref_fp)
    
    sims = [p_sims[i]*r_sims[i] for i in range(len(reactions))]
    return max(sims)

@timeout_decorator.timeout(20)
def suggest_corrections(rxns, templates, temp2h, temp2reac, template_s):
    """suggest_corrections.

    Suggest reactions for unmapped or unbalanced entries

    Args:
        rxns: List of reactions
        templates: reaction templates
        temp2h: dictionary with hydrogen counts
        temp2reac: dictionary with reactants
        template_s: sorted templates

    Returns:
        List(str): List of possible reactions
    """
    all_outcomes=[]
    for reaction in rxns:
        final_outcome=None
        final_sim=-1
        srct = reaction.split(">")[0]
        rct = rdchiralReactants(srct)
        for j, (rxn, sma) in enumerate(zip(templates,template_s)):
            if rxn is None:
                continue
            outcomes = rdchiralRun(rxn, rct, combine_enantiomers=False)
            outcomes = ['.'.join([x]+['[H+]']*temp2h[sma]) for x in outcomes]
            outcomes = [x for x in outcomes if get_balance(srct+">>"+x)]
            sims = [compute_sim(x, srct, temp2reac[sma]) for x in outcomes]
            for outcome, sim in zip(outcomes,sims):
                if sim>final_sim:
                    final_outcome=outcome
                    final_sim=sim
        if final_outcome is not None:
            all_outcomes.append(srct+">>"+final_outcome)
        else:
            all_outcomes.append(None)
    return [x for x in all_outcomes if x is not None]

def get_groups(rules):
    """get_groups.

    Get group of reactions rules that catalyze the same reactions (in forward and reverse direction)

    Args:
        rules: Pandas dataframe of rules

    Returns:
        List of rule groups
    """
    groups = []
    done = set()
    for key in rules.index:
        if key not in done:
            update=set([key]+rules['reverse_ids'][key])
            for key2 in rules['reverse_ids'][key]:
                update.update(rules['reverse_ids'][key2])
            groups.append(update)
            done.update(update)
    return groups

def probably_reversible(df, rules):
    """probably_reversible.

    Suggest reversibility if reversible tag (reversible or irreversible) is not specified

    Args:
        df: Pandas dataframe
        rules: Pandas dataframe of rules

    Returns:
        List of suggested reversibility tags
    """
    prob_rev=[]
    for i in df.index:
        
        if df['reversible'][i] == 'r':
            prob_rev.append(['r']*len(df['rule_ids'][i]))
            continue
        elif df['reversible'][i] == 'ir':
            prob_rev.append(['ir']*len(df['rule_ids'][i]))
            continue
        
        #Try to see whether rule should be reversed
        l = list(set([item for sublist in df['rule_ids'].values for item in sublist if item != None]))

        p=[]
        for j in range(len(df['rule_ids'][i])):
            if list(set(l).intersection(rules['reverse_ids'][df['rule_ids'][i][j]])):
                p.append('p_r')
            else:
                p.append('p_ir') 
        prob_rev.append(p)
    return prob_rev

def standardize_mapping_rxn(rxn_smi):
    """standardize_mapping_rxn.

    Standardizes a reaction mapping

    Args:
        rxn_smi (str): reaction SMILES

    Returns:
        str: Updated reaction SMILES
    """
    r,_,p = rxn_smi.split('>')
    mol_r_um = Chem.MolFromSmiles(unmap(r))
    mol_p_m = Chem.MolFromSmiles(p)
    mol_r_m = Chem.MolFromSmiles(r)
    
    mapnum2idx = {}
    for atom in mol_p_m.GetAtoms():
        mapnum2idx[atom.GetAtomMapNum()] = atom.GetIdx()   
    match = mol_r_m.GetSubstructMatch(mol_r_um)
    for idx_um,idx_m in enumerate(match):
        atom_r_m = mol_r_m.GetAtomWithIdx(idx_m)
        atom_r_um = mol_r_um.GetAtomWithIdx(idx_um)
        new_mapnum = atom_r_um.GetIdx() + 1
        old_mapnum = atom_r_m.GetAtomMapNum()
        if old_mapnum != 0:
            atom_p = mol_p_m.GetAtomWithIdx(mapnum2idx[old_mapnum])
            atom_r_um.SetAtomMapNum(new_mapnum)
            atom_p.SetAtomMapNum(new_mapnum)
        else:
            if atom_r_m.GetSymbol() != 'H':
                raise ValueError("unmapped atom")
        
    r = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol_r_um)))
    p = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol_p_m)))
    
    return r + ">>" + p

def reverse(r):
    """reverse.

    Reverse reaction

    Args:
        r (str): reaction SMILES

    Returns:
        str: Reversed reaction SMILES
    """
    a,b,c=r.split('>')
    return '>'.join([c,b,a])

def make_final(df):
    """make_final.

    Construct final dataframe with new rows for each reaction, add reversible reactions

    Args:
        df: Pandas Dataframe

    Returns:
        New dataframe with reactions on separate rows instead of lists within cells
    """
    orig_rxn_text = []
    mapped = []
    unmapped = []
    rule = []
    rule_id = []
    steps = []
    source = []
    quality = []
    rxn_idx = []
    natural = []
    organism = []
    protein_refs = []
    protein_db = []

    for i in df.index:
        print(i,end='\r')
        for j in range(len(df['mapped_rxns'][i])):
            rxn_idx.append(i)
            orig_rxn_text.append(df['rxn_text'][i])
            result = df['mapped_rxns'][i][j]
            mapped.append(result)
            unmapped.append(unmap(result.split('>')[0])+">>"+unmap(result.split('>')[-1]))
            rule.append(df['rules'][i][j])
            rule_id.append(df['rule_ids'][i][j])
            quality.append(df['quality'][i][j])
            source.append(df['source'][i])
            steps.append(df['step'][i])
            natural.append(df['natural'][i])
            organism.append(df['organism'][i])
            protein_refs.append(df['protein_refs'][i])
            protein_db.append(df['protein_db'][i])
        
            if df['prob_rev'][i][j]=='r' or df['prob_rev'][i][j]=='p_r':
                rxn_idx.append(i)
                orig_rxn_text.append(df['rxn_text'][i])
                result = df['mapped_rxns'][i][j]
                mapped.append(reverse(result))
                unmapped.append(reverse(unmap(result.split('>')[0])+">>"+unmap(result.split('>')[-1])))
                rule.append(df['rules'][i][j])
                rule_id.append(df['rule_ids'][i][j])
                quality.append(df['quality'][i][j])
                if df['prob_rev'][i][j]=='r':
                    source.append(df['source'][i]+' reversed')
                else:
                    source.append(df['source'][i]+' reversed_suggested')
                steps.append(df['step'][i])
                natural.append(df['natural'][i])
                organism.append(df['organism'][i])
                protein_refs.append(df['protein_refs'][i])
                protein_db.append(df['protein_db'][i])
            
            if df['step'][i] == 'multi':
                # also add individuals
                for result in df['individuals'][i][j]:
                    rxn_idx.append(i)
                    orig_rxn_text.append(df['rxn_text'][i])
                    mapped.append(result)
                    unmapped.append(unmap(result.split('>')[0])+">>"+unmap(result.split('>')[-1]))
                    rule.append(df['rules'][i][j])
                    rule_id.append(df['rule_ids'][i][j])
                    quality.append(df['quality'][i][j])
                    source.append(df['source'][i])
                    steps.append("single from multi")
                    natural.append(df['natural'][i])
                    organism.append(df['organism'][i])
                    protein_refs.append(df['protein_refs'][i])
                    protein_db.append(df['protein_db'][i])
                
                    if df['prob_rev'][i][j]=='r' or df['prob_rev'][i][j]=='p_r':
                        rxn_idx.append(i)
                        orig_rxn_text.append(df['rxn_text'][i])
                        mapped.append(reverse(result))
                        unmapped.append(reverse(unmap(result.split('>')[0])+">>"+unmap(result.split('>')[-1])))
                        rule.append(df['rules'][i][j])
                        rule_id.append(df['rule_ids'][i][j])
                        quality.append(df['quality'][i][j])
                        steps.append("single from multi")    
                    
                        if df['prob_rev'][i][j]=='r':
                            source.append(df['source'][i]+' reversed')
                        else:
                            source.append(df['source'][i]+' reversed_suggested')
                        natural.append(df['natural'][i])
                        organism.append(df['organism'][i])
                        protein_refs.append(df['protein_refs'][i])
                        protein_db.append(df['protein_db'][i])

    if len(set([len(x) for x in [rxn_idx,  mapped, unmapped, orig_rxn_text, rule, rule_id, source, steps, quality]])) != 1:
        raise ValueError('Error in suggesting reversible reactions encountered.')
    df_final = pd.DataFrame(list(zip(rxn_idx, mapped, unmapped, orig_rxn_text, rule, rule_id, source, steps, quality, natural, organism, protein_refs, protein_db)),
                            columns=['rxn_idx', 'mapped', 'unmapped', 'orig_rxn_text', 'rule', 'rule_id', 'source', 'steps', 'quality', 'natural', 'organism', 'protein_refs', 'protein_db'])

    #Standardize mappings:
    mapped=[]
    cached=dict()
    for i in df_final.index:
        print(i,end='\r')
        if df_final['mapped'][i] not in cached.keys():
            standardized = standardize_mapping_rxn(df_final['mapped'][i])
            mapped.append(standardized)
            cached[df_final['mapped'][i]] = standardized
        else:
            mapped.append(cached[df_final['mapped'][i]])
    
    df_final['mapped'] = mapped

    return df_final

def map_isomerase(rxn):
    """map_isomerase.

    Check whether the rxn only contains stereochemical changes

    Args:
        rxn: Reaction SMILES

    Returns:
        Changed reaction SMILES or None
    """
    r,_,p=rxn.split(">")
    if achiral(r)==achiral(p):
        mol_r=Chem.MolFromSmiles(r)
        mol_p=Chem.MolFromSmiles(p)
        [a.SetAtomMapNum(a.GetIdx()+1) for a in mol_r.GetAtoms() if a.GetSymbol()!='H']
        if len(mol_p.GetSubstructMatches(mol_r))>1:
            print("possible different mappings, ignoring all but first",rxn)
        for i,j in enumerate(mol_p.GetSubstructMatch(mol_r)):
            mol_p.GetAtomWithIdx(j).SetAtomMapNum(i+1)           
        r=Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol_r)))
        p=Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol_p)))
        return r+">>"+p
    else:
        return None
