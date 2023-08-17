from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
from rdkit.Chem import AllChem
from copy import deepcopy
from rdkit.Chem.rdchem import ChiralType, BondType, BondDir, BondStereo
from wrapt_timeout_decorator import *
import re
import itertools
from itertools import chain

def get_smi(x):
    """get_smi.

    Returns a standardized SMILES from any SMILES string

    Args:
        x (str): SMILES string

    Returns:
        str: Standardized string (else None if not possible)
    """
    mol = Chem.MolFromSmiles(x)
    if mol:
        try:
            mol = neutralize_atoms(mol)
            mol = remove_isotope(mol)
        except:
            pass
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    except:
        return None

def unmap(smi):
    """unmap.

    Remove atom map from SMILES

    Args:
        smi (str): SMILES string

    Returns:
        str: Unmapped SMILES string
    """
    mol = Chem.MolFromSmiles(smi)
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    return smi

def neutralize_atoms(mol):
    """neutralize_atoms.

    Neutralize molecule if possible

    Args:
        mol (rdkit.Chem.mol): RDKit molecule

    Returns:
        rdkit.Chem.mol: Updated molecule object
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def remove_isotope(mol):
    """remove_isotope.

    Remove isotope information

    Args:
        mol (rdkit.Chem.mol): RDKit molecule

    Returns:
        rdkit.Chem.mol: Updated molecule object
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return Chem.MolFromSmiles(Chem.MolToSmiles(mol))

def count_nonzero_charges(mol):
    """count_nonzero_charges.

    Count nonzero charges in a molecule

    Args:
        mol (rdkit.Chem.mol): RDKit molecule

    Returns:
        int: Number of nonzero charges 
    """
    return sum([a.GetFormalCharge()!=0 for a in mol.GetAtoms()])

def achiral_mol(mol):
    """achiral_mol.

    Remove stereoinformation information

    Args:
        mol (rdkit.Chem.mol): RDKit molecule

    Returns:
        rdkit.Chem.mol: Updated molecule object
    """
    mol_achiral = deepcopy(mol)
    [a.SetChiralTag(ChiralType.CHI_UNSPECIFIED) for a in mol_achiral.GetAtoms()]
    [(b.SetStereo(BondStereo.STEREONONE), b.SetBondDir(BondDir.NONE)) for b in mol_achiral.GetBonds()]
    return Chem.MolToSmiles(mol_achiral)

def count_chiral(mol):
    """count_chiral.

    Count stereocenters in a molecule

    Args:
        mol (rdkit.Chem.mol): RDKit molecule

    Returns:
        int: Number of stereocenters 
    """
    chiral_atoms = sum([a.GetChiralTag()!= ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()])
    chiral_bonds = sum([b.GetStereo()!= BondStereo.STEREONONE for b in mol.GetBonds()])
    return chiral_atoms + chiral_bonds
    
def get_more_chiral(smis):
    """get_more_chiral.

    From a list of smiles, find the one producing the molecule with most stereocenters specified

    Args:
        smis (List(str)): List of SMILES strings

    Returns:
        int: List of SMILES strings with most stereocenters specified
    """
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    achiral_smis = [achiral_mol(mol) for mol in mols]
    chiral_smis = []
    for ac in set(achiral_smis):
        current = [mols[i] for i in range(len(mols)) if achiral_smis[i] == ac]      
        if len(set(current)) > 1:
            counts = [count_chiral(mol) for mol in current]
            indices = [i for i, x in enumerate(counts) if x == max(counts)]
            chiral_smis.extend([Chem.MolToSmiles(current[i]) for i in indices])           
        else:
            chiral_smis.extend([Chem.MolToSmiles(x) for x in current])            
    return chiral_smis 

def get_tautomers(smi):
    """get_tautomers.

    Basic tautomerizations, e.g. phosphate groups or some nitrogen. If a more elaborate tautomerization is required, change this function.

    Args:
        smi (str): SMILES string

    Returns:
        str: Updated SMILES string
    """
    smarts_rules=['[C:1](-[OH1:2])=[N:3]>>[C:1](=[OH0:2])-[N:3]',
                  '[#6H0:1](-[OH1:2]):[#7H0:3]>>[#6H0:1](=[OH0:2]):[#7H1:3]',
                  '[#6H0:1](-[OH1:2]):[#7H0:3]:[#6H0:4]=[NH1:5]>>[#6H0:1](=[OH0:2]):[#7H0:3]:[#6H0:4]-[NH2:5]',
                 '[c:1]1[n:2][c:3][nH:4][c:5]2[n:6][c:7][nH0:8][c:9]1-2>>[c:1]1[n:2][c:3][nH0:4][c:5]2[n:6][c:7][nH:8][c:9]12',
                 '[N:1]-[CH0:2](-[NH2:3])=[NH0:4]>>[N:1]-[CH0:2](=[NH:3])-[NH1:4]',
                  '[O:1]-[PH1:2](=[OH0:3])(=[OH0:4])-[O:5]>>[O:1]-[PH0:2](=[OH0:3])([OH1:4])-[O:5]',
                 ]
    current = Chem.MolFromSmiles(smi)
    for smarts in smarts_rules:
        rxn = AllChem.ReactionFromSmarts(smarts)
        tautomer = current
        while current:
            try:
                current = rxn.RunReactants([current])[0][0]
                current = Chem.MolFromSmiles(Chem.MolToSmiles(current))
                if current:
                    tautomer = current
            except:
                current = None
        current=tautomer
    return Chem.MolToSmiles(tautomer)

def combine_enantiomers_into_racemic(final_outcomes):
    '''
    If two products are identical except for an inverted CW/CCW or an
    opposite cis/trans, then just strip that from the product. Return
    the achiral one instead.
    
    Transform smiles to achiral canonical order first.
    
    This does not look to invert multiple stereocenters at once
    Args:
        final_outcomes: iterable to act upon
    Returns:
        list: modified final_outcomes
    '''

    #Generalize atom order:
    final_outcomes_corrected=[]
    for smiles in list(final_outcomes)[:]:
        #Generalize atom order:
        try:
            m=Chem.MolFromSmiles(smiles)
            m_achiral=Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(achiral(smiles))))
            m_renum = Chem.RenumberAtoms(m,m.GetSubstructMatch(m_achiral))
            smiles = Chem.MolToSmiles(m_renum,canonical=False)
            final_outcomes_corrected.append(smiles)
        except:  
            final_outcomes_corrected.append(smiles)
    final_outcomes = final_outcomes_corrected
        
    for smiles in list(final_outcomes)[:]:       
        # Look for @@ tetrahedral center
        for match in re.finditer(r'@@', smiles):
            smiles_inv = '%s@%s' % (smiles[:match.start()], smiles[match.end():])
            if smiles_inv in final_outcomes:
                if smiles in final_outcomes:
                    final_outcomes.pop(final_outcomes.index(smiles))
                final_outcomes.pop(final_outcomes.index(smiles_inv))
                # Re-parse smiles so that hydrogens can become implicit
                smiles = smiles[:match.start()] + smiles[match.end():]
                outcome = Chem.MolFromSmiles(smiles)
                if outcome is None:
                    raise ValueError('Horrible mistake when fixing duplicate!')
                smiles = '.'.join(sorted(Chem.MolToSmiles(outcome, True).split('.')))
                final_outcomes.append(smiles)

        # Look for // or \\ trans bond
        # where [^=\.] is any non-double bond or period or slash
        for match in chain(re.finditer(r'(\/)([^=\.\\\/]+=[^=\.\\\/]+)(\/)', smiles), 
                re.finditer(r'(\\)([^=\.\\\/]+=[^=\.\\\/]+)(\\)', smiles)):
            # See if cis version is present in list of outcomes
            opposite = {'\\': '/', '/': '\\'}
            smiles_cis1 = '%s%s%s%s%s' % (smiles[:match.start()], 
                match.group(1), match.group(2), opposite[match.group(3)],
                smiles[match.end():]) 
            smiles_cis2 = '%s%s%s%s%s' % (smiles[:match.start()], 
                opposite[match.group(1)], match.group(2), match.group(3),
                smiles[match.end():])
            # Also look for equivalent trans
            smiles_trans2 = '%s%s%s%s%s' % (smiles[:match.start()], 
                opposite[match.group(1)], match.group(2), 
                opposite[match.group(3)], smiles[match.end():])
            # Kind of weird remove conditionals...
            remove = False
            if smiles_cis1 in final_outcomes:
                final_outcomes.pop(final_outcomes.index(smiles_cis1))
                remove = True 
            if smiles_cis2 in final_outcomes:
                final_outcomes.pop(final_outcomes.index(smiles_cis2))
                remove = True 
            if smiles_trans2 in final_outcomes and smiles in final_outcomes:
                final_outcomes.pop(final_outcomes.index(smiles_trans2))
            if remove:
                final_outcomes.pop(final_outcomes.index(smiles))
                smiles = smiles[:match.start()] + match.group(2) + smiles[match.end():]
                outcome = Chem.MolFromSmiles(smiles)
                if outcome is None:
                    raise ValueError('Horrible mistake when fixing duplicate!')
                smiles = '.'.join(sorted(Chem.MolToSmiles(outcome, True).split('.')))
                final_outcomes.append(smiles)
                
    final_outcomes = [Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) for smiles in final_outcomes]
    return final_outcomes

def count_CNOPH(smi):
    """count_CNOPH.

    Count elements C, N, O, P and H

    Args:
        smi (str): SMILES string

    Returns:
        dict: Dictionary with element counts
    """
    counts={'C':0, 'N':0, 'O':0, 'P':0, 'H':0}
    mol=Chem.AddHs(Chem.MolFromSmiles(smi))
    for atom in mol.GetAtoms():
        try:
            counts[atom.GetSymbol()] += 1
        except:
            pass
    return counts

def diff_CNOPH(d1, d2):
    """count_CNOPH.

    Count difference in elements C, N, O, P and H

    Args:
        d1 (dict): Dictionary with element counts
        d2 (dict): Dictionary with element counts

    Returns:
        dict: Dictionary with difference in element counts
    """    
    d3 = {}
    for key in d1.keys():
        d3[key] = d2[key] - d1[key]
    return d3

def get_diff(rxn, reduced_to_oxidized, oxidized_to_reduced):
    """get_diff.

    Try to balance unbalanced reactions by adding missing cofactors.

    Args:
        rxn (str): Reaction SMILES
        reduced_to_oxidized (dict): Dictionary mapping reduced to oxidized cofactors
        oxidized_to_reduced (dict): Dictionary mapping oxidized to reduced cofactors

    Returns:
        str: Updated Reaction SMILES
    """  
    strip_list = list(reduced_to_oxidized.keys()) + list(oxidized_to_reduced.keys()) + ['[H+]']

    reac, _, prod = rxn.split(">")
    reac = reac.split('.')
    prod = prod.split('.')
    
    append_ox = []
    append_red = []
    for key in reduced_to_oxidized.keys():
        if key in reac:
            append_ox = [reduced_to_oxidized[key]]
            append_red = [key, '[H+]']
            break
    if append_ox == []:
        for key in oxidized_to_reduced.keys():
            if key in reac:
                append_ox = [key]
                append_red = [oxidized_to_reduced[key], '[H+]']
                break      
    
    reac = [x for x in reac if x not in strip_list]
    prod = [x for x in prod if x not in strip_list]
    
    d1 = count_CNOPH('.'.join(reac))
    d2 = count_CNOPH('.'.join(prod))

    diff = diff_CNOPH(d1, d2)
    if diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': 2}:
        return '.'.join(reac+append_red) + '>>' + '.'.join(prod+append_ox)
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -2}:
        return '.'.join(reac+append_ox) + '>>' + '.'.join(prod+append_red)
    #Double oxidation
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': 4}:
        return '.'.join(reac+append_red+append_red) + '>>' + '.'.join(prod+append_ox+append_ox)
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -4}:
        return '.'.join(reac+append_ox+append_ox) + '>>' + '.'.join(prod+append_red+append_red)    
    #Triple oxidation
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': 6}:
        return '.'.join(reac+append_red+append_red+append_red) + '>>' + '.'.join(prod+append_ox+append_ox+append_ox)
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -6}:
        return '.'.join(reac+append_ox+append_ox+append_ox) + '>>' + '.'.join(prod+append_red+append_red+append_red)
    
    return rxn

def get_diff_h(rxn): 
    """get_diff_h.

    Try to balance unbalanced reactions by adding protons.

    Args:
        rxn (str): Reaction SMILES

    Returns:
        str: Updated Reaction SMILES
    """
    reac, _, prod = rxn.split(">")
    
    d1 = count_CNOPH(reac)
    d2 = count_CNOPH(prod)

    diff = diff_CNOPH(d1, d2)
    if diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': 1}:
        return reac + ".[H+]" + '>>' + prod
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -1}:
        return reac + '>>' + prod + ".[H+]"
    return rxn

def get_diff_h2o2(rxn): 
    """get_diff_h2o2.

    Try to balance unbalanced reactions by adding hydrogen peroxide (use only with redoxreactions containing H2O2!).

    Args:
        rxn (str): Reaction SMILES

    Returns:
        str: Updated Reaction SMILES
    """
    reac, _, prod = rxn.split(">")
    h2o2 = reac.count('.OO')
    if h2o2 == 0:
        return rxn
    
    reac = reac.split('.')
    prod = prod.split('.')
    reac = [x for x in reac if x not in ['OO','[H+]','O']]
    prod = [x for x in prod if x not in ['OO','[H+]','O']]

    d1 = count_CNOPH('.'.join(reac))
    d2 = count_CNOPH('.'.join(prod))
    diff = diff_CNOPH(d1, d2)
    
    if diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -2}:
        h2o2 = 1
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -4}:
        h2o2 = 2
    elif diff == {'C': 0, 'N': 0, 'O': 0, 'P': 0, 'H': -6}:
        h2o2 = 3 
        
    return '.'.join(reac + h2o2 * ['OO']) + '>>' + '.'.join(prod + h2o2 * 2 * ['O'])
    

def get_balance(rxn):
    """get_balance.

    Compute whether a reaction is balanced (has the same sum formula for reactants and products)

    Args:
        rxn (str): Reaction SMILES

    Returns:
        bool: Boolean whether the reaction is balanced
    """
    reac, _, prod = rxn.split(">")
    r=Chem.MolFromSmiles(reac)
    p=Chem.MolFromSmiles(prod)
    xr=Chem.rdMolDescriptors.CalcMolFormula(r)
    xp=Chem.rdMolDescriptors.CalcMolFormula(p)
    return xr == xp

def find_multiple_bal(rxn):
    """find_multiple_bal.

    Try to balance unbalanced reactions by changing the stoichiometry.

    Args:
        rxn (str): Reaction SMILES

    Returns:
        List[str]: List of possible balanced Reaction SMILES
    """
    balanced=[]
    rs = list(set(rxn.split(">")[0].split('.')))
    nums_r = [Chem.AddHs(Chem.MolFromSmiles(r)).GetNumAtoms() for r in rs]
    ps = list(set(rxn.split(">")[-1].split('.')))
    nums_p = [Chem.AddHs(Chem.MolFromSmiles(p)).GetNumAtoms() for p in ps]
    multi_rs = [[1,2,3,4]]*len(rs)
    multi_rs = list(itertools.product(*multi_rs))
    multi_ps = [[1,2,3,4]]*len(ps)
    multi_ps = list(itertools.product(*multi_ps))
    lowest_num=9999
    accepted_perms=[]
    for perm_r in multi_rs:
        r = '.'.join([item for sublist in [[rs[i]]*perm_r[i] for i in range(len(rs))] for item in sublist])
        num_r = sum([item for sublist in [[nums_r[i]]*perm_r[i] for i in range(len(rs))] for item in sublist])
        for perm_p in multi_ps:
            skip=False
            for perm in accepted_perms:
                if multiple(perm,perm_r+perm_p):
                    skip=True
                    continue
            if skip:
                continue
            p = '.'.join([item for sublist in [[ps[i]]*perm_p[i] for i in range(len(ps))] for item in sublist])
            num_p = sum([item for sublist in [[nums_p[i]]*perm_p[i] for i in range(len(ps))] for item in sublist])
            if num_r == num_p and num_r<lowest_num:
                if get_balance(r+">>"+p):
                    balanced.append(r+">>"+p)
                    accepted_perms.append(perm_r+perm_p)
                    lowest_num = num_r
                    
    balanced=list(set(balanced))
    return balanced

@timeout(10)
def find_multiple_bal_optional_h(rxn):
    """find_multiple_bal_optional_h.

    Try to balance unbalanced reactions by changing the stoichiometry (disregarding hydrogens).

    Args:
        rxn (str): Reaction SMILES

    Returns:
        List[str]: List of possible balanced Reaction SMILES (disregarding hydrogens)
    """
    balanced = find_multiple_bal(rxn)
    if len(balanced) == 0:
        #Try deleting hydrogens
        if '.[H+]' in rxn:
            balanced = find_multiple_bal(rxn.replace('.[H+]',''))
    return balanced

def multiple(f1,f2):
    """multiple.

    Helper function to detect stoichiometries that are the same, eg 1A + 2B -> 1C is the same as 2A + 4B -> 2C

    Args:
        f1 (List(int)): List of integers
        f2 (List(int)): List of integers

    Returns:
        bool: Whether stoichiometries are the same
    """
    ratio=f2[0]/f1[0]
    if not False in [abs(f2[i]/ratio - f1[i])<0.01 for i in range(len(f1))]:
        return True
    return False

def delete_same_mols(rxns):
    """delete_same_mols.

    Delete molecules that occur on reactant and product site simultaneously

    Args:
        rxn (str): Reaction SMILES

    Returns:
        str: Updated reaction SMILES
    """
    balanced = []
    for rxn in rxns:
        rs = rxn.split(">")[0].split(".")
        ps = rxn.split(">")[-1].split(".")
        for r in rs:
            if r in ps:
                rs.pop(rs.index(r))
                ps.pop(ps.index(r))
        if rxn not in balanced:
            balanced.append(rxn)
    return balanced

def get_strip_list(compound_to_smiles):
    """get_strip_list.

    Makes reduced_to_oxidized and oxidized_to_reduced dictionaries needed for balancing via cofactors

    Args:
        compound_to_smiles (dict): Dictionary mapping trivial names to SMILES

    Returns:
        dict, dict: reduced_to_oxidized and oxidized_to_reduced dictionaries 
    """
    pairs=sorted([(x,compound_to_smiles[x]) for x in compound_to_smiles.keys() if 'NAD' in x and compound_to_smiles[x]!= []])

    reduced_to_oxidized = {}
    oxidized_to_reduced = {}
    for name, smiles in pairs:
        if "NADH" in name or "NADPH" in name: #reduced
            name2 = name.replace("NADPH", "NADP+").replace("NADH","NAD+")
            if name2 in [p[0] for p in pairs]:
                for s in smiles:
                    #I manually made sure that compound_to_smiles only has one entry for NAD species!
                    reduced_to_oxidized[compound_to_smiles[name][0]] = compound_to_smiles[name2][0]
        elif "NAD+" in name or "NADP+" in name:
            name2 = name.replace("NADP+", "NADPH").replace("NAD+","NADH")
            if name2 in [p[0] for p in pairs]:
                for s in smiles:
                    #I manually made sure that compound_to_smiles only has one entry for NAD species!
                    oxidized_to_reduced[compound_to_smiles[name][0]] = compound_to_smiles[name2][0]
    return reduced_to_oxidized, oxidized_to_reduced


def correct_reaction(rxns, rxn_text, reduced_to_oxidized, oxidized_to_reduced):
    """correct_reaction.

    Corrects the reactions of a single BRENDA entry, containing one or multiple possible reaction strings.

    Args:
        rxns (List(str)): List of reaction SMILES.
        rxn_text (str): Text of reaction entry
        reduced_to_oxidized (dict): Dictionary mapping reduced to oxidized cofactors
        oxidized_to_reduced (dict): Dictionary mapping oxidized to reduced cofactors

    Returns:
        List(str): List of corrected reactions
    """
    
    corrected_rxns = []
    for rxn in rxns:

        #Combine enantiomers into racemic compound:
        reac, _, prod = rxn.split(">")
        reac = ".".join(combine_enantiomers_into_racemic(reac.split(".")))
        reac = ".".join(combine_enantiomers_into_racemic(reac.split(".")))
        prod = ".".join(combine_enantiomers_into_racemic(prod.split(".")))
        prod = ".".join(combine_enantiomers_into_racemic(prod.split(".")))
        rxn = reac + ">>" + prod
        if get_balance(rxn):
            corrected_rxns.append(rxn)
            continue

        #Correct errors with unbalanced NAD cofactors:
        if 'NAD' in rxn_text:
            rxn = get_diff(rxn, reduced_to_oxidized, oxidized_to_reduced)
            if get_balance(rxn):
                corrected_rxns.append(rxn)
                continue

        #Correct for missing hydrogens:
        rxn = get_diff_h(rxn)
        if get_balance(rxn):
            corrected_rxns.append(rxn)
            continue

        #Correct for wrong H2O2:
        rxn = get_diff_h2o2(rxn)
        if get_balance(rxn):
            corrected_rxns.append(rxn)
            continue

        #Correct for wrong stoichiometry (takes long, therefore limit number of reactions and size of molecules:
        if len(rxns) <= 100 and Chem.MolFromSmiles(rxn.split(">")[0]).GetNumAtoms() <= 500 and Chem.MolFromSmiles(rxn.split(">")[-1]).GetNumAtoms() <= 500:
            try:
                #Can find more than one balanced option:
                balanced = find_multiple_bal_optional_h(rxn)
                balanced = delete_same_mols(balanced)
                if len(balanced)>0:
                    corrected_rxns.extend(balanced)
                    continue
            except TimeoutError:
                pass

    #Remove same mols from reac and prod side
    corrected_rxns = delete_same_mols(corrected_rxns)

    return corrected_rxns

        
def achiral(smi):
    """achiral.

    Get achiral SMILES string

    Args:
        smi (str): SMILES string

    Returns:
        str: Achiral SMILES string
    """    
    mol = Chem.MolFromSmiles(smi)
    [a.SetChiralTag(ChiralType.CHI_UNSPECIFIED) for a in mol.GetAtoms()]
    [(b.SetStereo(BondStereo.STEREONONE), b.SetBondDir(BondDir.NONE)) for b in mol.GetBonds()]
    return Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))          

def chiral(smi):
    """chiral.

    Get standardized chiral SMILES string

    Args:
        smi (str): SMILES string

    Returns:
        str: Standardized SMILES string
    """
    mol = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(mol)   

def bond_edit_stats(rsmi):
    """bond_edit_stats.

    Get number of changed atoms and bonds upon reaction

    Args:
        rsmi (str): reaction SMILES string

    Returns:
        dict: Dictionary containing change counts
    """
    r, s, p = rsmi.split('>')
    r=r.replace('.[H+]','')
    p=p.replace('.[H+]','')
    rmol = Chem.MolFromSmiles(r)    
    pmol = Chem.MolFromSmiles(p)
    
    pbonds = []
    for bond in pmol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a or b:
            pbonds.append(tuple(sorted([a, b])))
    
    rbonds = []
    for bond in rmol.GetBonds():
        a = bond.GetBeginAtom().GetAtomMapNum()
        b = bond.GetEndAtom().GetAtomMapNum()
        if a or b:
            rbonds.append(tuple(sorted([a, b])))
    
    r_changed = set(rbonds) - set(pbonds)
    p_changed = set(pbonds) - set(rbonds)
    
    atoms_changed = set()
    for ch in list(r_changed)+list(p_changed):
        atoms_changed.add(ch[0])
        atoms_changed.add(ch[1])
    atoms_changed -= set([0])
    return {
        'r_bonds': len(r_changed), 
        'p_bonds': len(p_changed), 
        'atoms': len(atoms_changed)
    }

def select_best(rxns, rules, ids, individuals):
    """select_best.

    From a list of possible reactions, select those with the lowest number of edits

    Args:
        rxns (List(str)): List of reaction SMILES string
        rules (List(str)): List of reaction rules
        ids (List(int)): List of rule ids
        individuals (List(List(str))): List of list of reaction SMILES of individual reactions for multistep reactions

    Returns:
        Updated rxns, rules, ids, individuals
    """
    bond_edit = [sum(bond_edit_stats(r).values()) for r in rxns]
    min_bond_edit = min(bond_edit)

    new_rxns = []
    new_rules = []
    new_ids = []
    new_indis = []
    for j in range(len(rxns)):
        if bond_edit[j] == min_bond_edit:
            new_rxns.append(rxns[j])
            new_rules.append(rules[j])
            new_ids.append(ids[j])
            new_indis.append(individuals[j])
    return new_rxns, new_rules, new_ids, new_indis
    

def put_h_last(rxn):
    reac,_,prod = rxn.split(">")
    reac_list = reac.split(".")
    prod_list = prod.split(".")
    reac = '.'.join([r for r in reac_list if r!="[H+]"]+["[H+]"]*reac_list.count("[H+]"))
    prod = '.'.join([r for r in prod_list if r!="[H+]"]+["[H+]"]*prod_list.count("[H+]"))
    return reac+">>"+prod
