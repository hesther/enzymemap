import re
from rdkit import Chem
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*')

import requests
import time
import shutil
import urllib.request as request
from contextlib import closing
from lxml import etree as ET
from collections import defaultdict
import logging
import os
import subprocess
from typing import List, Optional

from .helpers_rdkit import get_smi, get_tautomers, count_nonzero_charges, get_more_chiral

#Adappted from https://github.com/samgoldman97/brenda-parser

OPSIN_URL = "https://github.com/dan2097/opsin/releases/download/2.6.0/opsin-cli-2.6.0-jar-with-dependencies.jar"
#Old opsin location: OPSIN_URL = "https://bitbucket.org/dan2097/opsin/downloads/opsin-2.4.0-jar-with-dependencies.jar"
TEMP_OPSIN_INPUT = ".opsin_temp_input.txt"
TEMP_OPSIN_OUTPUT = ".opsin_temp_output.txt"

def query_opsin(mols: List[str],
                opsin_loc: Optional[str] = "opsin.jar") -> dict:
    """query_opsin.

    Download the .jar file for opsin and query it 

    Args:
        mols (List[str]): mols
        opsin_loc (Optional[str]): opsin_loc

    Returns:
        dict: mapping for all recovered outputs
    """

    # If no opsin, download it
    if not os.path.exists(opsin_loc):
        download_ftp(OPSIN_URL, opsin_loc)

    # Output mols to list
    with open(TEMP_OPSIN_INPUT, "w") as fp:
        fp.write("\n".join(mols))

    # Build command
    cmd = [
        "java", "-jar", opsin_loc, "-osmi", TEMP_OPSIN_INPUT, TEMP_OPSIN_OUTPUT
    ]
    # Run opsin
    try:
        cmd_out = subprocess.run(
            cmd,
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.DEVNULL)  
    except:
        logging.warning(f"Was unable to execute opsin command {' '.join(cmd)}")
        return {}

    # Parse results
    mapping = {}
    if os.path.exists(TEMP_OPSIN_OUTPUT):
        with open(TEMP_OPSIN_OUTPUT, "r") as fp:
            total_lines = fp.readlines()
            if len(total_lines) != len(mols):
                raise AssertionError(
                    f"{len(total_lines)} output by Opsin, expected {len(mols)}")
            for enum, line in enumerate(total_lines):
                mol = mols[enum]
                smiles = line.strip()
                # If smiles is not ""
                if smiles:
                    mapping[mol] = smiles

    # Deeete temp opsin files
    if os.path.exists(TEMP_OPSIN_OUTPUT): os.remove(TEMP_OPSIN_OUTPUT)
    if os.path.exists(TEMP_OPSIN_INPUT): os.remove(TEMP_OPSIN_INPUT)

    return mapping

def download_ftp(ftp_link: str, outfile: str):
    """download_ftp.

    Args:
        ftp_link (str): ftp_link
        outfile (str): outfile
    """

    with closing(request.urlopen(ftp_link)) as r:
        with open(outfile, 'wb') as f:
            shutil.copyfileobj(r, f)
            
def query_pubchem(ids: list,
                  query_type: str = "inchi",
                  save_file: str = "pubchem_save.txt",
                  rm_save: bool = True, 
                  encoding : str ="utf8",
                  return_single : bool = False) -> dict:
    """query_pubchem.

    Args:
        ids (list):
        query_type (str): inchi, chebi, or synonym
        save_file (str):
        rm_save (bool): If true, delete the file saved afterward
        encoding (str): Encoding to send request, defaults to utf-8
        return_single (bool): If true, only return the top hit for smiles

    Return:
        dict mapping ids to smiles lists
    """
    # Add options for query_type

    # 60 seconds
    WAIT_TIME = 60
    DOCHEADER = "<!DOCTYPE PCT-Data PUBLIC \"-//NCBI//NCBI PCTools/EN\" \"NCBI_PCTools.dtd\">"
    URL = "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi"
    REQUIRED_HEADINGS = [
        "PCT-Data_input", "PCT-InputData", "PCT-InputData_query", "PCT-Query",
        "PCT-Query_type", "PCT-QueryType", "PCT-QueryType_id-exchange",
        "PCT-QueryIDExchange"
    ]

    QUERY_SUBTREE_NAMES = ["PCT-QueryIDExchange_input", "PCT-QueryUids"]
    OUTPUT_HEADERS = [
        "PCT-QueryIDExchange_operation-type", "PCT-QueryIDExchange_output-type",
        "PCT-QueryIDExchange_output-method", "PCT-QueryIDExchange_compression"
    ]
    OUTPUT_VALUES = ["same", "smiles", "file-pair", "none"]

    # Start building the query tree
    root = ET.Element("PCT-Data")
    cur_pos = root

    for new_heading in REQUIRED_HEADINGS:
        cur_pos = ET.SubElement(cur_pos, new_heading)

    # Navigate down to where we add the inchis
    query_subtree = cur_pos
    for query_subtree_name in QUERY_SUBTREE_NAMES:
        query_subtree = ET.SubElement(query_subtree, query_subtree_name)

    # Now add the things SPECIFIC to inchi
    if query_type == "inchi":
        query_root, query_name = "PCT-QueryUids_inchis", "PCT-QueryUids_inchis_E"
        query_subtree = ET.SubElement(query_subtree, query_root)
        for id_ in ids:
            new_id = ET.SubElement(query_subtree, query_name)
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")

    elif query_type == "synonym":
        query_root, query_name = "PCT-QueryUids_synonyms", "PCT-QueryUids_synonyms_E"
        query_subtree = ET.SubElement(query_subtree, query_root)
        for id_ in ids:
            new_id = ET.SubElement(query_subtree, query_name)
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")

    elif query_type == "chebi":
        for i in ["PCT-QueryUids_source-ids", "PCT-RegistryIDs"]:
            query_subtree = ET.SubElement(query_subtree, i)
        source_id_name = ET.SubElement(query_subtree,
                                       "PCT-RegistryIDs_source-name")
        source_id_name.text = "ChEBI"

        query_subtree = ET.SubElement(query_subtree,
                                      "PCT-RegistryIDs_source-ids")
        for id_ in ids:
            new_id = ET.SubElement(query_subtree,
                                   "PCT-RegistryIDs_source-ids_E")
            # give this the id text
            try:
                new_id.text = id_
            except ValueError:
                logging.warning(f"Couldn't query {id_} due to bad encoding")
    else:
        raise NotImplemented

    # Go back up to to current position holder
    # Add the output specification
    for output_header, output_value in zip(OUTPUT_HEADERS, OUTPUT_VALUES):
        output_xml = ET.SubElement(cur_pos, output_header)
        output_xml.set("value", output_value)

    out_xml = ET.tostring(root,
                          encoding=encoding,
                          method="xml",
                          xml_declaration=True,
                          doctype=DOCHEADER).decode()

    # Post the request!
    resp = requests.post(URL, data=out_xml.encode('utf-8'))

    # Handle response and build a request to check on status
    resp_tree = ET.fromstring(resp.text)
    waiting_id = resp_tree.xpath("//PCT-Waiting_reqid")
    waiting_id = waiting_id[0].text if waiting_id else None

    STATUS_CHECK_HEADERS = [
        "PCT-Data_input", "PCT-InputData", "PCT-InputData_request",
        "PCT-Request"
    ]

    root = ET.Element("PCT-Data")
    cur_pos = root
    for header in STATUS_CHECK_HEADERS:
        cur_pos = ET.SubElement(cur_pos, header)
    req_id = ET.SubElement(cur_pos, "PCT-Request_reqid")
    req_id.text = waiting_id
    req_type = ET.SubElement(cur_pos, "PCT-Request_type")
    req_type.set("value", "status")
    query_xml = ET.tostring(root,
                            encoding=encoding,
                            method="xml",
                            xml_declaration=True,
                            doctype=DOCHEADER).decode()

    download_link = None
    waiting_time = 0

    # TODO: Add stop timeout condition?
    # Repeatedly query to see if the results are done, then sleep for WAITIME
    # in case they aren't
    while not download_link:
        resp = requests.post(URL, data=query_xml.encode('utf-8'))
        resp_tree = ET.fromstring(resp.text)
        download_link = resp_tree.xpath("//PCT-Download-URL_url")
        download_link = download_link[0].text if download_link else None
        time.sleep(WAIT_TIME)
        waiting_time += WAIT_TIME
        logging.warning(f"Waiting time: {waiting_time} seconds")

    # At conclusion, download the ftp file
    download_ftp(download_link, save_file)

    # Also parse this
    ret_dict = defaultdict(lambda: [])
    with open(save_file, "r") as fp:
        for linenum, line in enumerate(fp):
            line = line.strip()
            split_line = line.split("\t")
            if len(split_line) == 2:
                mol, smiles = split_line
                mol = mol.strip()
                smiles = smiles.strip()
                ret_dict[mol].append(smiles)
            else:
                logging.debug(f"No smiles mol found for {line}")

    # If we should only return a single item and not a list
    if return_single: 
        ret_dict = {k : v[0] for k,v in ret_dict.items() if len(v) > 0}

    # Remove temp file
    if os.path.exists(save_file) and rm_save:
        os.remove(save_file)

    return ret_dict

def parse_brenda_ligand_inchi(brenda_ligand_file: str) -> dict:
    """parse_brenda_ligand_inchi.

    Parse the BRENDA substrates file linking trivial names to InChi

    Args:
        brenda_ligand_file (str): Path to file

    Returns:
        dict: mapping between trivial names and InChis
    """
        
    ret = {}
    if brenda_ligand_file: 
        with open(brenda_ligand_file, "rb") as fh:
            for line in fh:
                try:
                    x = line.strip().decode('utf-8').split('\t')
                    name=x[0]
                    inchi=x[4]
                    
                    inchi = None if len(inchi.strip()) == 1 else inchi

                    # only add if chebi or inchi exists
                    if inchi:
                        ret[name.lower()] = inchi
                except UnicodeDecodeError:
                    pass
    return ret

def parse_brenda_ligand_chebi(brenda_ligand_file: str) -> dict:
    """parse_brenda_ligand_chebi.

    Parse the BRENDA substrates file linking trivial names to CHEBI keys

    Args:
        brenda_ligand_file (str): Path to file

    Returns:
        dict: mapping between trivial names and CHEBI keys
    """
    ret = {}
    if brenda_ligand_file: 
        with open(brenda_ligand_file, "rb") as fh:
            for line in fh:
                try:
                    x = line.strip().decode('utf-8').split('\t')
                    name=x[0]
                    chebi=x[5]
                    
                    chebi = None if len(chebi.strip()) == 1 else chebi

                    # only add if chebi or inchi exists
                    if chebi:
                        ret[name.lower()] = chebi
                except UnicodeDecodeError:
                    pass
    return ret

def swap_LS_DR(x):
    """swap_LS_DR.

    Swap L with S and D with R in string.

    Args:
        x (str): Trivial name

    Returns:
        str: Updated trivial name
    """
    SWAP_MAP = {"L": "S", "S": "L", "D": "R", "R": "D"}
    REPLACE_GROUP = r'([\(]{0,1})([LSRD])([\)]{0,1})([ -])'
    return re.sub(REPLACE_GROUP, lambda m: f"{m.group(1)+SWAP_MAP[m.group(2)]+m.group(3)+m.group(4)}", x)

def comma_for_space(x):
    """comma_for_space.

    Swap comma for space in string.

    Args:
        x (str): Trivial name

    Returns:
        str: Updated trivial name
    """
    return x.replace(" ",",")

def inchi_to_smiles_rdkit(compounds, alias, compound_to_inchi):
    """inchi_to_smiles_rdkit.

    Resolve name via InChi key to SMILES string via RDkit

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary
        compound_to_inchi (dict): Dictionary mapping names to InChi

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    for c,a in zip(compounds, alias):
        resolved[a] = []
        inchi = compound_to_inchi.get(c.lower(), None)
        if inchi:
            mol = Chem.MolFromInchi(inchi)
            smiles = None
            if mol:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            if smiles:
                resolved[a] = [smiles]
    return resolved

def inchi_to_smiles_pubchem(compounds, alias, compound_to_inchi):
    """inchi_to_smiles_pubchem.

    Resolve name via InChi key to SMILES string via PubChem Query

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary
        compound_to_inchi (dict): Dictionary mapping names to InChi

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    inchis = [compound_to_inchi.get(c.lower(), None) for c in compounds if compound_to_inchi.get(c.lower(), None) is not None]
    
    inchi_to_smiles = query_pubchem(inchis, query_type='inchi', save_file='tmp.txt')
        
    for c,a in zip(compounds, alias):
        resolved[a] = []
        inchi = compound_to_inchi.get(c.lower(), None)
        if inchi:
            smiles = inchi_to_smiles[inchi]
            if len(smiles) != 0:
                resolved[a] = list(set(smiles))

    return resolved

def name_to_smiles_opsin(compounds, alias, opsin_loc):
    """name_to_smiles_opsin.

    Resolve name via Opsin Query

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary
        opsin_loc (str): Path to Opsin JAR file

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    results = query_opsin([c.lower() for c in compounds], opsin_loc)
    
    for c,a in zip(compounds, alias):
        resolved[a] = []
        if c.lower() in results:
            smiles = results[c.lower()]
            if smiles:
                resolved[a] = [smiles]

    return resolved

def chebi_to_smiles_pubchem(compounds, alias, compound_to_chebi):
    """chebi_to_smiles_pubchem.

    Resolve name via CHEBI key to SMILES string via PubChem Query (direct)

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary
        compound_to_chebi (dict): Dictionary mapping names to CHEBI

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    chebis = [compound_to_chebi.get(c.lower(), None) for c in compounds if compound_to_chebi.get(c.lower(), None) is not None]
    
    chebi_to_smiles = query_pubchem(chebis, query_type='chebi', save_file='tmp.txt')
        
    for c,a in zip(compounds, alias):
        resolved[a] = []
        chebi = compound_to_chebi.get(c.lower(), None)
        if chebi:
            smiles = chebi_to_smiles[chebi]
            if len(smiles) != 0:
                resolved[a] = list(set(smiles))

    return resolved

def chebi_to_smiles_pubchem_syn(compounds, alias, compound_to_chebi):
    """chebi_to_smiles_pubchem_syn.

    Resolve name via CHEBI key to SMILES string via PubChem Query (synonyms)

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary
        compound_to_chebi (dict): Dictionary mapping names to CHEBI

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    chebis = [compound_to_chebi.get(c.lower(), None) for c in compounds if compound_to_chebi.get(c.lower(), None) is not None]
    
    chebi_to_smiles = query_pubchem(chebis, query_type='synonym', save_file='tmp.txt')
        
    for c,a in zip(compounds, alias):
        resolved[a] = []
        chebi = compound_to_chebi.get(c.lower(), None)
        if chebi:
            smiles = chebi_to_smiles[chebi]
            if len(smiles) != 0:
                resolved[a] = list(set(smiles))

    return resolved    

def name_to_smiles_pubchem(compounds, alias):
    """name_to_smiles_pubchem.

    Resolve name to SMILES string via PubChem Query

    Args:
        compounds (List[str]): List of trivial names
        alias (List[str]): List of names used in dictionary

    Returns:
        dict: mapping between trivial names and SMILES
    """
    resolved = {}
    
    results = query_pubchem([c.lower() for c in compounds], query_type='synonym', save_file='tmp.txt')
        
    for c,a in zip(compounds, alias):
        resolved[a] = []
        smiles = results[c.lower()]
        if len(smiles) != 0:
            resolved[a] = list(set(smiles))

    return resolved  

def resolve(compounds, f, name):
    """resolve.

    Resolve a list of trivial names via a specified method

    Args:
        compounds (List[str]): List of trivial names
        f: Lambda function with method to resolve trivial name to SMILES
        name (str): Name of method

    Returns:
        dict: mapping between trivial names and SMILES
    """
    print(name)
    resolved= f(compounds,compounds)
    print("... via identity:", len([c for c in compounds if resolved[c] != []]))

    unresolved = [c for c in compounds if resolved[c] == []]
    resolved_new = f([comma_for_space(c) for c in unresolved],unresolved)
    ctr=0
    for c in resolved_new:
        if resolved[c] == [] and resolved_new[c] != []:
            resolved[c] = resolved_new[c]
            ctr+=1
    print("... via comma for space:", ctr)

    unresolved = [c for c in compounds if resolved[c] == []]
    resolved_new = f([swap_LS_DR(c) for c in unresolved],unresolved)
    ctr=0
    for c in resolved_new:
        if resolved[c] == [] and resolved_new[c] != []:
            resolved[c] = resolved_new[c]
            ctr+=1
    print("... via swap:", ctr)
    
    return resolved

def resolve_all(compound_df, file_loc_inchi, file_loc_chebi):
    """resolve_all.

    Resolve a list of trivial names via different methods

    Args:
        compound_df (pandas.DataFrame): Pandas dataframe of compounds
        file_loc_inchi (str): Path to BRENDA ligands file mapping names to InChis
        file_loc_chebi (str): Path to BRENDA ligands file mapping names to CHEBI keys

    Returns:
        dict: Updated dataframe containing resolved SMILES columns
    """
    compounds = list(compound_df['compound'].values)

    compound_to_inchi=parse_brenda_ligand_inchi(file_loc_inchi)
    compound_to_chebi=parse_brenda_ligand_chebi(file_loc_chebi)

    resolve_methods = {
        'smiles_via_inchi_rdkit': lambda x, y: inchi_to_smiles_rdkit(x, y, compound_to_inchi),
        'smiles_via_inchi_pubchem': lambda x, y: inchi_to_smiles_pubchem(x, y, compound_to_inchi),
        'smiles_via_chebi_pubchem': lambda x, y: chebi_to_smiles_pubchem(x, y, compound_to_chebi),
        'smiles_via_chebi_pubchem_syn': lambda x, y: chebi_to_smiles_pubchem_syn(x, y, compound_to_chebi),
        'smiles_via_name_opsin': lambda x, y: name_to_smiles_opsin(x, y, "opsin.jar"),
        'smiles_via_name_pubchem': lambda x, y: name_to_smiles_pubchem(x, y),
    }

    for name in resolve_methods:
        f = resolve_methods[name]
        resolved = resolve(compounds, f, name)
        compound_df[name] = compound_df['compound'].map(resolved)

    resolved_columns = [x for x in compound_df.keys() if 'smiles_via' in x]
    compound_df['smiles_all'] = [[]]*len(compound_df)
    for i in resolved_columns:
        compound_df['smiles_all'] = compound_df['smiles_all'] + compound_df[i]
    compound_df['smiles_all'] = [list(set(x)) for x in compound_df['smiles_all'].values]

    print("Unresolved:",sum([x == [] for x in compound_df['smiles_all']]))
    print("Resolved:",sum([x != [] for x in compound_df['smiles_all']]))
    print("Resolved via INCHI RDKIT:",sum([x != [] for x in compound_df['smiles_via_inchi_rdkit']]))
    print("Resolved via INCHI PUBCHEM:",sum([x != [] for x in compound_df['smiles_via_inchi_pubchem']]))
    print("Resolved via CHEBI PUBCHEM:",sum([x != [] for x in compound_df['smiles_via_chebi_pubchem']])) 
    print("Resolved via CHEBI PUBCHEM SYN:",sum([x != [] for x in compound_df['smiles_via_chebi_pubchem_syn']]))
    print("Resolved via NAME PUBCHEM:",sum([x != [] for x in compound_df['smiles_via_name_pubchem']]))
    print("Resolved via NAME OPSIN:",sum([x != [] for x in compound_df['smiles_via_name_opsin']]))

    return compound_df

def standardize_compound_df(compound_df):
    """standardize_compound_df.

    Standardize SMILES strings

    Args:
        compound_df (pandas.DataFrame): Pandas dataframe of compounds and their resolved SMILES strings

    Returns:
        dict: Updated dataframe containing updated SMILES column
    """   
    valid = []
    for l in compound_df['smiles_all']:
        l2 = []
        for x in l:
            try:
                l2.append(Chem.MolToSmiles(Chem.MolFromSmiles(x)))
            except:
                pass
        valid.append(list(set(l2)))
    compound_df['smiles_all'] = valid

    neutral = []
    for i,l in enumerate(compound_df['smiles_all']):
        print(i,end='\r')
        if l != []:
            #Standardize, neutralize, canonicalize
            l2 = list(set([get_smi(x) for x in l]))
            l2 = [x for x in l2 if x]
        
            #Tautomerize a few common forms, try to use best SMILES option
            if l2 != []:
                l2 = list(set([get_tautomers(x) for x in l2]))
                if len(l2)>1:
                    #If multiple smiles options, use those with lowest count of molecules
                    counts = [len(smi.split(".")) for smi in l2]
                    indices = [i for i, x in enumerate(counts) if x == min(counts)]
                    l2 = [l2[i] for i in indices]
        
                    #If multiple smiles options, use those with lowest number of nonzero charges
                    counts = [count_nonzero_charges(Chem.MolFromSmiles(smi)) for smi in l2]
                    indices = [i for i, x in enumerate(counts) if x == min(counts)]
                    l2 = [l2[i] for i in indices]
            
                    #If multiple smiles options, use those with larger number of chiral atoms and bonds if achiral mols are the same
                    l2 = get_more_chiral(l2)     
        else:
            l2 = []
        neutral.append(l2)
    compound_df['smiles_neutral'] = neutral


    return compound_df

def manual_corrections_compounds(compound_df):
    """manual_corrections_compounds.

    Manual corrections of compounds that are wrong in BRENDA.

    Args:
        compound_df (pandas.DataFrame): Pandas dataframe of compounds and their resolved SMILES strings

    Returns:
        dict: Updated dataframe 
    """  
    corrected={}
    corrected['acceptor']=[]
    corrected['[oxidized NADPH-hemoprotein reductase]'] = ['*n1c2nc(=O)[nH]c(=O)c-2nc2cc(C)c(C)cc21']
    corrected['[oxidized NADPHhemoprotein reductase]'] = ['*n1c2nc(=O)[nH]c(=O)c-2nc2cc(C)c(C)cc21']
    corrected['oxidized NADPH-hemoprotein reductase'] = ['*n1c2nc(=O)[nH]c(=O)c-2nc2cc(C)c(C)cc21']
    corrected['oxidized NADPH-hemoprotein reductase]'] = ['*n1c2nc(=O)[nH]c(=O)c-2nc2cc(C)c(C)cc21']
    corrected['[oxidized NADH-hemoprotein reductase]'] = ['*n1c2nc(=O)[nH]c(=O)c-2nc2cc(C)c(C)cc21']
    corrected['[reduced NADPH-hemoprotein reductase]'] = ['*N1c2cc(C)c(C)cc2Nc2c1[nH]c(=O)[nH]c2=O']
    corrected['reduced NADPH-hemoprotein reductase'] = ['*N1c2cc(C)c(C)cc2Nc2c1[nH]c(=O)[nH]c2=O']
    corrected['reduced NADPH-hemoprotein reductase]'] = ['*N1c2cc(C)c(C)cc2Nc2c1[nH]c(=O)[nH]c2=O']
    corrected['[reduced NADH-hemoprotein reductase]'] = ['*N1c2cc(C)c(C)cc2Nc2c1[nH]c(=O)[nH]c2=O']
    corrected['oxidized ferredoxin'] = ['CS[Fe+3]1(SC)S[Fe+3](SC)(SC)S1']
    corrected['reduced ferredoxin'] = ['CS[Fe+2]1(SC)S[Fe+2](SC)(SC)S1']
    corrected['oxidized ferredoxin [iron-sulfur] cluster'] = ['CS[Fe+3]1(SC)S[Fe+3](SC)(SC)S1']
    corrected['reduced ferredoxin [iron-sulfur] cluster'] = ['CS[Fe+2]1(SC)S[Fe+2](SC)(SC)S1']
    corrected['thioredoxin disulfide'] = ['CC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H]1CSSC[C@H](NC(=O)[C@@H](C)Cc2c[nH]c3ccccc23)C(=O)NCC(=O)N2CCC[C@H]2C(=O)N1']
    corrected['thioredoxin'] = ['CC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CS)NC(=O)[C@@H]1CCCN1C(=O)CNC(=O)[C@H](CS)NC(=O)[C@@H](C)Cc1c[nH]c2ccccc12']
    corrected['phosphatidylcholine'] = ['CCCCC(=O)OC[C@@H](COP(=O)(O)OCC[N+](C)(C)C)OC(=O)CCCC']
    corrected['(S)-NADH-hydrate'] = ['NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['2-azido-NADH'] = ['[N-]=[N+]=Nc1nc(N)c2ncn([C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]4O[C@@H](N5C=CCC(C(N)=O)=C5)[C@H](O)[C@@H]4O)[C@@H](O)[C@H]3O)c2n1']
    corrected["3'-NADPH"] = ['NC(=O)C1=CN([C@@H]2O[C@@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3OP(=O)(O)O)[C@H](O)[C@@H]2O)C=CC1']
    corrected['NADHX'] = ['NC(=O)C1CCC(O)N([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](c4nc5c(N)ncnc5[nH]4)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C1']
    corrected['alpha-NADP+'] = ['NC(=O)c1ccc[n+]([C@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['alpha-NADPH'] = ['NC(=O)C1=CN([C@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['beta-thio-NADP+'] = ['NC(=S)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['beta-thio-NADPH'] = ['NC(=S)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['thio-NAD+'] = ['NC(=S)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['thio-NADH'] = ['NC(=S)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['thio-NADP+'] = ['NC(=S)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['thio-NADPH'] = ['NC(=S)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['thionicotinamide-NAD+'] = ['NC(=S)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['thionicotinamide-NADH'] = ['NC(=S)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['thionicotinamide-NADP+'] = ['NC(=S)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1']
    corrected['thionicotinamide-NADPH'] = ['NC(=S)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1']
    corrected['1,N6-ethanoadenine-NADPH'] = ['NC(=O)C1=CN([C@H]2O[C@@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c4N=CN4CCN=C54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@H](O)[C@@H]2O)C=CC1']

    for c in corrected.keys():
        index = compound_df[compound_df['compound']==c].index[0]
        compound_df.loc[index]['smiles_neutral']=corrected[c]

    return compound_df
