from enzymemap import map_group
from enzymemap.helpers_map import get_groups
import pandas as pd
import ast
import json
import sys

ec_num = sys.argv[1]

df=pd.read_csv("../data/raw_reactions.csv")
for s in ['SUBSTRATES', 'PRODUCTS', 'POSSIBLE_RXNS', 'BALANCED_RXNS']:
    df[s] = [ast.literal_eval(x) for x in df[s]]

rules = pd.read_pickle('../data/rules.pkl')

intermed = df[df['EC_NUM']==ec_num]
df = map_group(intermed, rules)

df.to_csv("../data/processed_reactions_"+ec_num+".csv",index=False)
