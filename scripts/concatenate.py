import pandas as pd

with open('ec_nums.csv') as f:
    lines = [line.rstrip() for line in f]

df = pd.DataFrame()
for line in lines:
    print(line)
    tmp = pd.read_csv("../data/processed_reactions_"+line+".csv")
    tmp['ec_num']=line
    df = pd.concat([df, tmp], ignore_index=True)

df.to_csv("../data/processed_reactions.csv",index=False)
