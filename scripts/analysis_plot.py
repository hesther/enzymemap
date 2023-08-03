import numpy as np
import pandas as pd
import ast
import matplotlib.pyplot as plt

def plot_numbers():
    df_sc = pd.read_csv("../data/processed_reactions.csv").drop_duplicates(subset=['mapped','ec_num'])
    raw=pd.read_csv("../../enzymemap/data/raw_reactions.csv").drop_duplicates(subset=['SUBSTRATES','PRODUCTS','EC_NUM'])

    d={}

    for s in range(1,7):
        x=raw[raw['EC_NUM'].str.startswith(str(s))]
        y=df_sc[df_sc['ec_num'].str.startswith(str(s))]
        d[s]={'total':len(x),'balanced':len(x[x['BALANCED_RXNS']!='[]']),'mapped':len(set(y['rxn_idx']))}
    print(d)
    sa=180
    plt.figure(figsize=(8,4))
    nc='lightgray'
    cl='lightgray'
    a1=0.75
    a2=0.5
    

    c = ['#35193e', '#701f57', '#ad1759', '#e13342', '#f37651', '#f6b48f']
    X=[d[x]['total'] for x in range(1,7)]
    plt.pie(X,
            labels=["EC "+str(x)+".x.x.x ("+str(d[x]['total'])+")" for x in range(1,7)],
            colors=c,
            startangle=sa,
            radius=1)

    centre_circle = plt.Circle((0,0),0.67,color='black', fc='white',linewidth=0)
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)
    
    X = [d[x]['balanced'] if y else d[x]['total']-d[x]['balanced'] for x in range(1,7) for y in [True,False]]
    n = plt.pie(X,
                colors=[x if y==True else nc for x in c for y in [True, False]],
                startangle=sa,
                radius=0.65,
                wedgeprops = {"alpha": a1})

    centre_circle = plt.Circle((0,0),0.52,color='black', fc='white',linewidth=0)
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    X = [d[x]['mapped'] if y else d[x]['total']-d[x]['mapped'] for x in range(1,7) for y in [True,False]]
    n = plt.pie(X,
                colors=[x if y==True else nc for x in c for y in [True, False]],
                startangle=sa,
                radius=0.5,
                wedgeprops = {"alpha": a2})

    centre_circle = plt.Circle((0,0),0.37,color='black', fc='white',linewidth=0)
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.text(1.35,0.65,"Overall",fontsize=11)
    plt.text(1.35,0.45,"Balanced",fontsize=11)
    plt.text(1.35,0.25,"Mapped",fontsize=11)
    plt.plot([0.4,2],[0.7,0.7],color=cl,linewidth=1)
    plt.plot([0.35,2],[0.5,0.5],color=cl,linewidth=1)
    plt.plot([0.3,2],[0.3,0.3],color=cl,linewidth=1)

    plt.savefig("composition.png",bbox_inches='tight', dpi=300)
    print("Saved to composition.png")

    
def plot_temprel():
    c = ['#611f53', '#cb1b4f', '#f58860']
    N=[1,3,5,10,50]
    organic_lower=[0.405,0.600,0.623,0.634,0.7]
    organic_higher=[0.552,0.746,0.805,0.869,0.93]
    results={}
    for system in ['brenda','rhea','metamdb']:
        with open("models/ml-fixed/evaluation_corrected_"+system+"/default_retro_accuracy_t.out") as f:
            d=f.readlines()
        results[system]=ast.literal_eval(d[1][0:-1])
        print(system, results)

    plt.figure(figsize=(6,3))
    plt.fill_between(N,organic_lower,organic_higher, color='gray',alpha=0.4,label='organic')
    plt.plot(N,results['rhea'], 'v-', label='RHEA', c=c[0])
    plt.plot(N,results['metamdb'], 's-', label='MetAMDB', c=c[1])
    plt.plot(N,results['brenda'], 'o-', label='EnzymeMap', c=c[2])
    plt.xlabel("$N$")
    plt.ylabel("Top-$N$ accuracy")
    plt.title("Retrosynthesis", fontsize='15')
    plt.xscale('log')
    plt.legend()
    plt.savefig("template_rel.png",bbox_inches='tight', dpi=300)
    print("saved to template_rel.png")

def plot_chemprop():
    c = ['#611f53', '#cb1b4f',  '#f58860','blue']
    regio_flat = {}
    regio_set = {}
    labels={'brendadirect':'EnzymeMap','rhea':'RHEA','metamdb':'MetAMDB'}
    for system in ['rhea','metamdb','brendadirect']:
        f={}
        s={}
        for i in ["_1000","_3000","_5000","_10000",""]:
            if i=="_10000" and system!='metamdb':
                continue
            if i == "":
                n=len(pd.read_csv("data_analysis/"+system+"_cgr_class_train_split0.csv"))
            else:
                n=int(i[1:])
            l_f=[]
            l_s=[]
            for ii in range(10):
                l_f.append(pd.read_csv("models/cgr-chemprop/"+system+str(ii)+i+"/flat_accuracy.csv",header=None)[0][0])
                l_s.append(pd.read_csv("models/cgr-chemprop/"+system+str(ii)+i+"/top_1_accuracy.csv",header=None)[0][0])
            f[n]=(np.mean(l_f),2.262/np.sqrt(10)*np.std(l_f))
            s[n]=(np.mean(l_s),2.262/np.sqrt(10)*np.std(l_s))
        regio_flat[labels[system]]=f
        regio_set[labels[system]]=s
    print(regio_flat)
    print(regio_set)
    fig, axs = plt.subplots(1, 2,figsize=(6,3), sharex=True, sharey=True)

    symbols=['v-','s-','o-','o-']
    for i,system in enumerate(regio_flat.keys()):
        n=list(regio_flat[system].keys())
        axs[0].plot(n,[regio_flat[system][nn][0] for nn in n], symbols[i], label=system, c=c[i])
        axs[0].errorbar(n,[regio_flat[system][nn][0] for nn in n],[regio_flat[system][nn][1] for nn in n], label="", c=c[i],capsize=3,elinewidth=1)

    for i,system in enumerate(regio_set.keys()):
        n=list(regio_set[system].keys())
        axs[1].plot(n,[regio_set[system][nn][0] for nn in n], symbols[i], label=system, c=c[i])
        axs[1].errorbar(n,[regio_set[system][nn][0] for nn in n],[regio_set[system][nn][1] for nn in n],  label="", c=c[i],capsize=3,elinewidth=1)

    axs[0].set_xscale("log")
    axs[0].set_xlabel("# Training reactions")
    axs[1].set_xlabel("# Training reactions")
    axs[0].set_ylabel("Accuracy")
    axs[0].set_title("Flat accuracy")
    axs[1].set_title("Top-1 accuracy in set")
    plt.suptitle("Regioselectivity in forward predictions", fontsize='15')
    fig.subplots_adjust(top=0.82)
    axs[0].legend(loc="lower right")
    plt.savefig("regiosel.png",bbox_inches='tight', dpi=300)
    print("Saved to regiosel.png")
    
if __name__ == '__main__':
    plot_numbers()
    plot_temprel()
    plot_chemprop()
