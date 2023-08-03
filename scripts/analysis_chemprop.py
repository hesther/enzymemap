import pandas as pd
import numpy as np
import chemprop

            
def cgr_model(system,n_cpus=48,ii=0):
    df = pd.read_csv("data_analysis/"+system+"_cgr_class_test_split"+str(ii)+".csv")
    for t in ["_1000", "_3000", "_5000", "_10000",""]:
        if t == "_10000" and system=='rhea':
            continue #this datapoint does not exist because rhea is small
        train = "data_analysis/"+system+"_cgr_class_train_split"+str(ii)+str(t)+".csv"
        print(train)
        arguments = [
            '--data_path', train,
            '--separate_val_path',"data_analysis/"+system+"_cgr_class_val_split"+str(ii)+".csv",
            '--separate_test_path',"data_analysis/"+system+"_cgr_class_test_split"+str(ii)+".csv",
            '--dataset_type', 'classification',
            '--save_dir', 'models/cgr-chemprop/'+system+str(ii)+t,
            '--reaction',
            '--smiles_column','smiles',
            '--target_columns','class',
            '--epochs', '100',
            '--num_workers', str(n_cpus),
            '--cache_cutoff', '50000',
            '--save_preds' 
        ]

        args = chemprop.args.TrainArgs().parse_args(arguments)
        mean_score, std_score = chemprop.train.cross_validate(args=args, train_func=chemprop.train.run_training)


        df['preds'] = pd.read_csv('models/cgr-chemprop/'+system+str(ii)+t+"/test_preds.csv")["class"]

        ok = 0
        notok = 0
        indices = list(set(df['index'].values))
        for i in indices:
            true = df[df['index']==i]['class'].values
            pred = df[df['index']==i]['preds'].values
            if true[np.argmax(pred)]==1:
                ok +=1
            else:
                notok +=1

        correct=0
        for i in df.index:
            if df['class'][i] == 1:
                if df['preds'][i] >=0.5:
                    correct+=1
            elif df['class'][i] == 0:
                if df['preds'][i] <0.5:
                    correct+=1
        print(train)
        print("flat accuracy:",correct/len(df))
        with open('models/cgr-chemprop/'+system+str(ii)+t+"/flat_accuracy.csv",'w') as f:
            print(correct/len(df),file=f)
        print("top-1-accuracy",ok,notok, ok/(notok+ok))
        with open('models/cgr-chemprop/'+system+str(ii)+t+"/top_1_accuracy.csv",'w') as f:
            print(ok/(notok+ok),file=f)

if __name__ == '__main__':
    for ii in range(0,10):
        cgr_model('brendadirectsingle',ii=ii)
        cgr_model('metamdb',ii=ii)
        cgr_model('rhea',ii=ii)


