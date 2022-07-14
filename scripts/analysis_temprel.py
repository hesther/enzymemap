import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from joblib import Parallel, delayed
import os
import tensorflow as tf
import sklearn
from functools import partial
from tensorflow.keras import backend as K
from tensorflow.keras.layers import Layer, Dense

def fingerprint_training_dataset(
    smiles, labels, batch_size=256, train=True,
    fp_length=2048, fp_radius=2, fp_use_features=False, fp_use_chirality=True,
    sparse_labels=False, shuffle_buffer=1024, nproc=8, cache=True, precompute=False
):
    smiles_ds = fingerprint_dataset_from_smiles(smiles, fp_length, fp_radius, fp_use_features, fp_use_chirality, nproc, precompute)
    labels_ds = labels_dataset(labels, sparse_labels)
    ds = tf.data.Dataset.zip((smiles_ds, labels_ds))
    ds = ds.shuffle(shuffle_buffer).batch(batch_size)
    if train:
        ds = ds.repeat()
    if cache:
        ds = ds.cache()
    ds = ds.prefetch(buffer_size=batch_size*3)
    return ds

def fingerprint_dataset_from_smiles(smiles, length, radius, useFeatures, useChirality, nproc=8, precompute=False):
    def smiles_tensor_to_fp(smi, length, radius, useFeatures, useChirality):
        smi = smi.numpy().decode('utf-8')
        length = int(length.numpy())
        radius = int(radius.numpy())
        useFeatures = bool(useFeatures.numpy())
        useChirality = bool(useChirality.numpy())
        fp_bit = smiles_to_fingerprint(smi, length, radius, useFeatures, useChirality)
        return np.array(fp_bit)
    def parse_smiles(smi):
        output = tf.py_function(
            smiles_tensor_to_fp, 
            inp=[smi, length, radius, useFeatures, useChirality], 
            Tout=tf.float32
        )
        output.set_shape((length,))
        return output
    if not precompute:
        ds = tf.data.Dataset.from_tensor_slices(smiles)
        ds = ds.map(map_func=parse_smiles, num_parallel_calls=nproc)
    else:
        if nproc!=0:
            fps = Parallel(n_jobs=nproc, verbose=1)(
                delayed(smiles_to_fingerprint)(smi, length, radius, useFeatures, useChirality) for smi in smiles
            )
        else:
            fps = [smiles_to_fingerprint(smi, length, radius, useFeatures, useChirality) for smi in smiles]
        fps = np.array(fps)
        ds = tf.data.Dataset.from_tensor_slices(fps)
    return ds

def sparse_categorical_crossentropy_from_logits(labels, logits):
    return tf.keras.losses.sparse_categorical_crossentropy(labels, logits, from_logits=True)

def top_k(k=1):
    partial_fn = partial(tf.keras.metrics.sparse_top_k_categorical_accuracy, k=k)
    partial_fn.__name__ = 'top_{}'.format(k)
    return partial_fn

def shuffle_arrays(a, b):
    p = np.random.permutation(len(a))
    return a[p], b[p]

def build_model(
    input_shape, output_shape, num_hidden, hidden_size,
    activation='relu', output_activation=None, dropout=0.0, clipnorm=None,
    optimizer=None, learning_rate=0.001, 
    compile_model=True, loss=None, metrics=None
):
    model = tf.keras.models.Sequential()
    model.add(tf.keras.layers.Input(input_shape))
    for _ in range(num_hidden):
        model.add(tf.keras.layers.Dense(hidden_size, activation=activation))
        if dropout:
            model.add(tf.keras.layers.Dropout(dropout))
    model.add(tf.keras.layers.Dense(output_shape, activation=output_activation))
    if optimizer is None or optimizer == 'adam':
        optimizer = tf.keras.optimizers.Adam(learning_rate)
    if clipnorm is not None:
        optimizer.clipnorm = clipnorm
    if compile_model:
        model.compile(
            optimizer=optimizer,
            loss=loss,
            metrics=metrics
        )
    return model

def relevance(**kwargs):
    loss = sparse_categorical_crossentropy_from_logits
    metrics = [
        top_k(k=1),
        top_k(k=3),
        top_k(k=5),
        top_k(k=10),
        top_k(k=50),
    ]
    options = {
        'loss': loss,
        'metrics': metrics
    }
    options.update(kwargs)
    return build_model(**options)

def labels_dataset(labels, sparse=False):
    if not sparse:
        return tf.data.Dataset.from_tensor_slices(labels)
    coo = labels.tocoo()
    indices = np.array([coo.row, coo.col]).T
    labels = tf.SparseTensor(indices, coo.data, coo.shape)
    labels_ds = tf.data.Dataset.from_tensor_slices(labels)
    labels_ds = labels_ds.map(map_func=tf.sparse.to_dense)
    return labels_ds

def smiles_to_fingerprint(smi, length=2048, radius=2, useFeatures=False, useChirality=True):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        raise ValueError('Cannot parse {}'.format(smi))
    fp_bit = AllChem.GetMorganFingerprintAsBitVect(
        mol=mol, radius=radius, nBits = length, 
        useFeatures=useFeatures, useChirality=useChirality
    )
    return np.array(fp_bit)

def top_k(k=1):
    partial_fn = partial(tf.keras.metrics.sparse_top_k_categorical_accuracy, k=k)
    partial_fn.__name__ = 'top_{}'.format(k)
    return partial_fn

def train_model(num_classes,train_smiles,train_labels,valid_smiles,valid_labels,test_smiles,test_labels,model_name1, model_name2,n_cpus,system):
    #Hyperparameters:
    fp_length=2048
    fp_radius=2
    weight_classes=True
    num_hidden=1
    hidden_size=2048
    dropout=0.2
    learning_rate=0.001
    activation='relu'
    batch_size=512
    clipnorm=None
    epochs=25
    early_stopping=3
    nproc=n_cpus
    precompute_fps=True

    train_smiles, train_labels = shuffle_arrays(train_smiles, train_labels)
    valid_smiles, valid_labels = shuffle_arrays(valid_smiles, valid_labels)

    train_ds = fingerprint_training_dataset(
        train_smiles, train_labels, batch_size=batch_size, train=True,
        fp_length=fp_length, fp_radius=fp_radius, nproc=nproc, precompute=precompute_fps
    )
    train_steps = np.ceil(len(train_smiles)/batch_size).astype(int)

    valid_ds = fingerprint_training_dataset(
        valid_smiles, valid_labels, batch_size=batch_size, train=False,
        fp_length=fp_length, fp_radius=fp_radius, nproc=nproc, precompute=precompute_fps
    )
    valid_steps = np.ceil(len(valid_smiles)/batch_size).astype(int)

    #Setup model details
    model = relevance(
        input_shape=(fp_length,),
        output_shape=num_classes,
        num_hidden=num_hidden,
        hidden_size=hidden_size,
        dropout=dropout,
        learning_rate=learning_rate,
        activation=activation,
        clipnorm=clipnorm
    )

    if not os.path.exists('models/'+system+'/training_'+model_name1):
        os.makedirs('models/'+system+'/training_'+model_name1)
    model_output = 'models/'+system+'/training_'+model_name1+'/'+model_name2+'-weights.hdf5'
    history_output = 'models/'+system+'/training_'+model_name1+'/'+model_name2+'-history.json'

    callbacks = []
    if early_stopping:
        callbacks.append(
            tf.keras.callbacks.EarlyStopping(
                patience=early_stopping,
                restore_best_weights=True
            )
        )
    callbacks.append(
        tf.keras.callbacks.ModelCheckpoint(
            model_output, monitor='val_loss', save_weights_only=True
        )
    )
    if weight_classes:
        class_weight = sklearn.utils.class_weight.compute_class_weight(
            'balanced', np.unique(train_labels), train_labels
        )
    else:
        class_weight = sklearn.utils.class_weight.compute_class_weight(
            None, np.unique(train_labels), train_labels
        )

    if nproc !=0 :
        multiproc = True
        nproc = nproc
    else:
        multiproc = False
        nproc = None
    
    length = len(np.unique(train_labels))
    class_weight_dict = {i : class_weight[i] for i in range(length)}

    #Train neural network
    history = model.fit(
        train_ds, epochs=epochs, steps_per_epoch=train_steps,
        validation_data=valid_ds, validation_steps=valid_steps,
        callbacks=callbacks, class_weight=class_weight_dict,
        use_multiprocessing=multiproc, workers=nproc
    )

    pd.DataFrame(history.history).to_json(history_output)



    #Predict on test set:
        
    if nproc!=0:
        test_fps = Parallel(n_jobs=nproc, verbose=1)(
            delayed(smiles_to_fingerprint)(smi, length=fp_length, radius=fp_radius) for smi in test_smiles
        )
    else:
        test_fps = [smiles_to_fingerprint(smi, length=fp_length, radius=fp_radius) for smi in test_smiles]
    test_fps = np.array(test_fps)

    acc = model.evaluate(test_fps, test_labels, batch_size=batch_size)
    print(model.metrics_names[1:])
    print(acc[1:])
    with open('models/'+system+'/evaluation_'+model_name1+'/'+model_name2+'_accuracy_t.out','w') as f:
        print(model.metrics_names[1:],file=f)
        print(acc[1:],file=f)

    return model,batch_size,test_fps

def ml_fixed_model(system,n_cpus=1):
    for name in ["corrected_"]:
        print("###",name)
        if not os.path.exists('models/ml-fixed/evaluation_'+name+system):
            os.makedirs('models/ml-fixed/evaluation_'+name+system)

        choices=["default"]
            
        for  i in choices:
            num_classes = int(np.load("data_analysis/"+system+"num_classes_"+name+"template_"+i+".npy"))
            print("Training model ("+i+" templates, retro) with",num_classes,"classes:")
            data_train=pd.read_csv("data_analysis/"+system+"_"+name+"template_"+i+"_train.csv")
            data_val=pd.read_csv("data_analysis/"+system+"_"+name+"template_"+i+"_val.csv")
            data_test=pd.read_csv("data_analysis/"+system+"_"+name+"template_"+i+"_test.csv")

            train_smiles=data_train['prod_smiles'].values
            train_labels=data_train[name+"template_"+i+"_id"].values
            valid_smiles=data_val['prod_smiles'].values
            valid_labels=data_val[name+"template_"+i+"_id"].values
            test_smiles=data_test['prod_smiles'].values
            test_labels=data_test[name+"template_"+i+"_id"].values

            model,batch_size,test_fps = train_model(num_classes,train_smiles,train_labels,valid_smiles,valid_labels,test_smiles,test_labels,name+system,i+"_retro",n_cpus,'ml-fixed')
            

if __name__ == '__main__':
    ml_fixed_model('brenda')
    ml_fixed_model('brendadirect')
    ml_fixed_model('metamdb')
    ml_fixed_model('rhea')



