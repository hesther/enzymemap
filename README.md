enzymemap
==============================

Python package to atom-map, correct and suggest enzymatic reactions

### Cite us

If you use EnzymeMap, please cite our publication "EnzymeMap: Curation, validation and data-driven prediction of enzymatic reactions" by E. Heid, D. Probst, W. H. Green and G. K. H. Madsen. Check our preprint on [ChemRxiv](https://doi.org/10.26434/chemrxiv-2023-jzw9w).

### News

* August 2023: VERSION 2: A new version of EnzymeMap is released including important bugfixes for isomerase reactions and some reactions containing protons, as well as the addition of protein information. This makes the raw and processed files rather large. For your application, if no protein information is required, you should delete the respective columns and then drop duplicates.


### Database

To simply use the EnzymeMap database, use data/processed_reactions.csv or download EnzymeMap from [Zenodo](https://doi.org/10.5281/zenodo.7841781). Within python (with a valid enzymemap installation) you can also run `enzymemap.get_data()`.

### Installation

Download enzymemap from Github:

```
git clone https://github.com/hesther/enzymemap.git
cd enzymemap
```

Set up a conda environment (or install the packages in `environment.yml` in any other way convenient to you):

```
conda env create -f environment.yml
conda activate enzymemap
```

Install the enzymemap package:

```
pip install -e .
```

### Reproduce our study: Recreate EnzymeMap

Extract BRENDA in the data folder (run `tar -xzvf brenda_2023_1.txt.tar.gz` in the data folder).

Go to the scripts folder and run

```
python make_raw.py
```

to produce `data/raw_reactions.csv`, `data/compound_to_smiles.json` and `ec_nums.csv`. This step processes BRENDA entries and resolves all trivial names to SMILES. You might need to download a new opsin.jar from the internet that is suitable for your system. We also provide the three processed files, so you can continue with the following steps without running `make_inital.py`

Then, for each EC number run `process.py`, for example to process EC number 1.1.3.2:

```
python process.py 1.1.3.2
```

This produces `data/processed_reactions_1.1.3.2.csv`. Run this for all EC numbers (it is best to parallelize this over many cores). You can also run this the individual calculations on different machines. Once all calculations are done, run

```
python concatenate.py
```

to make one dataframe containing all EC numbers. You now have recreated EnzymeMap.

### Reproduce our study: Train and evaluate machine learning models

Run the scripts `analysis_preprocess.py` (process data), `analysis_temprel.py` (train template relevance model, use conda environment from [templatecorr](https://github.com/hesther/templatecorr)), `analysis_chemprop.py`(train CGR-chemprop model, use conda environment from [chemprop](https://github.com/chemprop/chemprop)) and `analysis_plot` (plot results).


### Reproduce our study: Additional benchmarks

Follow the instructions in the additional_benchmarks folder to process KEGG and MetaCyc.


### Copyright

Copyright (c) 2023, Esther Heid


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
