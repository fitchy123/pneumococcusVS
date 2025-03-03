### Repository for "Artificial intelligence-guided phenotypic virtual screening against drug-resistant Streptococcus Pneumoniae"

A study on the best models for finding potential antibiotics effective against drug-resistant *Streptococcus pneumoniae*.

Dataset for this study is compiled from ChEMBL and PubChem.

#### Repository Structure
This repository contains:
- Pre-processed dataset used in the paper (processed_datasets/10uM_FP_clustered_resistant_pneumococcus_augmented_dataset.csv)
- Instructions for obtaining saved models (README.md)
- Code for evaluating models (models/)
- Instructions for running evaluation (README.md)
- Files for ensuring the correct packages are installed (environment.yml and requirements.txt)

#### Model Evaluation
There are three models in the paper which had results analysed. These are the Random Forest, ChemProp and MolFormer ensemblemodels seen in Figure 4.

The saved models can be downloaded from zenodo: https://zenodo.org/records/14960323

These models should be downloaded and saved in the repository, by default the paths are set to look for them in the "model_checkpoints" folder but you can also pass your own path as a command line argument. 

Below are the commands to run the evaluation for each model:

- Random Forest: `python models/random_forest/run_rf.py --dataset_path processed_datasets/10uM_FP_clustered_resistant_pneumococcus_augmented_dataset.csv --load_model --model_path model_checkpoints/tuned_rfc_morgan_seed3.joblib`
- ChemProp: `python models/chemprop/eval_chemprop.py --data_path processed_datasets/10uM_FP_clustered_resistant_pneumococcus_augmented_dataset.csv --algo chemprop --model_dirs model_checkpoints/_`
- MolFormer: `python models/molformer/lightning_predict.py --data_path "processed_datasets/10uM_FP_clustered_resistant_pneumococcus_augmented_dataset.csv" --agg_fn "median" --checkpoint_path "model_checkpoints/TunedCEMolFormerSeed1RandomVal*" --test_only`

#### Package Installation
- Set up conda environment with environment.yml: `conda env create -f environment.yml`
- Install pip packages: `pip install -r requirements.txt`
