# SELINA: Single-cell Assignment using Deep Adversarial Domain Adaptation Network with Large-scale References

## Introduction

Selina is a deep-learning based framework for single cell assignment with multiple references. The algorithm consists of three main steps: cell type balancing, pre-training and fine-tuning. The rare cell types in reference datasets are first oversampled using SMOTE(Synthetic Minority Oversampling Technique), and then trained with MADA(Multi-Adversarial Domain Adaptation), a supervised deep learning framework, to obtain a pre-trained model. An autoencoder is subsquently used to fine-tune the parameters of the pre-trained model. Finally, the labels from reference datasets are transferred to the query dataset based on the fully-trained model. Along with the annotation algorithm, we also collect 137 datasets which were uniformly processed and curated to provide users with comprehensive pre-trained models.

## Workflow

![image](https://github.com/pfren1998/SELINA/blob/main/docs/source/_images/Algorithm.png)

## Install

```
conda create -n Selina
conda activate Selina
conda install -c pfren selina
```

## Tutorial

### Preprocess of query data

This step is to fit your data to our pre-trained models. Make the gene match and do some normalization. We support 3 formats of input: plain,h5 and mtx. The plain format means the the expression file is a matrix where each row is a gene and each column is a cell. Below is an example of the command.

```
selina preprocess --format plain --matrix counts.txt --separator tab --gene-idtype symbol --assembly GRCh38 --count-cutoff 1000 --gene-cutoff 500 --cell-cutoff 10 --directory ./ --outprefix query --mode single
```

This step will generate two output files:

- `query_res.rds`: a rds file which stores a seurat object

- `query_{single/cluster}_expr.txt`: input of the prediction step

### Predict

Here you can use our pre-trained model to predict for your query data

```
selina predict --mode single --input query_single_expr.txt --model pre-trained_params.pt --outprefix query --plot True query_res.rds
```

This step will output three files:

- `query_predictions.txt`: predicted cell types for each cell in the query data(default choose the celltype corresponding to the max probablity as the prediction results)

- `query_probability.txt`: probablity of cells predicted as each of the reference cell types

- `query_pred.png`: a umap png file with cell type annotation on it

## Documentation

https://selina.readthedocs.io/en/latest/index.html

## Citation

Pengfei Ren, Xiaoying Shi, Taiwen Li, Chenfei Wang. SELINA: Single-cell Assignment using Deep Adversarial Domain Adaptation Network with Large-scale References

## Contacts

pfren@tongji.edu.cn
