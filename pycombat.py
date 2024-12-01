%load_ext autoreload

autoreload 2

# ---- #

import logging

logging.basicConfig(level=logging.INFO)

# ---- #

import inmoose

import pandas as pd

import matplotlib.pyplot as plt

# ---- #

dataset_1 = pd.read_pickle("data/GSE18520.pickle")

dataset_2 = pd.read_pickle("data/GSE66957.pickle")

dataset_3 = pd.read_pickle("data/GSE69428.pickle")

# ---- #

#dataset_1.to_csv("data/GSE18520.tsv", sep="\t")
#
#dataset_2.to_csv("data/GSE66957.tsv", sep="\t")
#
#dataset_3.to_csv("data/GSE69428.tsv", sep="\t")

# ---- #

df_expression = pd.concat([dataset_1, dataset_2, dataset_3], join="inner", axis=1)

# ---- #

#plt.boxplot(df_expression)
#
#plt.show()

# ---- #

df_corrected = inmoose.pycombat.pycombat_norm(
    df_expression,
    [
        j
        for j, ds in enumerate([dataset_1, dataset_2, dataset_3])
        for _ in range(len(ds.columns))
    ],
)

# ---- #

#plt.boxplot(df_corrected)
#
#plt.show()
