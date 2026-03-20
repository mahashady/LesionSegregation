#$ -S /usr/bin/python

import argparse
import glob
import string
import gzip
import re
import os
import subprocess
import statistics
import json
import numpy as np
import random
import pandas as pd
from hmmlearn import hmm
from scipy.stats.stats import pearsonr   
from collections import Counter
from pathlib import Path

def run_HMM(values_vector, HMM_type):
    """
    Runs HMM on the vector of values.
    
    Parameters
    ----------
    values_vector : vector
        Vector of values to run HMM inference
    HMM_type : str
        Could be HMM_as (for WC asymmetry), HMM_multi(for multiallelic sites) or HMM_ploidy (HMMas for clonal and subclonal mutations). 
        This parameter determines initial matrices of emission and transition and number of possible states (3 for HMM_multi and 5 for HMM_as).
    """    

    if HMM_type in ["HMM_as", "HMM_ploidy"]:
        category_order = ["T>N", "A>N"]
    elif HMM_type == "HMM_multi":
        category_order = ["B", "M"]   # example
    else:
        raise ValueError("Unknown HMM type")

    categories = np.asarray(values_vector).flatten()    
    n_categories = len(np.unique(categories))

    category_map = {val: i for i, val in enumerate(category_order)}
    categories_idx = np.array([category_map[c] for c in categories])
    unique_categories = np.array(category_order)

    # create one-hot encoded matrix
    X_onehot = np.zeros((len(categories_idx), len(unique_categories)), dtype=int)
    X_onehot[np.arange(len(categories_idx)), categories_idx] = 1
    print(X_onehot)

    if HMM_type == "HMM_multi":
        hmm_curr = hmm.MultinomialHMM(n_components=3,init_params="s")
        hmm_curr.transmat_ = [[598/600,1/600,1/600],[1/600,598/600,1/600],[1/600,1/600,598/600]]
        hmm_curr.emissionprob_ = [[0.98, 0.02],[0.92, 0.08],[0.84,0.16]]
        state_names = ["M1", "M2", "M3"]
    elif HMM_type == "HMM_as":
        hmm_curr = hmm.MultinomialHMM(n_components=5,init_params="s")
        hmm_curr.transmat_ = [[0.992,0.002,0.002,0.002,0.002],[0.002,0.992,0.002,0.002,0.002],[0.002,0.002,0.992,0.002,0.002],[0.002,0.002,0.002,0.992,0.002],[0.002,0.002,0.002,0.002,0.992]]
        hmm_curr.emissionprob_ = [[0.9, 0.1],[0.7, 0.3],[0.5,0.5],[0.3,0.7],[0.1,0.9]]
        state_names = ["A1", "A2", "A3", "A4", "A5"]
    elif HMM_type == "HMM_ploidy":
        hmm_curr = hmm.MultinomialHMM(n_components=3,init_params="s")
        hmm_curr.transmat_ = [[0.996,0.002,0.002],[0.002,0.996,0.002],[0.002,0.002,0.996]]
        hmm_curr.emissionprob_ = [[0.9, 0.1],[0.5,0.5],[0.1,0.9]]
        state_names = ["A1", "A2", "A3"]
    else:
        print("Incorrect HMM type provided!")

    hmm_curr.fit(X_onehot)
    ll_curr = hmm_curr.score(X_onehot)
    # Predict the hidden states corresponding to observed X.
    Z = hmm_curr.predict(X_onehot)
    # print(Z)
    states = pd.unique(Z)
    # print(states)
    #print("Learned emission probs:")
    #print(hmm_curr.emissionprob_)

    #print("Learned transition matrix:")
    #print(hmm_curr.transmat_)
    Z_named = np.array([state_names[z] for z in Z])
    return states, hmm_curr.emissionprob_, hmm_curr.transmat_, Z_named, unique_categories