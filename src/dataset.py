"""
Dataclass for handling the calculated compound-target pair dataset 
and related data.
"""

from dataclasses import dataclass

import pandas as pd


@dataclass
class Dataset:
    """
    Calculated compound-target pairs dataset (df_results) and related data.
    
    - df_result:                  Pandas DataFrame with the full dataset
    - drug_mechanism_pairs_set:   Set of compound-target pairs in the drug_mechanism table, \
                                used for DTI assignments
    - drug_mechanism_targets_set: Set of targets in the drug_mechanism table, \
                                used for DTI assigments
    - df_sizes_all:               Pandas DataFrame of intermediate sizes of the dataset, \
                                used for debugging
    - df_sizes_pchembl:           Pandas DataFrame of intermediate sizes of the dataset, \
                                restricted to entries with a pchembl value, \
                                used for debugging
    """

    df_result: pd.DataFrame
    drug_mechanism_pairs_set: set
    drug_mechanism_targets_set: set
    df_sizes_all: pd.DataFrame
    df_sizes_pchembl: pd.DataFrame
