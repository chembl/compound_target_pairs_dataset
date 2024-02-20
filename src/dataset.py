from dataclasses import dataclass

import pandas as pd


@dataclass
class Dataset:
    """
    df_result:                  Pandas DataFrame with the full dataset
    drug_mechanism_pairs_set:   Set of compound-target pairs in the drug_mechanism table,
                                used for DTI assignments
    drug_mechanism_targets_set: Set of targets in the drug_mechanism table,
                                used for DTI assigments
    df_sizes_all:               List of intermediate sized of the dataset used for debugging
    df_sizes_pchembl:           List of intermediate sized of the dataset used for debugging
    """

    df_result: pd.DataFrame
    drug_mechanism_pairs_set: set
    drug_mechanism_targets_set: set
    df_sizes_all: list[int]
    df_sizes_pchembl: list[int]
