from dataclasses import dataclass

import pandas as pd


@dataclass
class Dataset:
    """
    df_result:                  Pandas DataFrame with the full dataset
    df_sizes_all:               List of intermediate sized of the dataset used for debugging
    df_sizes_pchembl:           List of intermediate sized of the dataset used for debugging
    drug_mechanism_pairs_set:   Set of compound-target pairs in the drug_mechanism table
    drug_mechanism_targets_set: Set of targets in the drug_mechanism table
    df_cpd_props:               Pandas DataFrame with compound properties and
                                structures for all compound ids in ChEMBL
    atc_levels:                 Pandas DataFrame with ATC annotations in ChEMBL
    target_classes_level1:      Pandas DataFrame with mapping from target id to level 1 target class
    target_classes_level2:      Pandas DataFrame with mapping from target id to level 2 target class
    """

    df_result: pd.DataFrame
    df_cpd_props: pd.DataFrame
    atc_levels: pd.DataFrame
    target_classes_level1: pd.DataFrame
    target_classes_level2: pd.DataFrame
    drug_mechanism_pairs_set: set
    drug_mechanism_targets_set: set
    df_sizes_all: list[int]
    df_sizes_pchembl: list[int]
