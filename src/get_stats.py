"""
Get statistics of dataset for final results and debugging.
"""

import logging
import pandas as pd

from dataset import Dataset


##### Logging Stats #####
def get_stats_columns() -> tuple[list[str], list[str]]:
    """
    Get the relevant columns for which stats should be calculated
    and a list of descriptions corresponding to the columns.
    """
    df_columns = [
        "parent_molregno",
        "tid",
        "tid_mutation",
        "cpd_target_pair",
        "cpd_target_pair_mutation",
    ]
    columns_descs = [
        "compound ID",
        "target ID",
        "target ID with mutation annotations",
        "compound-target pair",
        "compound-target pair with mutation annotations",
    ]
    return df_columns, columns_descs


def get_stats_for_column(
    df: pd.DataFrame,
    column: str,
    columns_desc: str,
) -> list[list[str, str, int]]:
    """
    Calculate the number of unique values in df[column] and various subsets of df.

    :param df: Pandas Dataframe for which the number of unique values should be calculated
    :type df: pd.DataFrame
    :param column: Column of df that the values should be calculated for
    :type column: str
    :param columns_desc: Description of the column
    :type columns_desc: str
    :return: List of results in the format [column_name, subset_type, size]
    :rtype: list[list[str, str, int]]
    """
    return [
        [column, columns_desc, "all", df[column].nunique()],
        [
            column,
            columns_desc,
            "comparators",
            df[df["DTI"].isin(["DT"])][column].nunique(),
        ],
        [column, columns_desc, "drugs", df[df["DTI"] == "D_DT"][column].nunique()],
        [
            column,
            columns_desc,
            "candidates",
            df[df["DTI"].isin(["C0_DT", "C1_DT", "C2_DT", "C3_DT"])][column].nunique(),
        ],
        [
            column,
            columns_desc,
            "candidates_phase_3",
            df[df["DTI"] == "C3_DT"][column].nunique(),
        ],
        [
            column,
            columns_desc,
            "candidates_phase_2",
            df[df["DTI"] == "C2_DT"][column].nunique(),
        ],
        [
            column,
            columns_desc,
            "candidates_phase_1",
            df[df["DTI"] == "C1_DT"][column].nunique(),
        ],
        [
            column,
            columns_desc,
            "candidates_phase_0",
            df[df["DTI"] == "C0_DT"][column].nunique(),
        ],
    ]


##### Debugging Stats #####
def get_dataset_sizes(df: pd.DataFrame, label: str) -> pd.DataFrame:
    """
    Calculate the number of unique compounds, targets and pairs
    for df and df limited to drugs.

    :param df: Pandas DataFrame for which the dataset sizes should be calculated.
    :type df: pd.DataFrame
    :param label: Description of pipeline step (e.g., initial query).
    :type label: str
    :return: Pandas DataFrame with calculated unique counts.
    :rtype: pd.DataFrame
    """
    stats = {"step": label}

    if "DTI" in df.columns:
        # drugs = compounds of a compound-target pair with a known interaction
        df_drugs = df[df["DTI"] == "D_DT"]
    else:
        df_drugs = df[df["max_phase"] == 4]

    df_columns, _ = get_stats_columns()
    for column in df_columns:
        stats[f"{column}_all"] = df[column].nunique()
        stats[f"{column}_drugs"] = df_drugs[column].nunique()

    df_stats = pd.DataFrame([stats])
    return df_stats


def add_dataset_sizes(
    dataset: Dataset,
    df: pd.DataFrame,
    label: str,
):
    """
    Count and add representative counts of df used for debugging to the dataset.

    :param dataset: Dataset with compound-target pairs and debugging sizes.
    :type dataset: Dataset
    :param df: Pandas DataFrame with current compound-target pairs
    :type df: pd.DataFrame
    :param label: Description of pipeline step (e.g., initial query).
    :type label: str
    """
    df_stats = get_dataset_sizes(df, label)

    dataset.df_sizes_all = pd.concat([dataset.df_sizes_all, df_stats])

    # restrict to data with any pchembl value (any data with a pchembl,
    # even if it is based on only functional data)
    # these statistics are purely based on removing
    # compound-target pairs without pchembl information,
    # i.e., the subset of the dataset is determined by the given df and not recalculated
    df_copy = df.copy()
    df_pchembl = df_copy.dropna(
        subset=[x for x in df_copy.columns if x.startswith("pchembl_value")], how="all"
    )
    df_stats = get_dataset_sizes(df_pchembl, label)
    dataset.df_sizes_pchembl = pd.concat([dataset.df_sizes_pchembl, df_stats])


def add_debugging_info(
    dataset: Dataset,
    df: pd.DataFrame,
    label: str,
):
    """
    Wrapper for add_dataset_sizes.
    Handles logging level.
    """
    if logging.DEBUG >= logging.root.level:
        add_dataset_sizes(dataset, df, label)
