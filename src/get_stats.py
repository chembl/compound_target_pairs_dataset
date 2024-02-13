import pandas as pd


##### Debugging Stats #####
def calculate_dataset_sizes(df: pd.DataFrame) -> list[int]:
    """
    Calculate the number of unique compounds, targets and pairs
    for df and df limited to drugs.

    :param df: Pandas DataFrame for which the dataset sizes should be calculated.
    :type df: pd.DataFrame
    :return: List of calculated unique counts.
    :rtype: list[int]
    """
    now_mols = df["parent_molregno"].nunique()
    now_targets = df["tid"].nunique()
    now_targets_mutation = df["tid_mutation"].nunique()
    now_pairs = df["cpd_target_pair"].nunique()
    now_pairs_mutation = df["cpd_target_pair_mutation"].nunique()

    if "DTI" in df.columns:
        # drugs = compounds of a compound-target pair with a known interaction
        df_drugs = df[df["DTI"] == "D_DT"]
    else:
        df_drugs = df[df["max_phase"] == 4]

    now_drugs = df_drugs["parent_molregno"].nunique()
    now_drug_targets = df_drugs["tid"].nunique()
    now_drug_targets_mutation = df_drugs["tid_mutation"].nunique()
    now_drug_pairs = df_drugs["cpd_target_pair"].nunique()
    now_drug_pairs_mutation = df_drugs["cpd_target_pair_mutation"].nunique()

    return [
        now_mols,
        now_drugs,
        now_targets,
        now_drug_targets,
        now_targets_mutation,
        now_drug_targets_mutation,
        now_pairs,
        now_drug_pairs,
        now_pairs_mutation,
        now_drug_pairs_mutation,
    ]


def add_dataset_sizes(
    df: pd.DataFrame, label: str, df_sizes: list[list[int], list[int]]
):
    """
    Count and add representative counts of df to the list df_sizes used for debugging.

    :param df: Pandas DataFrame with current compound-target pairs
    :type df: pd.DataFrame
    :param label: Description of pipeline step (e.g., initial query).
    :type label: str
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    """
    df_copy = df.copy()
    df_sizes[0].append([label] + calculate_dataset_sizes(df_copy))

    # restrict to data with any pchembl value (any data with a pchembl,
    # even if it is based on only functional data)
    # these statistics are purely based on removing
    # compound-target pairs without pchembl information,
    # i.e., the subset of the dataset is determined by the given df and not recalculated
    df_pchembl = df_copy.dropna(
        subset=[x for x in df_copy.columns if x.startswith("pchembl_value")], how="all"
    )
    df_sizes[1].append([label] + calculate_dataset_sizes(df_pchembl))


##### Logging Stats #####
def get_stats_for_column(
    df: pd.DataFrame, column: str, columns_desc: str
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
