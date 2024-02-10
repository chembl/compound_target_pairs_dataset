import logging
import os

import pandas as pd

import write_subsets


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


def output_debug_sizes(
    df_sizes: list[list[int], list[int]],
    output_path: str,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
):
    """
    Output counts at various points during calculating the final dataset for debugging.

    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :param output_path: Path to write the dataset counts to
    :type output_path: str
    :param write_to_csv: True if counts should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if counts should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
    """
    column_names = [
        "type",
        "#mols",
        "#drugs",
        "#targets",
        "#drug_ targets",
        "#targets_ mutation",
        "#drug_ targets_mutation",
        "#cpd_tid_ pairs",
        "#drug_tid_ pairs",
        "#cpd_ tid_mutation_ pairs",
        "#drug_ tid_mutation_ pairs",
    ]

    logging.debug("Size of full dataset at different points.")
    full_df_sizes = pd.DataFrame(df_sizes[0], columns=column_names)
    logging.debug(full_df_sizes)
    name_full_df_sizes = os.path.join(output_path, "debug_full_df_sizes")
    write_subsets.write_output(
        full_df_sizes, name_full_df_sizes, write_to_csv, write_to_excel, delimiter
    )

    logging.debug("Size of dataset with any pchembl values at different points.")
    logging.debug(
        "This includes data for which we only have pchembl data \
            for functional assays but not for binding assays."
    )
    df_pchembl_sizes = pd.DataFrame(df_sizes[1], columns=column_names)
    logging.debug(df_pchembl_sizes)
    name_pchembl_df_sizes = os.path.join(output_path, "debug_pchembl_df_sizes")
    write_subsets.write_output(
        full_df_sizes, name_pchembl_df_sizes, write_to_csv, write_to_excel, delimiter
    )


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


def output_stats(
    df: pd.DataFrame,
    output_file: str,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
):
    """
    Summarise and output the number of unique values in the following columns:

    - parent_molregno (compound ID)
    - tid (target ID)
    - tid_mutation (target ID + mutation annotations)
    - cpd_target_pair (compound-target pairs)
    - cpd_target_pair_mutation (compound-target pairs including mutation annotations)

    :param df: Pandas Dataframe for which the stats should be calculated
    :type df: pd.DataFrame
    :param output_file: Path and filename to write the dataset stats to
    :type output_file: str
    :param write_to_csv: True if stats should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if stats should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
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

    logging.debug("Stats for %s", output_file)
    stats = []
    for column, columns_desc in zip(df_columns, columns_descs):
        logging.debug("Stats for column %s:", column)
        column_stats = get_stats_for_column(df, column, columns_desc)
        stats += column_stats
        for colum_stat in column_stats:
            # TODO remove
            logging.debug(
                "%20s %20s %20s %s",
                colum_stat[0],
                colum_stat[1],
                colum_stat[2],
                colum_stat[3],
            )
            # logging.debug(f"{colum_stat[1] : <40} {colum_stat[3]}")

    df_stats = pd.DataFrame(
        stats, columns=["column", "column_description", "subset_type", "counts"]
    )
    write_subsets.write_output(
        df_stats, output_file, write_to_csv, write_to_excel, delimiter
    )
