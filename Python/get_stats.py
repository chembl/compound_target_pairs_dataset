import logging
import os
import pandas as pd

import write_subsets


def get_stats_for_column(df: pd.DataFrame, column: str) -> list[list[str, str, int]]:
    """
    Calculate the number of unique values in df[column] and various subsets of df.

    :param df: Pandas Dataframe for which the number of unique values should be calculated
    :type df: pd.DataFrame
    :param column: Column of df that the values should be calculated for 
    :type column: str
    :return: List of results in the format [column_name, subset_type, size]
    :rtype: list[list[str, str, int]]
    """
    return [
    [column, "all", df[column].nunique()], 
    [column, "comparators", df[df['DTI'].isin(['DT'])][column].nunique()], 
    [column, "drugs", df[df['DTI'] == 'D_DT'][column].nunique()], 
    [column, "candidates", df[df['DTI'].isin(['C0_DT', 'C1_DT', 'C2_DT', 'C3_DT'])][column].nunique()], 
    [column, "candidates_phase_3", df[df['DTI'] == 'C3_DT'][column].nunique()],
    [column, "candidates_phase_2", df[df['DTI'] == 'C2_DT'][column].nunique()],
    [column, "candidates_phase_1", df[df['DTI'] == 'C1_DT'][column].nunique()],
    [column, "candidates_phase_0", df[df['DTI'] == 'C0_DT'][column].nunique()]]


def output_stats(df: pd.DataFrame, output_path: str, write_to_csv: bool, write_to_excel: bool, delimiter: str):
    """
    Summarise and output the number of unique values in the following columns:

    - parent_molregno (compound ID)
    - tid (target ID)
    - tid_mutation (target ID + mutation annotations)
    - cpd_target_pair (compound-target pairs)
    - cpd_target_pair_mutation (compound-target pairs including mutation annotations)

    :param df: Pandas Dataframe for which the stats should be calculated
    :type df: pd.DataFrame
    :param output_path: Path to write the dataset stats to
    :type output_path: str
    :param write_to_csv: True if stats should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if stats should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
    """
    columns = ["parent_molregno", "tid", "tid_mutation", "cpd_target_pair", "cpd_target_pair_mutation"]
    stats = []
    for column in columns:
        logging.debug(f"Stats for column {column}:")
        column_stats = get_stats_for_column(df, column)
        stats += column_stats
        for colum_stat in column_stats:
            logging.debug(f"{colum_stat[1] : <40} {colum_stat[2]}")

    df_stats = pd.DataFrame(stats, columns=["column", "subset_type", "size"])
    name_state = os.path.join(output_path, "full_dataset_stats")
    write_subsets.write_output(df_stats, name_state, write_to_csv, write_to_excel, delimiter)

