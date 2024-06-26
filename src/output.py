"""
Write the dataset, subsets and related statistics to files 
and to the command line.
"""

import logging
import os
import pandas as pd
import sanity_checks

from arguments import OutputArgs, CalculationArgs
from dataset import Dataset
import get_stats


##### Writing Output #####
def write_output(
    df: pd.DataFrame,
    filename: str,
    out: OutputArgs,
) -> list[str]:
    """
    Write DataFrame df to output file named <filename>.

    :param df: Pandas Dataframe to write to output file.
    :type df: pd.DataFrame
    :param filename: Filename to write the output to
    :type filename: bool
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    :return: Returns list of types of files that was written to (csv and/or xlsx)
    :rtype: list[str]
    """
    file_type_list = []
    if out.write_to_csv:
        df.to_csv(f"{filename}.csv", sep=out.delimiter, index=False)
        file_type_list.append("csv")
    if out.write_to_excel:
        try:
            with pd.ExcelWriter(f"{filename}.xlsx", engine="xlsxwriter") as writer:
                writer.book.use_zip64()
                df.to_excel(writer, index=False)
            file_type_list.append("xlsx")
        except ValueError as e:  # full dataset may be too large to write to excel
            # remove empty file in case of error to avoid confusion
            if os.path.exists(f"{filename}.xlsx"):
                os.remove(f"{filename}.xlsx")
            print(e)
    return file_type_list


def output_stats(
    df: pd.DataFrame,
    output_file: str,
    out: OutputArgs,
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
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    logging.debug("Stats for %s", output_file)
    stats = []
    df_columns, columns_descs = get_stats.get_stats_columns()
    for column, columns_desc in zip(df_columns, columns_descs):
        logging.debug("Stats for column %s:", column)
        column_stats = get_stats.get_stats_for_column(df, column, columns_desc)
        stats += column_stats
        for colum_stat in column_stats:
            logging.debug("%20s %s", colum_stat[2], colum_stat[3])

    df_stats = pd.DataFrame(
        stats, columns=["column", "column_description", "subset_type", "counts"]
    )
    write_output(
        df_stats,
        output_file,
        out,
    )


def write_and_check_output(
    df: pd.DataFrame,
    filename: str,
    assay_type: str,
    args: CalculationArgs,
    out: OutputArgs,
):
    """
    Write df to file and check that writing was successful.

    :param df: Pandas Dataframe to write to output file.
    :type df: pd.DataFrame
    :param filename: Filename to write the output to (should not include the file extension)
    :type filename: bool
    :param assay_type: Types of assays current_df contains information about. \
        Options: "BF" (binding+functional), 
        "B" (binding), 
        "all" (contains both BF and B information)
    :type assay_type: str
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    file_type_list = write_output(df, filename, out)
    sanity_checks.test_equality(
        df, filename, assay_type, file_type_list, args.calculate_rdkit
    )
    output_stats(df, f"{filename}_stats", out)


##### Output Specific Results #####
def write_full_dataset_to_file(
    dataset: Dataset,
    args: CalculationArgs,
    out: OutputArgs,
):
    """
    If write_full_dataset, write df_combined with filtering columns to output_path.

    :param dataset: Dataset with compound-target pairs.
    :type dataset: Dataset
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    desc = "all"
    if out.write_full_dataset:
        name_all = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_CTI_{args.limited_flag}_full_dataset",
        )
        write_and_check_output(dataset.df_result, name_all, desc, args, out)


def write_debug_sizes(
    dataset: Dataset,
    out: OutputArgs,
):
    """
    Output counts at various points during calculating the final dataset for debugging.

    :param dataset: Dataset with compound-target pairs and debugging sizes.
    :type dataset: Dataset
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    # Size of full dataset at different points.
    name_full_df_sizes = os.path.join(out.output_path, "debug_full_df_sizes")
    write_output(
        dataset.df_sizes_all,
        name_full_df_sizes,
        out,
    )

    # Size of dataset with any pchembl values at different points.
    # This includes data for which we only have pchembl data
    # for functional assays but not for binding assays.
    name_pchembl_df_sizes = os.path.join(out.output_path, "debug_pchembl_df_sizes")
    write_output(
        dataset.df_sizes_pchembl,
        name_pchembl_df_sizes,
        out,
    )
