import logging
import os
import pandas as pd
import sanity_checks

import get_stats
from arguments import OutputArgs, CalculationArgs
from dataset import Dataset


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
    :param filename: Filename to write the output to
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


def output_debug_sizes(
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
    full_df_sizes = pd.DataFrame(dataset.df_sizes_all, columns=column_names)
    logging.debug(full_df_sizes)
    name_full_df_sizes = os.path.join(out.output_path, "debug_full_df_sizes")
    write_output(
        full_df_sizes,
        name_full_df_sizes,
        out,
    )

    logging.debug("Size of dataset with any pchembl values at different points.")
    logging.debug(
        "This includes data for which we only have pchembl data \
            for functional assays but not for binding assays."
    )
    df_pchembl_sizes = pd.DataFrame(dataset.df_sizes_pchembl, columns=column_names)
    logging.debug(df_pchembl_sizes)
    name_pchembl_df_sizes = os.path.join(out.output_path, "debug_pchembl_df_sizes")
    write_output(
        full_df_sizes,
        name_pchembl_df_sizes,
        out,
    )


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
        column_stats = get_stats.get_stats_for_column(df, column, columns_desc)
        stats += column_stats
        for colum_stat in column_stats:
            logging.debug(
                "%20s %s",
                colum_stat[2],
                colum_stat[3],
            )

    df_stats = pd.DataFrame(
        stats, columns=["column", "column_description", "subset_type", "counts"]
    )
    write_output(
        df_stats,
        output_file,
        out,
    )


def output_all_stats(dataset: Dataset, args: CalculationArgs, out: OutputArgs):
    """
    Output stats for all datasets and subsets calculated.

    :param dataset: Dataset with compound-target pairs.
    :type dataset: Dataset
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    output_file = os.path.join(
        out.output_path,
        f"ChEMBL{args.chembl_version}_CTI_{args.limited_flag}_full_dataset_stats",
    )

    output_stats(dataset.df_result, output_file, out)

    if out.write_bf:
        output_file = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"BF_{args.min_nof_cpds_bf}_c_dt_d_dt_stats",
        )
        output_stats(
            dataset.df_result[dataset.df_result["BF_100_c_dt_d_dt"]],
            output_file,
            out,
        )

    if out.write_b:
        output_file = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"B_{args.min_nof_cpds_b}_c_dt_d_dt_stats",
        )
        output_stats(
            dataset.df_result[dataset.df_result["B_100_c_dt_d_dt"]],
            output_file,
            out,
        )
