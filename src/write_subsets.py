import logging
import os
import pandas as pd
import sanity_checks

import get_stats


def write_output(
    df: pd.DataFrame,
    filename: str,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
) -> list[str]:
    """
    Write DataFrame df to output file named <filename>.

    :param df: Pandas Dataframe to write to output file.
    :type df: pd.DataFrame
    :param filename: Filename to write the output to
    :type filename: bool
    :param write_to_csv: True if output should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if output should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
    :return: Returns list of types of files that was written to (csv and/or xlsx)
    :rtype: list[str]
    """
    file_type_list = []
    if write_to_csv:
        df.to_csv(f"{filename}.csv", sep=delimiter, index=False)
        file_type_list.append("csv")
    if write_to_excel:
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
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
    assay_type: str,
    calculate_RDKit: bool,
):
    """
    Write df to file and check that writing was successful.

    :param df: Pandas Dataframe to write to output file.
    :type df: pd.DataFrame
    :param filename: Filename to write the output to
    :type filename: bool
    :param write_to_csv: True if output should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if output should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
    :param assay_type: Types of assays current_df contains information about. \
        Options: "BF" (binding+functional), "B" (binding), "all" (contains both BF and B information)
    :type assay_type: str
    :param calculate_RDKit: If True, current_df contains RDKit-based columns
    :type calculate_RDKit: bool
    """
    file_type_list = write_output(df, filename, write_to_csv, write_to_excel, delimiter)
    sanity_checks.test_equality(
        df, filename, assay_type, file_type_list, calculate_RDKit
    )


def get_data_subsets(
    data: pd.DataFrame, min_nof_cpds: int, desc: str
) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Calculate and return the different subsets of interest.

    :param data: Pandas DataFrame with compound-target pairs
    :type data: pd.DataFrame
    :param min_nof_cpds: Miminum number of compounds per target
    :type min_nof_cpds: int
    :param desc: Types of assays current_df contains information about. \
        Options: "BF" (binding+functional), "B" (binding)
    :type desc: str
    :return: - data: Pandas DataFrame with compound-target pairs without the annotations for the opposite desc, \
            e.g. if desc = "BF", the average pchembl value based on binding data only is dropped
        - df_enough_cpds: Pandas DataFrame with targets with at least <min_nof_cpds> compounds with a pchembl value, 
        - df_c_dt_d_dt: As df_enough_cpds but with \
            at least one compound-target pair labelled as 'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT' (i.e., known interaction), 
        - df_d_dt: As df_enough_cpds but with \
            at least one compound-target pair labelled as 'D_DT' (i.e., known drug-target interaction)
    :rtype: (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame)
    """
    if desc == "B":
        drop_desc = "BF"
    else:
        drop_desc = "B"
    data = data.drop(
        columns=[
            f"pchembl_value_mean_{drop_desc}",
            f"pchembl_value_max_{drop_desc}",
            f"pchembl_value_median_{drop_desc}",
            f"first_publication_cpd_target_pair_{drop_desc}",
            f"first_publication_cpd_target_pair_w_pchembl_{drop_desc}",
            f"LE_{drop_desc}",
            f"BEI_{drop_desc}",
            f"SEI_{drop_desc}",
            f"LLE_{drop_desc}",
        ]
    ).drop_duplicates()

    # Restrict the dataset to targets with at least *min_nof_cpds* compounds with a pchembl value.
    comparator_counts = (
        data[data[f"pchembl_value_mean_{desc}"].notnull()]
        .groupby(["tid_mutation"])["parent_molregno"]
        .count()
    )
    targets_w_enough_cpds = comparator_counts[
        comparator_counts >= min_nof_cpds
    ].index.tolist()
    df_enough_cpds = data.query("tid_mutation in @targets_w_enough_cpds")

    # Restrict the dataset further to targets
    # with at least one compound-target pair labelled as 'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT',
    # i.e., compound-target pairs with a known interactions.
    c_dt_d_dt_targets = set(
        df_enough_cpds[
            df_enough_cpds["DTI"].isin(["D_DT", "C3_DT", "C2_DT", "C1_DT", "C0_DT"])
        ].tid_mutation.to_list()
    )
    df_c_dt_d_dt = df_enough_cpds.query("tid_mutation in @c_dt_d_dt_targets")

    # Restrict the dataset further to targets with at least one compound-target pair labelled as 'D_DT',
    # i.e., known drug-target interactions.
    d_dt_targets = set(
        df_enough_cpds[df_enough_cpds["DTI"] == "D_DT"].tid_mutation.to_list()
    )
    df_d_dt = df_enough_cpds.query("tid_mutation in @d_dt_targets")

    return data, df_enough_cpds, df_c_dt_d_dt, df_d_dt


def write_BF_to_file(
    df_combined: pd.DataFrame,
    chembl_version: str,
    min_nof_cpds_BF: int,
    output_path: str,
    write_BF: bool,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
    limited_flag: str,
    calculate_RDKit: bool,
    df_sizes: list[list[int], list[int]],
) -> pd.DataFrame:
    """
    Calculate relevant subsets for the portion of df_combined that is based on binding+functional data.
    If write_BF the subsets are written to output_path.
    Independent of write_BF, filtering columns for BF are added to df_combined and returned.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_version: Version of ChEMBL for output files
    :type chembl_version: str
    :param min_nof_cpds_BF: Miminum number of compounds per target
    :type min_nof_cpds_BF: int
    :param output_path: Path to write the output to
    :type output_path: str
    :param write_BF: Should the subsets be written to files?
    :type write_BF: bool
    :param write_to_csv: Should the subsets be written to csv?
    :type write_to_csv: bool
    :param write_to_excel: Should the subsets be written to excel?
    :type write_to_excel: bool
    :param delimiter: Delimiter for csv output
    :type delimiter: str
    :param limited_flag: Document suffix indicating whether the dataset was limited to literature sources
    :type limited_flag: str
    :param calculate_RDKit: Does df_combined include RDKit-based columns?
    :type calculate_RDKit: bool
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :return: Pandas DataFrame with additional filtering columns for BF subsets
    :rtype: pd.Dataframe
    """
    # consider binding and functional assays
    # assay description = binding+functional
    desc = "BF"
    # df_combined with additional filtering columns
    df_combined_annotated = df_combined.copy()
    # df_combined without binding only data
    df_combined_BF = df_combined.copy()
    (
        df_combined_BF,
        df_combined_BF_enough_cpds,
        df_combined_BF_c_dt_d_dt,
        df_combined_BF_d_dt,
    ) = get_data_subsets(df_combined_BF, min_nof_cpds_BF, desc)

    # add filtering columns to df_combined_annotated
    for df, col_name in zip(
        [
            df_combined_BF_enough_cpds,
            df_combined_BF_c_dt_d_dt,
            df_combined_BF_d_dt,
        ],
        [
            f"BF_{min_nof_cpds_BF}",
            f"BF_{min_nof_cpds_BF}_c_dt_d_dt",
            f"BF_{min_nof_cpds_BF}_d_dt",
        ],
    ):
        df_combined_annotated[col_name] = False
        df_combined_annotated.loc[
            (df_combined_annotated.index.isin(df.index)), col_name
        ] = True
        # check that filtering works
        assert df_combined_annotated[df_combined_annotated[col_name] == True][
            df.columns
        ].equals(df), f"Filtering is not accurate for {col_name}."

    if write_BF:
        # NOTE: This is almost identical to the full dataset which will be saved later on.
        # However, the binding-related columns are dropped
        name_BF = os.path.join(
            output_path, f"ChEMBL{chembl_version}_CTI_{limited_flag}_BF"
        )
        write_and_check_output(
            df_combined_BF,
            name_BF,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_BF_100 = os.path.join(
            output_path,
            f"ChEMBL{chembl_version}_CTI_{limited_flag}_BF_{min_nof_cpds_BF}",
        )
        write_and_check_output(
            df_combined_BF_enough_cpds,
            name_BF_100,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_BF_100_c_dt_d_dt = os.path.join(
            output_path,
            f"ChEMBL{chembl_version}_CTI_{limited_flag}_BF_{min_nof_cpds_BF}_c_dt_d_dt",
        )
        write_and_check_output(
            df_combined_BF_c_dt_d_dt,
            name_BF_100_c_dt_d_dt,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_BF_100_d_dt = os.path.join(
            output_path,
            f"ChEMBL{chembl_version}_CTI_{limited_flag}_BF_{min_nof_cpds_BF}_d_dt",
        )
        write_and_check_output(
            df_combined_BF_d_dt,
            name_BF_100_d_dt,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined_BF, "binding + functional", df_sizes)
        get_stats.add_dataset_sizes(df_combined_BF_enough_cpds, "BF, >= 100", df_sizes)
        get_stats.add_dataset_sizes(
            df_combined_BF_c_dt_d_dt, "BF, >= 100, c_dt and d_dt", df_sizes
        )
        get_stats.add_dataset_sizes(df_combined_BF_d_dt, "BF, >= 100, d_dt", df_sizes)

    return df_combined_annotated


def write_B_to_file(
    df_combined: pd.DataFrame,
    df_combined_annotated: pd.DataFrame,
    chembl_version: str,
    min_nof_cpds_B: int,
    output_path: str,
    write_B: bool,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
    limited_flag: str,
    calculate_RDKit: bool,
    df_sizes: list[list[int], list[int]],
) -> pd.DataFrame:
    """
    Calculate relevant subsets for the portion of df_combined that is based on binding data.
    If write_B the subsets are written to output_path.
    Independent of write_B, filtering columns for B are added to df_combined_annotated.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param df_combined_annotated: Pandas DataFrame with additional filtering columns
    :type df_combined_annotated: pd.DataFrame
    :param chembl_version: Version of ChEMBL for output files
    :type chembl_version: str
    :param min_nof_cpds_BF: Miminum number of compounds per target
    :type min_nof_cpds_BF: int
    :param output_path: Path to write the output to
    :type output_path: str
    :param write_BF: Should the subsets be written to files?
    :type write_BF: bool
    :param write_to_csv: Should the subsets be written to csv?
    :type write_to_csv: bool
    :param write_to_excel: Should the subsets be written to excel?
    :type write_to_excel: bool
    :param delimiter: Delimiter for csv output
    :type delimiter: str
    :param limited_flag: Document suffix indicating whether the dataset was limited to literature sources
    :type limited_flag: str
    :param calculate_RDKit: Does df_combined include RDKit-based columns?
    :type calculate_RDKit: bool
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :return: Pandas DataFrame with additional filtering columns for B subsets
    :rtype: pd.Dataframe
    """
    # consider only binding assays
    # assay description = binding
    desc = "B"
    df_combined_B = df_combined[df_combined["keep_for_binding"] == True].copy()
    (
        df_combined_B,
        df_combined_B_enough_cpds,
        df_combined_B_c_dt_d_dt,
        df_combined_B_d_dt,
    ) = get_data_subsets(df_combined_B, min_nof_cpds_B, desc)

    # add filtering columns to df_combined_annotated
    for df, col_name in zip(
        [df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt],
        [
            f"B_{min_nof_cpds_B}",
            f"B_{min_nof_cpds_B}_c_dt_d_dt",
            f"B_{min_nof_cpds_B}_d_dt",
        ],
    ):
        df_combined_annotated[col_name] = False
        df_combined_annotated.loc[
            (df_combined_annotated.index.isin(df.index)), col_name
        ] = True
        # check that filtering works
        assert df_combined_annotated[df_combined_annotated[col_name] == True][
            df.columns
        ].equals(df), f"Filtering is not accurate for {col_name}."

    if write_B:
        name_B = os.path.join(
            output_path, f"ChEMBL{chembl_version}_CTI_{limited_flag}_B"
        )
        write_and_check_output(
            df_combined_B,
            name_B,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_B_100 = os.path.join(
            output_path, f"ChEMBL{chembl_version}_CTI_{limited_flag}_B_{min_nof_cpds_B}"
        )
        write_and_check_output(
            df_combined_B_enough_cpds,
            name_B_100,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_B_100_c_dt_d_dt = os.path.join(
            output_path,
            f"ChEMBL{chembl_version}_CTI_{limited_flag}_B_{min_nof_cpds_B}_c_dt_d_dt",
        )
        write_and_check_output(
            df_combined_B_c_dt_d_dt,
            name_B_100_c_dt_d_dt,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

        name_B_100_d_dt = os.path.join(
            output_path,
            f"ChEMBL{chembl_version}_CTI_{limited_flag}_B_{min_nof_cpds_B}_d_dt",
        )
        write_and_check_output(
            df_combined_B_d_dt,
            name_B_100_d_dt,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined_B, "binding", df_sizes)
        get_stats.add_dataset_sizes(df_combined_B_enough_cpds, "B, >= 100", df_sizes)
        get_stats.add_dataset_sizes(
            df_combined_B_c_dt_d_dt, "B, >= 100, c_dt and d_dt", df_sizes
        )
        get_stats.add_dataset_sizes(df_combined_B_d_dt, "B, >= 100, d_dt", df_sizes)

    return df_combined_annotated


def write_full_dataset_to_file(
    df_combined: pd.DataFrame,
    chembl_version: str,
    output_path: str,
    write_full_dataset: bool,
    write_to_csv: bool,
    write_to_excel: bool,
    delimiter: str,
    limited_flag: str,
    calculate_RDKit: bool,
):
    """
    If write_full_dataset, write df_combined with filtering columns to output_path.

    :param df_combined: Pandas DataFrame with compound-target pairs and filtering columns
    :type df_combined: pd.DataFrame
    :param chembl_version: Version of ChEMBL for output files
    :type chembl_version: str
    :param output_path: Path to write the output to
    :type output_path: str
    :param write_full_dataset: Should the subsets be written to files?
    :type write_full_dataset: bool
    :param write_to_csv: Should the subsets be written to csv?
    :type write_to_csv: bool
    :param write_to_excel: Should the subsets be written to excel?
    :type write_to_excel: bool
    :param delimiter: Delimiter for csv output
    :type delimiter: str
    :param limited_flag: Document suffix indicating whether the dataset was limited to literature sources
    :type limited_flag: str
    :param calculate_RDKit: Does df_combined include RDKit-based columns?
    :type calculate_RDKit: bool
    """
    desc = "all"
    if write_full_dataset:
        name_all = os.path.join(
            output_path, f"ChEMBL{chembl_version}_CTI_{limited_flag}_full_dataset"
        )
        write_and_check_output(
            df_combined,
            name_all,
            write_to_csv,
            write_to_excel,
            delimiter,
            desc,
            calculate_RDKit,
        )
