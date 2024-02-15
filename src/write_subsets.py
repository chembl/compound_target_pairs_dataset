import logging
import os
import pandas as pd
import sanity_checks

import get_stats
from arguments import OutputArgs, CalculationArgs


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


def get_data_subsets(
    data: pd.DataFrame, min_nof_cpds: int, desc: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Calculate and return the different subsets of interest.

    :param data: Pandas DataFrame with compound-target pairs
    :type data: pd.DataFrame
    :param min_nof_cpds: Miminum number of compounds per target
    :type min_nof_cpds: int
    :param desc: Types of assays current_df contains information about. \
        Options: "BF" (binding+functional), "B" (binding)
    :type desc: str
    :return: 
        - data: Pandas DataFrame with compound-target pairs 
            without the annotations for the opposite desc, \
            e.g. if desc = "BF", the average pchembl value based on 
            binding data only is dropped
        - df_enough_cpds: Pandas DataFrame with targets 
            with at least <min_nof_cpds> compounds with a pchembl value, 
        - df_c_dt_d_dt: As df_enough_cpds but with \
            at least one compound-target pair labelled as 
            'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT' (i.e., known interaction), 
        - df_d_dt: As df_enough_cpds but with \
            at least one compound-target pair labelled as 
            'D_DT' (i.e., known drug-target interaction)
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
        + [  # exclude columns related to the other assay types
            col for col in data.columns if col.startswith("B_") or col.startswith("BF_")
        ]  # exclude filtering columns
    ).drop_duplicates()

    # Restrict the dataset to targets with at least *min_nof_cpds* compounds with a pchembl value.
    comparator_counts = (
        data[data[f"pchembl_value_mean_{desc}"].notnull()]
        .groupby(["tid_mutation"])["parent_molregno"]
        .count()
    )
    # pylint: disable-next=unused-variable
    targets_w_enough_cpds = comparator_counts[
        comparator_counts >= min_nof_cpds
    ].index.tolist()
    df_enough_cpds = data.query("tid_mutation in @targets_w_enough_cpds")

    # Restrict the dataset further to targets
    # with at least one compound-target pair labelled as
    # 'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT',
    # i.e., compound-target pairs with a known interactions.
    # pylint: disable-next=unused-variable
    c_dt_d_dt_targets = set(
        df_enough_cpds[
            df_enough_cpds["DTI"].isin(["D_DT", "C3_DT", "C2_DT", "C1_DT", "C0_DT"])
        ].tid_mutation.to_list()
    )
    df_c_dt_d_dt = df_enough_cpds.query("tid_mutation in @c_dt_d_dt_targets")

    # Restrict the dataset further to targets with
    # at least one compound-target pair labelled as 'D_DT',
    # i.e., known drug-target interactions.
    # pylint: disable-next=unused-variable
    d_dt_targets = set(
        df_enough_cpds[df_enough_cpds["DTI"] == "D_DT"].tid_mutation.to_list()
    )
    df_d_dt = df_enough_cpds.query("tid_mutation in @d_dt_targets")

    return data, df_enough_cpds, df_c_dt_d_dt, df_d_dt


def write_subset_to_file(
    df_combined_subset: pd.DataFrame,
    df_combined: pd.DataFrame,
    desc: str,
    args: CalculationArgs,
    out: OutputArgs,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Write BF or B subsets to file.

    :param df_combined_subset: Subset with binding+functional (BF) or binding (B) assay-based data
        in df_combined
    :type df_combined_subset: pd.DataFrame
    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param desc: Assay description,
        either "BF" (binding+functional) or "B" (binding)
    :type desc: str
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    :return: List of calculated subsets
    :rtype: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]
    """
    (
        df_combined_subset,
        df_combined_subset_enough_cpds,
        df_combined_subset_c_dt_d_dt,
        df_combined_subset_d_dt,
    ) = get_data_subsets(
        df_combined_subset,
        args.min_nof_cpds_bf if desc == "BF" else args.min_nof_cpds_b,
        desc,
    )

    # add filtering columns to df_combined
    for df, col_name in zip(
        [
            df_combined_subset_enough_cpds,
            df_combined_subset_c_dt_d_dt,
            df_combined_subset_d_dt,
        ],
        [
            f"{desc}_{args.min_nof_cpds_bf}",
            f"{desc}_{args.min_nof_cpds_bf}_c_dt_d_dt",
            f"{desc}_{args.min_nof_cpds_bf}_d_dt",
        ],
    ):
        df_combined[col_name] = False
        df_combined.loc[(df_combined.index.isin(df.index)), col_name] = True
        # check that filtering works
        assert df_combined[df_combined[col_name] == True][df.columns].equals(
            df
        ), f"Filtering is not accurate for {col_name}."

    if (desc == "BF" and out.write_bf) or (desc == "B" and out.write_b):
        # NOTE: For BF this is almost identical to the full dataset
        # which will be saved later on.
        # However, the binding-related columns are dropped.
        name_subset = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_CTI_{args.limited_flag}_{desc}",
        )
        write_and_check_output(
            df_combined_subset,
            name_subset,
            desc,
            args,
            out,
        )

        name_subset_100 = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"{desc}_{args.min_nof_cpds_bf}",
        )
        write_and_check_output(
            df_combined_subset_enough_cpds,
            name_subset_100,
            desc,
            args,
            out,
        )

        name_subset_100_c_dt_d_dt = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"{desc}_{args.min_nof_cpds_bf}_c_dt_d_dt",
        )
        write_and_check_output(
            df_combined_subset_c_dt_d_dt,
            name_subset_100_c_dt_d_dt,
            desc,
            args,
            out,
        )

        name_subset_100_d_dt = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"{desc}_{args.min_nof_cpds_bf}_d_dt",
        )
        write_and_check_output(
            df_combined_subset_d_dt, name_subset_100_d_dt, desc, args, out
        )

    return (
        df_combined,
        df_combined_subset,
        df_combined_subset_enough_cpds,
        df_combined_subset_c_dt_d_dt,
        df_combined_subset_d_dt,
    )


def write_bf_to_file(
    df_combined: pd.DataFrame,
    df_sizes: list[list[int], list[int]],
    args: CalculationArgs,
    out: OutputArgs,
) -> pd.DataFrame:
    """
    Calculate relevant subsets for the portion of df_combined
    that is based on binding+functional data.
    If write_bf the subsets are written to output_path.
    Independent of write_bf, filtering columns for BF are added to df_combined and returned.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    :return: Pandas DataFrame with additional filtering columns for BF subsets
    :rtype: pd.Dataframe
    """
    # consider binding and functional assays
    # assay description = binding+functional
    desc = "BF"
    # df_combined without binding only data
    df_combined_subset = df_combined.copy()
    (
        df_combined,
        df_combined_bf,
        df_combined_bf_enough_cpds,
        df_combined_bf_c_dt_d_dt,
        df_combined_bf_d_dt,
    ) = write_subset_to_file(
        df_combined_subset,
        df_combined,
        desc,
        args,
        out,
    )

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined_bf, "binding + functional", df_sizes)
        get_stats.add_dataset_sizes(df_combined_bf_enough_cpds, "BF, >= 100", df_sizes)
        get_stats.add_dataset_sizes(
            df_combined_bf_c_dt_d_dt, "BF, >= 100, c_dt and d_dt", df_sizes
        )
        get_stats.add_dataset_sizes(df_combined_bf_d_dt, "BF, >= 100, d_dt", df_sizes)

    return df_combined


def write_b_to_file(
    df_combined: pd.DataFrame,
    df_sizes: list[list[int], list[int]],
    args: CalculationArgs,
    out: OutputArgs,
) -> pd.DataFrame:
    """
    Calculate relevant subsets for the portion of df_combined that is based on binding data.
    If write_b the subsets are written to output_path.
    Independent of write_b, filtering columns for B are added to df_combined_annotated.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    :return: Pandas DataFrame with additional filtering columns for B subsets
    :rtype: pd.Dataframe
    """
    # consider only binding assays
    # assay description = binding
    desc = "B"
    df_combined_subset = df_combined[df_combined["keep_for_binding"] == True].copy()
    (
        df_combined,
        df_combined_b,
        df_combined_b_enough_cpds,
        df_combined_b_c_dt_d_dt,
        df_combined_b_d_dt,
    ) = write_subset_to_file(
        df_combined_subset,
        df_combined,
        desc,
        args,
        out,
    )

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined_b, "binding", df_sizes)
        get_stats.add_dataset_sizes(df_combined_b_enough_cpds, "B, >= 100", df_sizes)
        get_stats.add_dataset_sizes(
            df_combined_b_c_dt_d_dt, "B, >= 100, c_dt and d_dt", df_sizes
        )
        get_stats.add_dataset_sizes(df_combined_b_d_dt, "B, >= 100, d_dt", df_sizes)

    return df_combined


def write_full_dataset_to_file(
    df_combined: pd.DataFrame,
    args: CalculationArgs,
    out: OutputArgs,
):
    """
    If write_full_dataset, write df_combined with filtering columns to output_path.

    :param df_combined: Pandas DataFrame with compound-target pairs and filtering columns
    :type df_combined: pd.DataFrame
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
        write_and_check_output(df_combined, name_all, desc, args, out)


def output_debug_sizes(
    df_sizes: list[list[int], list[int]],
    out: OutputArgs,
):
    """
    Output counts at various points during calculating the final dataset for debugging.

    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
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
    full_df_sizes = pd.DataFrame(df_sizes[0], columns=column_names)
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
    df_pchembl_sizes = pd.DataFrame(df_sizes[1], columns=column_names)
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


def output_all_stats(
    df_combined_annotated: pd.DataFrame, args: CalculationArgs, out: OutputArgs
):
    """
    Output stats for all datasets and subsets calculated.

    :param df_combined_annotated: Pandas DataFrame with additional filtering columns
    :type df_combined_annotated: pd.DataFrame
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    output_file = os.path.join(
        out.output_path,
        f"ChEMBL{args.chembl_version}_CTI_{args.limited_flag}_full_dataset_stats",
    )

    output_stats(df_combined_annotated, output_file, out)

    if out.write_bf:
        output_file = os.path.join(
            out.output_path,
            f"ChEMBL{args.chembl_version}_"
            f"CTI_{args.limited_flag}_"
            f"BF_{args.min_nof_cpds_bf}_c_dt_d_dt_stats",
        )
        output_stats(
            df_combined_annotated[df_combined_annotated["BF_100_c_dt_d_dt"]],
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
            df_combined_annotated[df_combined_annotated["B_100_c_dt_d_dt"]],
            output_file,
            out,
        )
