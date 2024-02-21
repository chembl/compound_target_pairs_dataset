"""
Add filtering columns for obtaining the different subsets to the dataset.
"""

import logging
import os

import pandas as pd

from arguments import CalculationArgs, OutputArgs
from dataset import Dataset
import get_stats
import output


def get_data_subsets(data: pd.DataFrame, min_nof_cpds: int, desc: str) -> tuple[
    tuple[pd.DataFrame, str],
    tuple[pd.DataFrame, str],
    tuple[pd.DataFrame, str],
    tuple[pd.DataFrame, str],
]:
    """
    Calculate and return the different subsets of interest.

    :param data: Pandas DataFrame with compound-target pairs
    :type data: pd.DataFrame
    :param min_nof_cpds: Miminum number of compounds per target
    :type min_nof_cpds: int
    :param desc: Types of assays current_df contains information about.
        Options: "BF" (binding+functional), "B" (binding)
    :type desc: str
    :return: List of dataset subsets and the string describing them
        - data: Pandas DataFrame with compound-target pairs
            without filtering columns and without
            the annotations for the opposite desc,
            e.g. if desc = "BF", the average pchembl value based on
            binding data only is dropped
        - df_enough_cpds: Pandas DataFrame with targets
            with at least <min_nof_cpds> compounds with a pchembl value,
        - df_c_dt_d_dt: As df_enough_cpds but with
            at least one compound-target pair labelled as
            'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT' (i.e., known interaction),
        - df_d_dt: As df_enough_cpds but with
            at least one compound-target pair labelled as
            'D_DT' (i.e., known drug-target interaction)
    :rtype: tuple[tuple[pd.DataFrame, str],
           tuple[pd.DataFrame, str],
           tuple[pd.DataFrame, str],
           tuple[pd.DataFrame, str]]
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

    return [
        [data, f"{desc}"],
        [df_enough_cpds, f"{desc}_{min_nof_cpds}"],
        [df_c_dt_d_dt, f"{desc}_{min_nof_cpds}_c_dt_d_dt"],
        [df_d_dt, f"{desc}_{min_nof_cpds}_d_dt"],
    ]


def add_subset_filtering_columns(
    df_combined_subset: pd.DataFrame,
    dataset: Dataset,
    desc: str,
    args: CalculationArgs,
    out: OutputArgs,
):
    """
    Add filtering column for binding + functional vs binding

    :param df_combined_subset: Subset with binding+functional (BF) or binding (B) assay-based data
        in df_combined
    :type df_combined_subset: pd.DataFrame
    :param dataset: Dataset with compound-target pairs.
        Will be updated to only include filtering columns.
    :type dataset: Dataset
    :param desc: Assay description,
        either "BF" (binding+functional) or "B" (binding)
    :type desc: str
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    subsets = get_data_subsets(
        df_combined_subset,
        args.min_nof_cpds_bf if desc == "BF" else args.min_nof_cpds_b,
        desc,
    )

    # write subsets if required
    if (desc == "BF" and out.write_bf) or (desc == "B" and out.write_b):
        for [df_subset, subset_desc] in subsets:
            name_subset = os.path.join(
                out.output_path,
                f"ChEMBL{args.chembl_version}_"
                f"CTI_{args.limited_flag}_"
                f"{subset_desc}",
            )
            output.write_and_check_output(
                df_subset,
                name_subset,
                desc,
                args,
                out,
            )

    # add filtering columns to df_combined
    # do not add a filtering column for BF / B (-> [1:])
    for [df, col_name] in subsets[1:]:
        dataset.df_result[col_name] = False
        dataset.df_result.loc[(dataset.df_result.index.isin(df.index)), col_name] = True
        # check that filtering works
        assert dataset.df_result[dataset.df_result[col_name]][df.columns].equals(
            df
        ), f"Filtering is not accurate for {col_name}."

    if logging.DEBUG >= logging.root.level:
        for [df_subset, subset_desc] in subsets:
            get_stats.add_debugging_info(dataset, df_subset, subset_desc)


def add_filtering_columns(
    dataset: Dataset,
    args: CalculationArgs,
    out: OutputArgs,
):
    """
    Add filtering columns to main dataset and save subsets if required.

    :param dataset: Dataset with compound-target pairs.
        Will be updated to only include filtering columns.
    :type dataset: Dataset
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    # consider binding and functional assays
    # assay description = binding+functional
    desc = "BF"
    # df_combined without binding only data
    df_combined_subset = dataset.df_result.copy()
    add_subset_filtering_columns(
        df_combined_subset,
        dataset,
        desc,
        args,
        out,
    )

    # consider only binding assays
    # assay description = binding
    desc = "B"
    df_combined_subset = dataset.df_result[dataset.df_result["keep_for_binding"]].copy()
    add_subset_filtering_columns(
        df_combined_subset,
        dataset,
        desc,
        args,
        out,
    )
