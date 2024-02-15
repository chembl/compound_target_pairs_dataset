import logging
import os

import pandas as pd

from arguments import CalculationArgs, OutputArgs
import get_stats
import write_subsets


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


def add_subset_filtering_columns(
    df_combined_subset: pd.DataFrame,
    df_combined: pd.DataFrame,
    desc: str,
    args: CalculationArgs,
    out: OutputArgs,
    df_sizes,
) -> pd.DataFrame:
    # TODO update documentation
    """
    Add filtering column for binding + functional vs binding

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
    :return: List of calculated subsets
    :rtype: pd.DataFrame
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

    # write subsets if required
    if (desc == "BF" and out.write_bf) or (desc == "B" and out.write_b):
        for df_subset, subset_desc in zip(
            [
                df_combined_subset,
                df_combined_subset_enough_cpds,
                df_combined_subset_c_dt_d_dt,
                df_combined_subset_d_dt,
            ],
            [
                f"{desc}",
                f"{desc}_{args.min_nof_cpds_bf}",
                f"{desc}_{args.min_nof_cpds_bf}_c_dt_d_dt",
                f"{desc}_{args.min_nof_cpds_bf}_d_dt",
            ],
        ):
            name_subset = os.path.join(
                out.output_path,
                f"ChEMBL{args.chembl_version}_"
                f"CTI_{args.limited_flag}_"
                f"{subset_desc}",
            )
            write_subsets.write_and_check_output(
                df_subset,
                name_subset,
                desc,
                args,
                out,
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

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(
            df_combined_subset, "binding + functional", df_sizes
        )
        get_stats.add_dataset_sizes(
            df_combined_subset_enough_cpds, "BF, >= 100", df_sizes
        )
        get_stats.add_dataset_sizes(
            df_combined_subset_c_dt_d_dt, "BF, >= 100, c_dt and d_dt", df_sizes
        )
        get_stats.add_dataset_sizes(
            df_combined_subset_d_dt, "BF, >= 100, d_dt", df_sizes
        )

    return df_combined


def add_filtering_columns(
    df_combined,
    df_sizes,
    args,
    out,
):
    # TODO: documentation
    # consider binding and functional assays
    # assay description = binding+functional
    desc = "BF"
    # df_combined without binding only data
    df_combined_subset = df_combined.copy()
    df_combined = add_subset_filtering_columns(
        df_combined_subset,
        df_combined,
        desc,
        args,
        out,
        df_sizes,
    )

    # consider only binding assays
    # assay description = binding
    desc = "B"
    df_combined_subset = df_combined[df_combined["keep_for_binding"] == True].copy()
    df_combined = add_subset_filtering_columns(
        df_combined_subset,
        df_combined,
        desc,
        args,
        out,
        df_sizes,
    )

    return df_combined
