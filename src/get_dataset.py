import logging
import sqlite3

import get_activity_ct_pairs
import get_drug_mechanism_ct_pairs
import add_dti_annotations
import add_chembl_compound_properties
import clean_dataset
import add_chembl_target_class_annotations
import add_rdkit_compound_descriptors
import sanity_checks
import write_subsets
import get_stats
from arguments import OutputArgs, CalculationArgs


def get_ct_pair_dataset(
    chembl_con: sqlite3.Connection, args: CalculationArgs, out: OutputArgs
):
    """
    Calculate and output the compound-target pair dataset.

    :param chembl_con: Sqlite3 connection to ChEMBL database
    :type chembl_con: sqlite3.Connection
    :param args: Arguments related to how to calculate the dataset
    :type args: CalculationArgs
    :param out: Arguments related to how to output the dataset
    :type out: OutputArgs
    """
    # list with sizes of full dataset and dataset subset with pchembl values for debugging
    df_sizes = [[], []]

    logging.info("get_aggregated_activity_ct_pairs")
    df_combined = get_activity_ct_pairs.get_aggregated_activity_ct_pairs(
        chembl_con, args.limit_to_literature
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "activity ct-pairs", df_sizes)

    logging.info("add_cti_from_drug_mechanisms")
    df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set = (
        get_drug_mechanism_ct_pairs.add_drug_mechanism_ct_pairs(df_combined, chembl_con)
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "dm ct-pairs", df_sizes)

    logging.info("add_cti_annotations")
    df_combined = add_dti_annotations.add_dti_annotations(
        df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "DTI annotations", df_sizes)

    logging.info("add_all_chembl_compound_properties")
    df_combined, df_cpd_props, atc_levels = (
        add_chembl_compound_properties.add_all_chembl_compound_properties(
            df_combined, chembl_con, args.limit_to_literature
        )
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "ChEMBL props", df_sizes)

    logging.info("remove_compounds_without_smiles_and_mixtures")
    df_combined = clean_dataset.remove_compounds_without_smiles_and_mixtures(
        df_combined, chembl_con
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "removed smiles", df_sizes)

    logging.info("add_chembl_target_class_annotations")
    df_combined, target_classes_level1, target_classes_level2 = (
        add_chembl_target_class_annotations.add_chembl_target_class_annotations(
            df_combined,
            chembl_con,
            args,
            out,
        )
    )
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "tclass annotations", df_sizes)

    logging.info("add_rdkit_compound_descriptors")
    if args.calculate_rdkit:
        df_combined = add_rdkit_compound_descriptors.add_rdkit_compound_descriptors(
            df_combined
        )
        if logging.DEBUG >= logging.root.level:
            get_stats.add_dataset_sizes(df_combined, "RDKit props", df_sizes)

    logging.info("clean_dataset")
    df_combined = clean_dataset.clean_dataset(df_combined, args.calculate_rdkit)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "clean df", df_sizes)

    logging.info("sanity_checks")
    sanity_checks.sanity_checks(
        df_combined,
        df_cpd_props,
        atc_levels,
        target_classes_level1,
        target_classes_level2,
        args.calculate_rdkit,
    )

    logging.info("write_BF_to_file")
    df_combined = write_subsets.write_bf_to_file(
        df_combined,
        df_sizes,
        args,
        out,
    )

    logging.info("write_B_to_file")
    df_combined = write_subsets.write_b_to_file(
        df_combined,
        df_sizes,
        args,
        out,
    )

    logging.info("write_full_dataset_to_file")
    write_subsets.write_full_dataset_to_file(
        df_combined,
        args,
        out,
    )

    logging.info("output_stats")
    write_subsets.output_all_stats(df_combined, args, out)

    if logging.DEBUG >= logging.root.level:
        write_subsets.output_debug_sizes(df_sizes, out)
