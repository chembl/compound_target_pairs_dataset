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
import add_filtering_columns


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
    logging.info("get_aggregated_activity_ct_pairs")
    dataset = get_activity_ct_pairs.get_aggregated_activity_ct_pairs(
        chembl_con, args.limit_to_literature
    )
    get_stats.add_debugging_info(dataset, dataset.df_result, "activity ct-pairs")

    logging.info("add_cti_from_drug_mechanisms")
    get_drug_mechanism_ct_pairs.add_drug_mechanism_ct_pairs(dataset, chembl_con)
    get_stats.add_debugging_info(dataset, dataset.df_result, "dm ct-pairs")

    logging.info("add_cti_annotations")
    add_dti_annotations.add_dti_annotations(dataset)
    get_stats.add_debugging_info(dataset, dataset.df_result, "DTI annotations")

    logging.info("add_all_chembl_compound_properties")
    add_chembl_compound_properties.add_all_chembl_compound_properties(
        dataset, chembl_con, args.limit_to_literature
    )
    get_stats.add_debugging_info(dataset, dataset.df_result, "ChEMBL props")

    logging.info("remove_compounds_without_smiles_and_mixtures")
    clean_dataset.remove_compounds_without_smiles_and_mixtures(dataset, chembl_con)
    get_stats.add_debugging_info(dataset, dataset.df_result, "removed smiles")

    logging.info("add_chembl_target_class_annotations")
    add_chembl_target_class_annotations.add_chembl_target_class_annotations(
        dataset,
        chembl_con,
        args,
        out,
    )
    get_stats.add_debugging_info(dataset, dataset.df_result, "tclass annotations")

    if args.calculate_rdkit:
        logging.info("add_rdkit_compound_descriptors")
        add_rdkit_compound_descriptors.add_rdkit_compound_descriptors(dataset)
        get_stats.add_debugging_info(dataset, dataset.df_result, "RDKit props")

    logging.info("clean_dataset")
    clean_dataset.clean_dataset(dataset, args.calculate_rdkit)
    get_stats.add_debugging_info(dataset, dataset.df_result, "clean df")

    logging.info("sanity_checks")
    sanity_checks.sanity_checks(dataset)

    logging.info("add_filtering_columns")
    add_filtering_columns.add_filtering_columns(dataset, args, out)

    logging.info("write_full_dataset_to_file")
    write_subsets.write_full_dataset_to_file(dataset, args, out)

    logging.info("output_stats")
    write_subsets.output_all_stats(dataset, args, out)

    if logging.DEBUG >= logging.root.level:
        write_subsets.write_debug_sizes(dataset, out)
