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


def get_ct_pair_dataset(chembl_con: sqlite3.Connection,
                        chembl_version: str,
                        output_path: str,
                        limit_to_literature: bool,
                        calculate_RDKit: bool,
                        write_to_csv: bool, 
                        write_to_excel: bool,
                        delimiter: str,
                        write_full_dataset: bool, 
                        write_BF: bool, 
                        write_B: bool):
    """
    Calculate and output the compound-target pair dataset.

    :param chembl_con: Sqlite3 connection to ChEMBL database
    :type chembl_con: sqlite3.Connection
    :param chembl_version: Version of ChEMBL for output file names
    :type chembl_version: str
    :param output_path: Path to write output files to
    :type output_path: str
    :param limit_to_literature: Include only literature sources if True. Include all available sources otherwise.
    :type limit_to_literature: bool
    :param calculate_RDKit: True if RDKit-based compound properties should be calculated
    :type calculate_RDKit: bool
    :param write_to_csv: True if output should be written to csv
    :type write_to_csv: bool
    :param write_to_excel: True if output should be written to excel
    :type write_to_excel: bool
    :param delimiter: Delimiter in csv-output
    :type delimiter: str
    :param write_full_dataset: True if the full dataset should be written to output
    :type write_full_dataset: bool
    :param write_BF: True if subsets based on binding+functional data should be written to output
    :type write_BF: bool
    :param write_B: True if subsets based on binding data only should be written to output
    :type write_B: bool
    """
    # list with sizes of full dataset and dataset subset with pchembl values for debugging
    df_sizes = [[], []]

    logging.info("get_aggregated_acticity_ct_pairs")
    df_combined = get_activity_ct_pairs.get_aggregated_acticity_ct_pairs(chembl_con, limit_to_literature, df_sizes)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "activity ct-pairs", df_sizes)
    
    logging.info("add_cti_from_drug_mechanisms")
    df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set = get_drug_mechanism_ct_pairs.add_drug_mechanism_ct_pairs(
        df_combined, chembl_con)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "dm ct-pairs", df_sizes)

    logging.info("add_cti_annotations")
    df_combined = add_dti_annotations.add_dti_annotations(
        df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "DTI annotations", df_sizes)

    logging.info("add_all_chembl_compound_properties")
    df_combined, df_cpd_props, atc_levels = add_chembl_compound_properties.add_all_chembl_compound_properties(
        df_combined, chembl_con, limit_to_literature)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "ChEMBL props", df_sizes)

    logging.info("remove_compounds_without_smiles_and_mixtures")
    df_combined = clean_dataset.remove_compounds_without_smiles_and_mixtures(df_combined, chembl_con)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "removed smiles", df_sizes)

    logging.info("add_chembl_target_class_annotations")
    df_combined, target_classes_level1, target_classes_level2 = add_chembl_target_class_annotations.add_chembl_target_class_annotations(
        df_combined, chembl_con, output_path, write_to_csv, write_to_excel, delimiter)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "tclass annotations", df_sizes)

    logging.info("add_rdkit_compound_descriptors")
    if calculate_RDKit:
        df_combined = add_rdkit_compound_descriptors.add_rdkit_compound_descriptors(
            df_combined)
        if logging.DEBUG >= logging.root.level:
            get_stats.add_dataset_sizes(df_combined, "RDKit props", df_sizes)

    logging.info("clean_dataset")
    df_combined = clean_dataset.clean_dataset(df_combined, calculate_RDKit)
    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_combined, "clean df", df_sizes)

    logging.info("sanity_checks")
    sanity_checks.sanity_checks(df_combined, df_cpd_props, atc_levels,
                                target_classes_level1, target_classes_level2, calculate_RDKit)

    if limit_to_literature: 
        limited_flag = "literature_only"
    else:
        limited_flag = "all_sources"
    logging.info("write_BF_to_file")
    min_nof_cpds_BF = 100
    df_combined_annotated = write_subsets.write_BF_to_file(df_combined, 
                                                           chembl_version, min_nof_cpds_BF,
                                                           output_path, write_BF, write_to_csv, write_to_excel, delimiter,
                                                           limited_flag, calculate_RDKit, 
                                                           df_sizes)

    logging.info("write_B_to_file")
    min_nof_cpds_B = 100
    df_combined_annotated = write_subsets.write_B_to_file(df_combined, df_combined_annotated,
                                                          chembl_version, min_nof_cpds_B,
                                                          output_path, write_B, write_to_csv, write_to_excel, delimiter,
                                                          limited_flag, calculate_RDKit, 
                                                          df_sizes)

    logging.info("write_full_dataset_to_file")
    write_subsets.write_full_dataset_to_file(df_combined_annotated, 
                                             chembl_version,
                                             output_path, write_full_dataset, write_to_csv, write_to_excel, delimiter,
                                             limited_flag, calculate_RDKit)

    logging.info("output_stats")
    get_stats.output_stats(df_combined_annotated, output_path, write_to_csv, write_to_excel, delimiter)

    if logging.DEBUG >= logging.root.level:
        get_stats.output_debug_sizes(df_sizes, output_path, write_to_csv, write_to_excel, delimiter)


