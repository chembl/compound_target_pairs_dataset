# import test_utils.add_dataset_sizes
import get_activity_ct_pairs
import get_drug_mechanism_ct_pairs
import add_dti_annotations
import add_chembl_compound_properties
import add_chembl_target_class_annotations
import add_rdkit_compound_descriptors
import clean_dataset
import sanity_checks
import write_subsets
import get_stats


def print_df_combined_stats(df_combined):
    print("df len:", len(df_combined))
    print("#cols: ", len(df_combined.columns))


########### ###########
def get_dataset(chembl_con,
                chembl_version,
                output_path,
                limit_to_literature,
                calculate_RDKit,
                write_to_csv, delimiter,
                write_to_excel,
                write_full_dataset, write_BF, write_B):
    # # TODO: include?
    # # list with sizes of full dataset
    # all_lengths = []
    # # list with sizes of dataset with pchembl values
    # # these statistics are purely based on removing compound-target pairs without pchembl information
    # # i.e., the subset of the dataset is determined by the given data parameter and not recalculated (see below)
    # all_lengths_pchembl = []

    print("get_aggregated_acticity_ct_pairs")
    df_combined = get_activity_ct_pairs.get_aggregated_acticity_ct_pairs(chembl_con, limit_to_literature)
    print_df_combined_stats(df_combined)

    print("add_cti_from_drug_mechanisms")
    df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set = get_drug_mechanism_ct_pairs.add_drug_mechanism_ct_pairs(
        df_combined, chembl_con)
    print_df_combined_stats(df_combined)

    print("add_cti_annotations")
    df_combined = add_dti_annotations.add_dti_annotations(
        df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set)
    print_df_combined_stats(df_combined)

    print("add_all_chembl_compound_properties")
    df_combined, df_cpd_props, atc_levels = add_chembl_compound_properties.add_all_chembl_compound_properties(
        df_combined, chembl_con, limit_to_literature)
    print_df_combined_stats(df_combined)

    print("remove_irrelevant_compounds")
    df_combined = clean_dataset.remove_compounds_without_smiles_and_mixtures(df_combined, chembl_con)
    print_df_combined_stats(df_combined)

    print("add_chembl_target_class_annotations")
    df_combined, target_classes_level1, target_classes_level2 = add_chembl_target_class_annotations.add_chembl_target_class_annotations(
        df_combined, chembl_con, output_path, write_to_csv, write_to_excel, delimiter)
    print_df_combined_stats(df_combined)

    print("add_rdkit_compound_descriptors")
    if calculate_RDKit:
        df_combined = add_rdkit_compound_descriptors.add_rdkit_compound_descriptors(
            df_combined)
    print_df_combined_stats(df_combined)

    print("clean_dataset")
    df_combined = clean_dataset.clean_dataset(df_combined, calculate_RDKit)
    print_df_combined_stats(df_combined)

    print("sanity_checks")
    sanity_checks.sanity_checks(df_combined, df_cpd_props, atc_levels,
                                target_classes_level1, target_classes_level2, calculate_RDKit)

    print("write_BF")
    min_nof_cpds_BF = 100
    df_combined_annotated = write_subsets.write_BF_to_file(df_combined, 
                                                           chembl_version, min_nof_cpds_BF,
                                                           output_path, write_BF, write_to_csv, write_to_excel, delimiter,
                                                           calculate_RDKit)
    print("full df:")
    print_df_combined_stats(df_combined)
    print("annotated:")
    print_df_combined_stats(df_combined_annotated)

    print("write_B")
    min_nof_cpds_B = 100
    df_combined_annotated = write_subsets.write_B_to_file(df_combined, df_combined_annotated,
                                                          chembl_version, min_nof_cpds_B,
                                                          output_path, write_B, write_to_csv, write_to_excel, delimiter,
                                                          calculate_RDKit)
    print("full df:")
    print_df_combined_stats(df_combined)
    print("annotated:")
    print_df_combined_stats(df_combined_annotated)

    print("write_full_dataset_to_file")
    write_subsets.write_full_dataset_to_file(df_combined_annotated, 
                                             chembl_version,
                                             output_path, write_full_dataset, write_to_csv, write_to_excel, delimiter,
                                             calculate_RDKit)
    print("full df:")
    print_df_combined_stats(df_combined)
    print("annotated:")
    print_df_combined_stats(df_combined_annotated)

    print("print_stats")
    get_stats.print_stats(df_combined_annotated)
