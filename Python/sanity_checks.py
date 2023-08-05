import pandas as pd

########### Sanity checks for the dataset ###########
def check_null_values(df_combined: pd.DataFrame):
    """
    Check if any columns contain nan or null which aren't recognised as null values. 
    """
    for col in df_combined.columns:
        col_as_str = set(df_combined[df_combined[col].notnull()][col].astype(str))
        assert ('nan' not in col_as_str), f"Problem with unrecognised nan value in column {col}"
        assert ('null' not in col_as_str), f"Problem with unrecognised null value in column {col}"


def check_for_mixed_types(df_combined: pd.DataFrame):
    """
    Check that there are no mixed types in columns with dtype=object.
    """
    for col, dtype in df_combined.dtypes.to_dict().items():
        if dtype == object:
            col_original = set(df_combined[df_combined[col].notnull()][col])
            col_as_str = set(df_combined[df_combined[col].notnull()][col].astype(str))
            # is there a difference in the two sets (ignoring null values)
            assert (len(col_original-col_as_str) == 0), f"Mixed types in colum {col}: {col_original-col_as_str}"
            assert (len(col_as_str-col_original) == 0), f"Mixed types in colum {col}: {col_as_str-col_original}"


def check_pairs_without_pchembl_are_in_drug_mechanisms(df_combined: pd.DataFrame):
    """
    Check that rows without a pchembl value based on binding+functional assays (pchembl_x_BF) are in the drug_mechanism table.  
    Note that this is not true for the pchembl_x_B columns which are based on binding data only. 
    They may be in the table because there is data based on functional assays but no data based on binding assays. 
    All pchembl_value_x_BF columns without a pchembl should be in the dm table.
    """
    for pchembl_col in ['pchembl_value_mean_BF', 'pchembl_value_max_BF', 'pchembl_value_median_BF']:
        assert (df_combined[(df_combined[pchembl_col].isnull())].equals(
            df_combined[(df_combined['in_dm_table'] == True) & (df_combined[pchembl_col].isnull())])), \
            f"Missing pchembl value in column {pchembl_col}"


def check_ligand_efficiency_metrics(df_combined: pd.DataFrame):
    """
    Check that ligand efficiency metrics are only null 
    when at least one of the values used to calculate them is null.
    Ligand efficiency metrics are only null when at least one of the values used to calculate them is null.
    """
    for suffix in ['BF', 'B']:
        assert (df_combined[(df_combined['LE_'+suffix].isnull())].equals(
            df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull())
                        | (df_combined['heavy_atoms'].isnull())
                        | (df_combined['heavy_atoms'] == 0)])), f"Missing LE value in LE_{suffix}"

        assert (df_combined[(df_combined['BEI_'+suffix].isnull())].equals(
            df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull())
                        | (df_combined['mw_freebase'].isnull())
                        | (df_combined['mw_freebase'] == 0)])), f"Missing BEI value in BEI_{suffix}"

        assert (df_combined[(df_combined['SEI_'+suffix].isnull())].equals(
            df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull())
                        | (df_combined['psa'].isnull())
                        | (df_combined['psa'] == 0)])), f"Missing SEI value in SEI_{suffix}"

        assert (df_combined[(df_combined['LLE_'+suffix].isnull())].equals(
            df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull())
                        | (df_combined['alogp'].isnull())])), f"Missing LLE value in LLE_{suffix}"


def check_compound_props(df_combined: pd.DataFrame, df_cpd_props: pd.DataFrame):
    """
    Check that compound props are only null if 

    - the property in the parent_molregno is not in df_cpd_props
    - or if the value in the compound props table is null.
    """
    # missing values because the parent_molregno is not in the compound props table
    no_cpd_prop_info = len(df_combined[~df_combined['parent_molregno'].isin(
        set(df_cpd_props['parent_molregno']))])

    for col in df_cpd_props.columns:
        if col != 'parent_molregno':
            # missing values because the compound props query returns null (exists but is null)
            missing_values = len(df_combined[df_combined['parent_molregno'].isin(
                set(df_cpd_props[df_cpd_props[col].isnull()]['parent_molregno']))])
            null_values = no_cpd_prop_info+missing_values
            assert (null_values == len(
                df_combined[df_combined[col].isnull()])), f"Too many null values in {col}"


def check_atc_and_target_classes(df_combined: pd.DataFrame,
                                 atc_levels: pd.DataFrame,
                                 target_classes_level1: pd.DataFrame,
                                 target_classes_level2: pd.DataFrame):
    """
    Check that atc_level1 and target class information is only null 
    if the parent_molregno / target id is not in the respective table.
    """
    assert (df_combined[(df_combined['atc_level1'].isnull())].equals(
        df_combined[~df_combined['parent_molregno'].isin(set(atc_levels['parent_molregno']))])), \
            "Null values in atc_level1 are not exclusively because the parent_molregno is not in the atc_classification table."

    assert (df_combined[(df_combined['target_class_l1'].isnull())].equals(
        df_combined[~df_combined['tid'].isin(set(target_classes_level1['tid']))])), \
            "Null values in target_class_l1 are not exclusively because the tid is not in the protein_classification table."

    assert (df_combined[(df_combined['target_class_l2'].isnull())].equals(
        df_combined[~df_combined['tid'].isin(set(target_classes_level2['tid']))])), \
            "Null values in target_class_l2 are not exclusively because the tid is not in the protein_classification table."


def check_rdkit_props(df_combined: pd.DataFrame):
    """
    Check that columns set by the RDKit are only null if there is no canonical SMILES for the molecule.
    Scaffolds are excluded from this test because they can be None if the molecule is acyclic.
    """
    for col in ['fraction_csp3', 'num_aliphatic_carbocycles', 'num_aliphatic_heterocycles', 'num_aliphatic_rings',
                'num_aromatic_carbocycles', 'num_aromatic_heterocycles', 'num_aromatic_rings',
                'num_heteroatoms', 'num_saturated_carbocycles', 'num_saturated_heterocycles',
                'num_saturated_rings', 'ring_count', 'num_stereocentres',
                'aromatic_atoms', 'aromatic_c', 'aromatic_n', 'aromatic_hetero']:
        assert (len(df_combined[df_combined[col].isnull()]) == len(df_combined[df_combined['canonical_smiles'].isnull()].copy())), \
            f"Missing value in {col} despite a smiles being available."


def sanity_checks(df_combined: pd.DataFrame,
                  df_cpd_props: pd.DataFrame, atc_levels: pd.DataFrame, target_classes_level1: pd.DataFrame, target_classes_level2: pd.DataFrame,
                  calculate_RDKit: bool):
    """
    Check basic assumptions about the finished dataset, specifically:

    - no columns contain nan or null values which aren't recognised as null values
    - there are no mixed types in columns with dtype=object
    - rows without a pchembl value based on binding+functional assays (pchembl_x_BF) are in the drug_mechanism table
    - ligand efficiency metrics are only null when at least one of the values used to calculate them is null
    - compound props are only null if the compound is not in df_cpd_props or the value in that table is null
    - atc_level1 and target class information is only null if the parent_molregno / target id is not in the respective table
    - columns set by the RDKit are only null if there is no canonical SMILES for the molecule (excluding scaffolds)

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param df_cpd_props: Pandas DataFrame with compound properties and structures for all compound ids in ChEMBL.
    :type df_cpd_props: pd.DataFrame
    :param atc_levels: Pandas DataFrame with ATC annotations in ChEMBL
    :type atc_levels: pd.DataFrame
    :param target_classes_level1: Pandas DataFrame with mapping from target id to level 1 target class
    :type target_classes_level1: pd.DataFrame
    :param target_classes_level2: Pandas DataFrame with mapping from target id to level 2 target class
    :type target_classes_level2: pd.DataFrame
    :param calculate_RDKit: True if the DataFrame contains RDKit-based compound properties
    :type calculate_RDKit: bool
    """
    check_null_values(df_combined)
    check_for_mixed_types(df_combined)
    check_pairs_without_pchembl_are_in_drug_mechanisms(df_combined)
    check_ligand_efficiency_metrics(df_combined)
    check_compound_props(df_combined, df_cpd_props)
    check_atc_and_target_classes(df_combined, atc_levels, target_classes_level1, target_classes_level2)
    if calculate_RDKit:
        check_rdkit_props(df_combined)

