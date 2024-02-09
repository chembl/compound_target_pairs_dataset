import logging
import pandas as pd
import sqlite3


def remove_compounds_without_smiles_and_mixtures(
    df_combined: pd.DataFrame, chembl_con: sqlite3.Connection
) -> pd.DataFrame:
    """
    Remove

    - compounds without a smiles
    - compounds with smiles containing a dot (mixtures and salts).

    Since compound information is aggregated for the parents of salts,
    the number of smiles with a dot is relatively low.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: Pandas DataFrame with compound-target pairs with a smiles that does not contain a '.'
    :rtype: pd.DataFrame
    """
    # Double-check that rows with a SMILES containing a '.' are the parent structures,
    # i.e., there was no error in using salt information instead of parent information.
    sql = """
    SELECT DISTINCT mh.molregno as salt_molregno, mh.parent_molregno
    FROM molecule_hierarchy mh
    """
    df_hierarchy = pd.read_sql_query(sql, con=chembl_con)

    smiles_with_dot = df_combined[
        df_combined["canonical_smiles"].notnull()
        & df_combined["canonical_smiles"].str.contains(".", regex=False)
    ][["canonical_smiles", "parent_molregno"]].drop_duplicates()

    for parent_molregno in set(smiles_with_dot["parent_molregno"]):
        assert (
            len(df_hierarchy[df_hierarchy["parent_molregno"] == parent_molregno]) > 0
        ), f"The compound with the parent_moregno {parent_molregno} does not occur \
                as a parent_molregno in the molecule_hierarchy."
        df_salt_molregno = df_hierarchy[
            df_hierarchy["salt_molregno"] == parent_molregno
        ]
        assert df_salt_molregno["salt_molregno"].equals(
            df_salt_molregno["parent_molregno"]
        ), f"The compound with the parent_moregno {parent_molregno} occurs as a salt \
                but has a different parent in the molecule_hierarchy."

    # Double-check that the SMILES is indeed the SMILES for the parent structure.
    sql = """
    SELECT DISTINCT mh.parent_molregno, struct.canonical_smiles
    FROM molecule_hierarchy mh
    INNER JOIN compound_structures struct
        ON mh.parent_molregno = struct.molregno
    """
    df_parent_smiles = pd.read_sql_query(sql, con=chembl_con)

    for parent_molregno in set(smiles_with_dot["parent_molregno"]):
        parent_smiles_in_chembl = df_parent_smiles[
            df_parent_smiles["parent_molregno"] == parent_molregno
        ]["canonical_smiles"].item()
        parent_smiles_in_df = smiles_with_dot[
            smiles_with_dot["parent_molregno"] == parent_molregno
        ]["canonical_smiles"].item()
        assert (
            parent_smiles_in_chembl == parent_smiles_in_df
        ), f"The smiles for the compound {parent_molregno} ({parent_smiles_in_df}) \
                in the dataframe is not the same as \
                the smiles for the compound in ChEMBL ({parent_smiles_in_chembl})."

    # Remove rows that contain a SMILES with a dot or that don't have a SMILES.
    len_missing_smiles = len(df_combined[df_combined["canonical_smiles"].isnull()])
    len_smiles_w_dot = len(
        df_combined[
            df_combined["parent_molregno"].isin(set(smiles_with_dot["parent_molregno"]))
        ]
    )
    logging.debug(f"#Compounds without a SMILES: {len_missing_smiles}")
    logging.debug(f"#SMILES with a dot: {len_smiles_w_dot}")

    df_combined = df_combined[
        (df_combined["canonical_smiles"].notnull())
        & ~(
            df_combined["parent_molregno"].isin(set(smiles_with_dot["parent_molregno"]))
        )
    ]

    return df_combined


def clean_none_values(df_combined):
    """
    Change nan values and empty strings to None for consistency.
    """
    # Change all None / nan values to None
    df_combined = df_combined.where(pd.notnull(df_combined), None)
    # replace empty strings with None
    df_combined = df_combined.replace("", None).reset_index(drop=True)

    return df_combined


def set_types_to_int(df_combined, calculate_RDKit):
    """
    Set the type of relevant columns to Int64.
    """
    df_combined = df_combined.astype(
        {
            "first_approval": "Int64",
            "usan_year": "Int64",
            "first_publication_cpd_target_pair_BF": "Int64",
            "first_publication_cpd_target_pair_w_pchembl_BF": "Int64",
            "first_publication_cpd_target_pair_B": "Int64",
            "first_publication_cpd_target_pair_w_pchembl_B": "Int64",
            "first_publication_cpd": "Int64",
            "hba": "Int64",
            "hbd": "Int64",
            "rtb": "Int64",
            "num_ro5_violations": "Int64",
            "aromatic_rings": "Int64",
            "heavy_atoms": "Int64",
            "hba_lipinski": "Int64",
            "hbd_lipinski": "Int64",
            "num_lipinski_ro5_violations": "Int64",
        }
    )

    if calculate_RDKit:
        df_combined = df_combined.astype(
            {
                "num_aliphatic_carbocycles": "Int64",
                "num_aliphatic_heterocycles": "Int64",
                "num_aliphatic_rings": "Int64",
                "num_aromatic_carbocycles": "Int64",
                "num_aromatic_heterocycles": "Int64",
                "num_aromatic_rings": "Int64",
                "num_heteroatoms": "Int64",
                "num_saturated_carbocycles": "Int64",
                "num_saturated_heterocycles": "Int64",
                "num_saturated_rings": "Int64",
                "ring_count": "Int64",
                "num_stereocentres": "Int64",
                "aromatic_atoms": "Int64",
                "aromatic_c": "Int64",
                "aromatic_n": "Int64",
                "aromatic_hetero": "Int64",
            }
        )

    return df_combined


def round_floats(df_combined, decimal_places=4):
    """
    Round float columns to <decimal_places> decimal places.
    This does not apply to max_phase.
    """
    for _, (col, dtype) in enumerate(df_combined.dtypes.to_dict().items()):
        if ((dtype == "float64") or (dtype == "Float64")) and col != "max_phase":
            df_combined[col] = df_combined[col].round(decimals=decimal_places)

    return df_combined


def reorder_columns(df_combined, calculate_RDKit):
    """
    Reorder the columns in the DataFrame.
    """
    len_columns_before = len(df_combined.columns)

    compound_target_pair_columns = [
        "parent_molregno",
        "parent_chemblid",
        "parent_pref_name",
        "max_phase",
        "first_approval",
        "usan_year",
        "black_box_warning",
        "prodrug",
        "oral",
        "parenteral",
        "topical",
        "tid",
        "mutation",
        "target_chembl_id",
        "target_pref_name",
        "target_type",
        "organism",
        "tid_mutation",
        "cpd_target_pair",
        "cpd_target_pair_mutation",
    ]
    aggregated_values = [
        "pchembl_value_mean_BF",
        "pchembl_value_max_BF",
        "pchembl_value_median_BF",
        "first_publication_cpd_target_pair_BF",
        "first_publication_cpd_target_pair_w_pchembl_BF",
        "pchembl_value_mean_B",
        "pchembl_value_max_B",
        "pchembl_value_median_B",
        "first_publication_cpd_target_pair_B",
        "first_publication_cpd_target_pair_w_pchembl_B",
    ]
    DTI_annotations = ["therapeutic_target", "DTI"]
    first_publication_cpd = ["first_publication_cpd"]
    chembl_compound_props = [
        "mw_freebase",
        "alogp",
        "hba",
        "hbd",
        "psa",
        "rtb",
        "ro3_pass",
        "num_ro5_violations",
        "cx_most_apka",
        "cx_most_bpka",
        "cx_logp",
        "cx_logd",
        "molecular_species",
        "full_mwt",
        "aromatic_rings",
        "heavy_atoms",
        "qed_weighted",
        "mw_monoisotopic",
        "full_molformula",
        "hba_lipinski",
        "hbd_lipinski",
        "num_lipinski_ro5_violations",
    ]
    chembl_structures = ["standard_inchi", "standard_inchi_key", "canonical_smiles"]
    ligand_efficieny_metrics = [
        "LE_B",
        "BEI_B",
        "SEI_B",
        "LLE_B",
        "LE_BF",
        "BEI_BF",
        "SEI_BF",
        "LLE_BF",
    ]
    chembl_target_annotations = ["atc_level1", "target_class_l1", "target_class_l2"]
    rdkit_columns = [
        "fraction_csp3",
        "ring_count",
        "num_aliphatic_rings",
        "num_aliphatic_carbocycles",
        "num_aliphatic_heterocycles",
        "num_aromatic_rings",
        "num_aromatic_carbocycles",
        "num_aromatic_heterocycles",
        "num_saturated_rings",
        "num_saturated_carbocycles",
        "num_saturated_heterocycles",
        "num_stereocentres",
        "num_heteroatoms",
        "aromatic_atoms",
        "aromatic_c",
        "aromatic_n",
        "aromatic_hetero",
        "scaffold_w_stereo",
        "scaffold_wo_stereo",
    ]
    filtering_columns = [
        "pair_mutation_in_dm_table",
        "pair_in_dm_table",
        "keep_for_binding",
    ]

    if calculate_RDKit:
        columns = (
            compound_target_pair_columns
            + aggregated_values
            + DTI_annotations
            + first_publication_cpd
            + chembl_compound_props
            + chembl_structures
            + ligand_efficieny_metrics
            + chembl_target_annotations
            + rdkit_columns
            + filtering_columns
        )
        df_combined = df_combined[columns]
    else:
        columns = (
            compound_target_pair_columns
            + aggregated_values
            + DTI_annotations
            + first_publication_cpd
            + chembl_compound_props
            + chembl_structures
            + ligand_efficieny_metrics
            + chembl_target_annotations
            + filtering_columns
        )
        df_combined = df_combined[columns]

    len_columns_after = len(df_combined.columns)
    assert (
        len_columns_before == len_columns_after
    ), f"Different number of columns after reordering (before: {len_columns_before}, after: {len_columns_after})."

    return df_combined


def clean_dataset(df_combined: pd.DataFrame, calculate_RDKit: bool) -> pd.DataFrame:
    """
    Clean the dataset by

    - changing nan values and empty strings to None
    - setting the type of relevant columns to Int64
    - rounding floats to 4 decimal places (with the exception of max_phase which is not rounded)
    - reordering columns
    - sorting rows by cpd_target_pair_mutation

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param calculate_RDKit: True if the DataFrame contains RDKit-based compound properties
    :type calculate_RDKit: bool
    :return: Cleaned pandas DataFrame with compound-target pairs
    :rtype: pd.DataFrame
    """
    df_combined = clean_none_values(df_combined)
    df_combined = set_types_to_int(df_combined, calculate_RDKit)
    df_combined = round_floats(df_combined, decimal_places=4)
    df_combined = reorder_columns(df_combined, calculate_RDKit)
    df_combined = df_combined.sort_values(by=["cpd_target_pair_mutation"]).reset_index(
        drop=True
    )
    return df_combined
