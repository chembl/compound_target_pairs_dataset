import logging
import sqlite3

import pandas as pd

from dataset import Dataset


########### Extract Drug-Target Interactions From the drug_mechanism Table ###########
def get_drug_mechanisms_interactions(chembl_con: sqlite3.Connection) -> pd.DataFrame:
    """
    Extract the known compound-target interactions from the ChEMBL drug_mechanisms table.
    Note: While the interactions are mostly between drugs and targets,
    the table also includes some known interactions between
    compounds with a max_phase < 4 and their targets.

    Only entries with a disease_efficacy of 1 are taken into account,
    i.e., the target is believed to play a role in the efficacy of the drug.

    *disease_efficacy: Flag to show whether the target assigned is believed
    to play a role in the efficacy of the drug in the indication(s)
    for which it is approved (1 = yes, 0 = no).*

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: Pandas DataFrame with compound-target pairs
        from the drug_mechanism table with disease relevance.
    :rtype: pd.DataFrame
    """
    sql = """
    SELECT DISTINCT mh.parent_molregno, dm.tid
    FROM drug_mechanism dm
    INNER JOIN molecule_hierarchy mh
        ON dm.molregno = mh.molregno
    INNER JOIN molecule_dictionary md
        ON mh.parent_molregno = md.molregno
    WHERE dm.disease_efficacy = 1
        and dm.tid is not null
    """

    df_dti = pd.read_sql_query(sql, con=chembl_con)

    return df_dti


def get_relevant_tid_mappings(chembl_con: sqlite3.Connection) -> pd.DataFrame:
    """
    Get DataFrame with mappings from target id to their related target ids
    based on the target_relations table.
    The following mappings are considered:

    +-------------------------------+-----------------------+----------------+
    |protein family                 | -[superset of]->      | single protein |
    +-------------------------------+-----------------------+----------------+
    |protein complex                | -[superset of]->      | single protein |
    +-------------------------------+-----------------------+----------------+
    |protein complex group          | -[superset of]->      | single protein |
    +-------------------------------+-----------------------+----------------+
    |single protein                 | -[equivalent to]->    | single protein |
    +-------------------------------+-----------------------+----------------+
    |chimeric protein               | -[superset of]->      | single protein |
    +-------------------------------+-----------------------+----------------+
    |protein-protein interaction    | -[superset of]->      | single protein |
    +-------------------------------+-----------------------+----------------+

    These mappings can be used to increase the number of target ids
    for which there is data in the drug_mechanisms table.
    For example, for *protein family -[superset of]-> single protein* this means:
    If there is a known relevant interaction between a compound and a protein family,
    interactions between the compound and single proteins of that protein family
    are considered to be known interactions as well.

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: Pandas DataFrame with mappings from tid to related tid
        for the defined subset of target relations.
    :rtype: pd.DataFrame
    """
    sql = """
    SELECT DISTINCT tr.tid, tr.relationship, tr.related_tid, 
        td1.pref_name as pref_name_1, td1.target_type as target_type_1, td1.organism as organism_1, 
        td2.pref_name as pref_name_2, td2.target_type as target_type_2, td2.organism as organism_2 
    FROM target_relations tr
    INNER JOIN target_dictionary td1
        ON tr.tid = td1.tid
    INNER JOIN target_dictionary td2
        ON tr.related_tid = td2.tid
    """
    df_related_targets = pd.read_sql_query(sql, con=chembl_con)

    protein_family_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "PROTEIN FAMILY")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "SUPERSET OF")
    ]

    protein_complex_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "PROTEIN COMPLEX")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "SUPERSET OF")
    ]

    protein_complex_group_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "PROTEIN COMPLEX GROUP")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "SUPERSET OF")
    ]

    single_protein_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "SINGLE PROTEIN")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "EQUIVALENT TO")
    ]

    chimeric_protein_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "CHIMERIC PROTEIN")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "SUPERSET OF")
    ]

    ppi_mapping = df_related_targets[
        (df_related_targets["target_type_1"] == "PROTEIN-PROTEIN INTERACTION")
        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
        & (df_related_targets["relationship"] == "SUPERSET OF")
    ]

    relevant_tid_mappings = pd.concat(
        [
            protein_family_mapping,
            protein_complex_mapping,
            protein_complex_group_mapping,
            single_protein_mapping,
            chimeric_protein_mapping,
            ppi_mapping,
        ]
    )

    return relevant_tid_mappings


def add_annotations_to_drug_mechanisms_cti(
    chembl_con: sqlite3.Connection, cpd_target_pairs: pd.DataFrame
) -> pd.DataFrame:
    """
    Add additional information to the compound-target pairs from the drug_mechanisms table
    to match the information that is present in the compound-target pairs table based on activities.

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param cpd_target_pairs: Pandas DataFrame with compound-target pairs
        from the drug_mechanism table.
    :type cpd_target_pairs: pd.DataFrame
    :return: Updated pandas DataFrame with the additional annotations.
    :rtype: pd.DataFrame
    """
    ##### Set columns existing in the df_results table. #####
    # None of the targets from the drug mechanism table have any mutation annotation,
    # hence tid_mutation = tid
    cpd_target_pairs["tid_mutation"] = cpd_target_pairs["tid"].astype("str")

    cpd_target_pairs["cpd_target_pair"] = [
        f"{a}_{b}"
        for a, b in zip(cpd_target_pairs["parent_molregno"], cpd_target_pairs["tid"])
    ]
    cpd_target_pairs["cpd_target_pair_mutation"] = [
        f"{a}_{b}"
        for a, b in zip(
            cpd_target_pairs["parent_molregno"], cpd_target_pairs["tid_mutation"]
        )
    ]

    # New column: is the compound target pair in the drug_mechanism table?
    cpd_target_pairs["pair_mutation_in_dm_table"] = True
    cpd_target_pairs["pair_in_dm_table"] = True

    ##### Query and combine compound information with compound-target pairs #####
    sql = """
    SELECT md.molregno as parent_molregno, 
        md.chembl_id as parent_chemblid, md.pref_name as parent_pref_name,
        md.max_phase, md.first_approval, md.usan_year, md.black_box_warning, 
        md.prodrug, md.oral, md.parenteral, md.topical
    FROM molecule_dictionary md
    """

    df_compound_info = pd.read_sql_query(sql, con=chembl_con)
    cpd_target_pairs = cpd_target_pairs.merge(
        df_compound_info, on="parent_molregno", how="left"
    )

    ##### Query and combine target information with compound-target pairs #####
    sql = """
    SELECT td.tid, td.chembl_id as target_chembl_id, td.pref_name as target_pref_name, td.target_type, td.organism
    FROM target_dictionary td
    """

    df_target_info = pd.read_sql_query(sql, con=chembl_con)
    # Fix problems with null not being recognised as None
    df_target_info.loc[df_target_info["organism"].astype(str) == "null", "organism"] = (
        None
    )
    cpd_target_pairs = cpd_target_pairs.merge(df_target_info, on="tid", how="left")

    return cpd_target_pairs


def get_drug_mechanism_ct_pairs(chembl_con: sqlite3.Connection) -> pd.DataFrame:
    """
    Get compound-target pairs from the drug_mechanism table
    with all the columns that are present in the compound-target pairs based on activities.
    Relevant mappings of target ids to related target ids are taken into account.

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: Pandas DataFrame with compound-target interactions from the drug_mechanism table.
    :rtype: pd.DataFrame
    """
    # get known compound-target interactions (CTI) from the drug_mechanisms table
    df_dti = get_drug_mechanisms_interactions(chembl_con)

    # Query target_relations for related target ids
    # to increase the number of target ids for which there is data in the drug_mechanisms table.
    relevant_tid_mappings = get_relevant_tid_mappings(chembl_con)
    # table with mapped target ids
    df_dti_mapped_targets = df_dti.merge(relevant_tid_mappings, on="tid", how="inner")

    # combine CTIs from drug_mechanism table with mapped CTIs
    cpd_target_pairs = pd.concat(
        [
            df_dti[["parent_molregno", "tid"]],
            df_dti_mapped_targets[["parent_molregno", "related_tid"]].rename(
                columns={"related_tid": "tid"}
            ),
        ]
    ).drop_duplicates()

    cpd_target_pairs = add_annotations_to_drug_mechanisms_cti(
        chembl_con, cpd_target_pairs
    )

    return cpd_target_pairs


########### Add Compounds From the drug_mechanism Table to the Dataset ###########
def add_drug_mechanism_ct_pairs(dataset: Dataset, chembl_con: sqlite3.Connection):
    """
    Add compound-target pairs from the drug_mechanism table
    that are not in the dataset based on the initial ChEMBL query.
    These are compound-target pairs for which there is no associated pchembl value data.
    Since the pairs are known interactions,
    they are added to the dataset despite not having a pchembl value.
    Add the set of compound-target pairs in the drug_mechanism table and
    the set of targets in the drug_mechanism table to the dataset.

    :param dataset: Pandas Dataframe with compound-target pairs based on ChEMBL activity data
    :type dataset: Dataset
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    """
    cpd_target_pairs = get_drug_mechanism_ct_pairs(chembl_con)
    dataset.drug_mechanism_pairs_set = set(
        f"{a}_{b}"
        for a, b in zip(cpd_target_pairs["parent_molregno"], cpd_target_pairs["tid"])
    )

    dataset.drug_mechanism_targets_set = set(cpd_target_pairs["tid"])

    # Add a new column *pair_mutation_in_dm_table* which is set to True if the compound target pair
    # (taking mutation annotations into account) is in the drug_mechanism table.
    dataset.df_result["pair_mutation_in_dm_table"] = False
    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair_mutation"].isin(
                dataset.drug_mechanism_pairs_set
            )
        ),
        "pair_mutation_in_dm_table",
    ] = True

    # Add a new column *pair_in_dm_table* which is set to True if the compound target pair
    # (NOT taking mutation annotations into account) is in the drug_mechanism table.
    dataset.df_result["pair_in_dm_table"] = False
    dataset.df_result.loc[
        (dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)),
        "pair_in_dm_table",
    ] = True

    ##### Limit the drug_mechanism pairs to the ones that are not yet in the dataset. #####
    # Mutation annotations are taken into account.
    # Therefore, *(cpd A, target B without mutation)* will be added
    # if a pchembl is present for *(cpd A, target B with mutation C)*
    # but not for *(cpd A, target B without mutation)*.
    cpd_target_pairs = cpd_target_pairs[
        ~(
            cpd_target_pairs["cpd_target_pair_mutation"].isin(
                set(dataset.df_result["cpd_target_pair_mutation"])
            )
        )
    ].copy()

    logging.debug(
        "#Pairs not yet present based on binding or functional assays: %s",
        len(cpd_target_pairs),
    )

    # Combined data of existing query with new compound-target pairs.
    dataset.df_result = pd.concat([dataset.df_result, cpd_target_pairs])

    # Add a new column *keep_for_binding* which is set to True if the row should be kept
    # if you want to limit the dataset to only data based on binding assays.
    # Rows are kept if
    # - there is a binding data-based pchembl value or
    # - the compound-target pair (including mutation info) is in the drug_mechanism table
    dataset.df_result["keep_for_binding"] = False
    dataset.df_result.loc[
        (
            (dataset.df_result["pchembl_value_mean_B"].notnull())
            | (dataset.df_result["pair_mutation_in_dm_table"] == True)
        ),
        "keep_for_binding",
    ] = True
