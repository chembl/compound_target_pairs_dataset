import logging
import numpy as np
import pandas as pd
import sqlite3

import get_stats

########### Get Initial Compound-Target Data From ChEMBL ###########
def get_compound_target_pairs_with_pchembl(chembl_con: sqlite3.Connection, limit_to_literature: bool, df_sizes: list[list[int], list[int]]) -> pd.DataFrame:
    """
    Query ChEMBL activities and related assay for compound-target pairs with an associated pchembl value.  
    Compound-target pairs are required to have a pchembl value.  
    Salt forms of compounds are mapped to their parent form.  
    If limit_to_literature is true, only literature sources will be considered. Otherwise, all sources are included.  
    Includes information about targets, mutations and year of publication (based on docs).  

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Include only literature sources if True. Include all available sources otherwise.
    :type limit_to_literature: bool
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :return: Pandas DataFrame with compound-target pairs with a pchembl value.
    :rtype: pd.DataFrame
    """
    # NOTE: DO NOT USE DISTINCT
    # This query does not capture act.activity_id. 
    # There can be mulitple activities with different activity_ids but the same queries values. 
    # For accurate mean, median, and max pchembl these additional rows are important.
    sql = '''
    SELECT act.pchembl_value, 
        md.molregno as parent_molregno, md.chembl_id as parent_chemblid, md.pref_name as parent_pref_name,
        md.max_phase, md.first_approval, md.usan_year, md.black_box_warning, 
        md.prodrug, md.oral, md.parenteral, md.topical, 
        ass.assay_type, ass.tid, 
        vs.mutation,
        td.chembl_id as target_chembl_id, td.pref_name as target_pref_name, td.target_type, td.organism, 
        docs.year
    FROM activities act
    INNER JOIN molecule_hierarchy mh 
        ON act.molregno = mh.molregno         -- act.molregno = salt_molregno
    INNER JOIN molecule_dictionary md
        ON mh.parent_molregno = md.molregno   -- compound information based on parent compound
    INNER JOIN assays ass 
        ON  act.assay_id = ass.assay_id
    LEFT JOIN variant_sequences vs
        ON ass.variant_id = vs.variant_id
    INNER JOIN target_dictionary td
        ON ass.tid = td.tid
    LEFT JOIN docs
        ON act.doc_id = docs.doc_id
    WHERE act.pchembl_value is not null
        and act.potential_duplicate = 0
        and act.standard_relation = '='
        and act.data_validity_comment is null
        and td.tid <>22226                    -- exclude unchecked targets
        and td.target_type like '%PROTEIN%'
    '''
    if limit_to_literature:
        sql += '''    and docs.src_id = 1'''

    df_mols = pd.read_sql_query(sql, con=chembl_con)

    # Set relevant combinations of columns for easier processing later
    # target_id_mutation
    df_mols['tid_mutation'] = np.where(df_mols['mutation'].notnull(),
                                       df_mols['tid'].astype(
                                           'str')+'_'+df_mols['mutation'],
                                       df_mols['tid'].astype('str'))
    # compound-target association
    df_mols['cpd_target_pair'] = df_mols.agg('{0[parent_molregno]}_{0[tid]}'.format, axis=1)
    df_mols['cpd_target_pair_mutation'] = df_mols.agg('{0[parent_molregno]}_{0[tid_mutation]}'.format, axis=1)

    if logging.DEBUG >= logging.root.level:
        get_stats.add_dataset_sizes(df_mols,  "initial query", df_sizes)

    return df_mols



########### Calculate Mean, Median, and Max *pchembl* Values for Each Compound-Target Pair ###########
def get_average_info(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
    """
    Aggregate the information about compound-target pairs for which there is more than one entry into one entry. 
    Compound-target pairs are considered equal if parent_molregno (internal compound ID) and tid_mutation (target ID + mutation annotations) are equal.

    The following values are aggregated:  

    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+
    | pchembl_value_mean                            | mean pchembl value for a compound-target pair                                                 |
    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+
    | pchembl_value_max                             | maximum pchembl value for a compound-target pair                                              |
    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+
    | pchembl_value_median                          | median pchembl value for a compound-target pair                                               |
    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+
    | first_publication_cpd_target_pair             | first publication in ChEMBL with this compound-target pair                                    |
    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+
    | first_publication_cpd_target_pair_w_pchembl   | first publication in ChEMBL with this compound-target pair and an associated pchembl value    |
    +-----------------------------------------------+-----------------------------------------------------------------------------------------------+

    :param df: Pandas DataFrame with compound-target pairs for which the information should be aggregated.
    :type df: pd.DataFrame
    :param suffix: Suffix indicating the type of the given DataFrame, e.g., _B for binding assays, _BF for binding+functional assays.
    :type suffix: str
    :return: Pandas DataFrame with 'parent_molregno', 'tid_mutation', and the aggregated columns.
    :rtype: pd.DataFrame
    """
    # pchembl mean, max, median
    df[f"pchembl_value_mean_{suffix}"] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('mean')
    df[f"pchembl_value_max_{suffix}"] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('max')
    df[f"pchembl_value_median_{suffix}"] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('median')

    # first publication of pair
    df[f"first_publication_cpd_target_pair_{suffix}"] = df.groupby(['parent_molregno', 'tid_mutation'])['year'].transform('min')

    # first publication of pair with pchembl value
    df_mols_all_first_publication_pchembl = df[df['pchembl_value'].notnull()] \
        .groupby(['parent_molregno', 'tid_mutation'])['year'].min().reset_index() \
        .rename(columns={'year': f"first_publication_cpd_target_pair_w_pchembl_{suffix}"})
    df = df.merge(df_mols_all_first_publication_pchembl, on=['parent_molregno', 'tid_mutation'], how='left')

    # return relevant summarised information without duplicates
    df = df[['parent_molregno', 'tid_mutation',
            f"pchembl_value_mean_{suffix}", f"pchembl_value_max_{suffix}", f"pchembl_value_median_{suffix}",
            f"first_publication_cpd_target_pair_{suffix}", f"first_publication_cpd_target_pair_w_pchembl_{suffix}"]].drop_duplicates()
    
    return df



########### Get Aggregated Compound-Target Pair Information ###########
def get_aggregated_activity_ct_pairs(chembl_con: sqlite3.Connection, limit_to_literature: bool, df_sizes: list[list[int], list[int]]) -> pd.DataFrame:
    """
    Get dataset of compound target-pairs with an associated pchembl value 
    with pchembl and publication dates aggregated into one entry per pair.

    Values are aggregated for 

    - a subset of the initial dataset based on binding and functional assays (suffix '_BF') and 
    - a subset of the initial dataset set on only binding assays (suffix '_B'). 

    Therefore, there are two columns for pchembl_value_mean, _max, _median, first_publication_cpd_target_pair and first_publication_cpd_target_pair_w_pchembl, 
    one with the suffix '_BF' based on binding + functional data and one with the suffix '_B' based on only binding data.

    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Include only literature sources if True. Include all available sources otherwise.
    :type limit_to_literature: bool
    :param df_sizes: List of intermediate sized of the dataset used for debugging.
    :type df_sizes: list[list[int], list[int]]
    :return: Pandas Dataframe with compound-target pairs based on ChEMBL activity data aggregated into one entry per compound-target pair.
    :rtype: pd.DataFrame
    """
    df_mols = get_compound_target_pairs_with_pchembl(chembl_con, limit_to_literature, df_sizes)

    # Summarise the information for binding and functional assays
    suffix = 'BF'
    df_mols_BF = df_mols[(df_mols['assay_type'] == 'B') | (df_mols['assay_type'] == 'F')].copy()
    df_mols_BF = get_average_info(df_mols_BF, suffix)

    # Summarise the information for only binding assays
    suffix = 'B'
    df_mols_B = df_mols[df_mols['assay_type'] == 'B'].copy()
    df_mols_B = get_average_info(df_mols_B, suffix)

    # Combine both into one table with two columns per value
    # (one with suffix '_BF' for binding+functional and one with suffix '_B' for binding).
    # df_mols_B is a subset of the compound-target pairs of df_mols_BF
    df_combined = df_mols_BF.merge(df_mols_B,
                                   on=['parent_molregno', 'tid_mutation'], how='left')
    # Merge with other information from df_mols
    # left merge because df_mols may contain assays that are of other types than binding / functional
    df_combined = df_combined.merge(df_mols.drop(columns=['pchembl_value', 'year', 'assay_type']).drop_duplicates(),
                                    on=['parent_molregno', 'tid_mutation'], how='left')

    return df_combined
