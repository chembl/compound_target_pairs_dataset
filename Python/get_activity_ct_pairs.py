import numpy as np
import pandas as pd

########### Get Initial Compound-Target Data From ChEMBL ###########
def get_initial_dataset(chembl_con, limit_to_literature):
    """
    Initial query for activities + related assay, mutation, target and docs information. 
    Compound-target pairs are required to have a pchembl value.

    :param chembl_con: _description_
    :type chembl_con: _type_
    :param limit_to_literature: _description_
    :type limit_to_literature: _type_
    :return: _description_
    :rtype: _type_
    """
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
        and data_validity_comment is null
        and td.tid <>22226                    -- exclude unchecked targets
        and td.target_type like '%PROTEIN%'
    '''
    if limit_to_literature:
        sql += '''    and docs.src_id = 1'''

    df_mols = pd.read_sql_query(sql, con=chembl_con)
    # target_id_mutation
    df_mols['tid_mutation'] = np.where(df_mols['mutation'].notnull(), 
                                    df_mols['tid'].astype('str')+'_'+df_mols['mutation'], 
                                    df_mols['tid'].astype('str'))
    # compound-target association
    df_mols['cpd_target_pair'] = df_mols.agg('{0[parent_molregno]}_{0[tid]}'.format, axis=1)
    df_mols['cpd_target_pair_mutation'] = df_mols.agg('{0[parent_molregno]}_{0[tid_mutation]}'.format, axis=1)
    # TODO: include this?
    # test_utils.add_dataset_sizes(df_mols, "init", all_lengths, all_lengths_pchembl)
    return df_mols



########### Calculate Mean, Median, and Max *pchembl* Values for Each Compound-Target Pair ###########
def get_average_info(df, suffix):
    # pchembl mean, max, median
    df['pchembl_value_mean_'+suffix] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('mean')
    df['pchembl_value_max_'+suffix] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('max')
    df['pchembl_value_median_'+suffix] = df.groupby(['parent_molregno', 'tid_mutation'])['pchembl_value'].transform('median')
    
    # first publication of pair
    df['first_publication_cpd_target_pair_'+suffix] = df.groupby(['parent_molregno', 'tid_mutation'])['year'].transform('min')
    
    # first publication of pair with pchembl value
    df_mols_all_first_publication_pchembl = df[df['pchembl_value'].notnull()] \
            .groupby(['parent_molregno', 'tid_mutation'])['year'].min().reset_index() \
            .rename(columns={'year': 'first_publication_cpd_target_pair_w_pchembl_'+suffix})
    df = df.merge(df_mols_all_first_publication_pchembl, on=['parent_molregno', 'tid_mutation'], how = 'left')
    
    # return relevant summarised information without duplicates
    df = df[['parent_molregno', 'tid_mutation', 
            'pchembl_value_mean_'+suffix, 'pchembl_value_max_'+suffix, 'pchembl_value_median_'+suffix, 
            'first_publication_cpd_target_pair_'+suffix, 'first_publication_cpd_target_pair_w_pchembl_'+suffix]].drop_duplicates()
    return df



def get_aggregated_acticity_ct_pairs(chembl_con, limit_to_literature):
    """
    The following values are set to summarise the information for compound-target pairs:  

    |||
    | :----------- | :----------- |
    | *pchembl_value_mean* | mean pchembl value for a compound-target pair|
    | *pchembl_value_max*| maximum pchembl value for a compound-target pair|
    | *pchembl_value_median*| median pchembl value for a compound-target pair|
    | *first_publication_cpd_target_pair* | first publication in ChEMBL with this compound-target pair |
    | *first_publication_cpd_target_pair_w_pchembl* | first publication in ChEMBL with this compound-target pair and an associated pchembl value |

    The values are set for 
    - a subset of the dataset based on binding and functional assays (suffix '_BF') and 
    - a subset of the dataset set on only binding assays (suffix '_B'). 

    Therefore, there are two columns for each of the values above, one with the suffix '_BF' based on binding + functional data and one with the suffix '_B' based on only binding data.

    :param df_mols: _description_
    :type df_mols: _type_
    :return: _description_
    :rtype: _type_
    """
    df_mols = get_initial_dataset(chembl_con, limit_to_literature)

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
                                    on=['parent_molregno', 'tid_mutation'], how = 'left')
    # left merge because df_mols may contain assays that are of other types than binding / functional
    df_combined = df_combined.merge(df_mols.drop(columns=['pchembl_value', 'year', 'assay_type']).drop_duplicates(), 
                                    on=['parent_molregno', 'tid_mutation'], how = 'left')
    
    # TODO: include this?
    # test_utils.add_dataset_sizes(df_combined, "pre dm table", all_lengths, all_lengths_pchembl)
    return df_combined

