import pandas as pd

########### Extract Drug-Target Interactions With Disease Relevance From the drug_mechanism Table ###########
def get_drug_mechanisms_interactions(chembl_con):
    """
    Extract the known drug-target interactions from ChEMBL (these include some interactions between compounds with a max_phase < 4 and targets). These will be used to determine if compound-target pairs from the activities query above are known compound-target interactions. 

    Note: Compound-target pairs can be in the drug_mechanisms table even though the compound is not a drug (max_phase < 4). For ease of writing, these will be referred to as drug-target interactions as well rather than compound-target pairs with a known disease-relevant interaction. 

    Only entries with a disease_efficacy of 1 are taken into account, i.e., the target is believed to play a role in the efficacy of the drug.  
    *disease_efficacy: Flag to show whether the target assigned is believed to play a role in the efficacy of the drug in the indication(s) for which it is approved (1 = yes, 0 = no).*

    :param chembl_con: _description_
    :type chembl_con: _type_
    :return: _description_
    :rtype: _type_
    """
    sql = '''
    SELECT DISTINCT mh.parent_molregno, dm.tid
    FROM drug_mechanism dm
    INNER JOIN molecule_hierarchy mh
        ON dm.molregno = mh.molregno
    INNER JOIN molecule_dictionary md
        ON mh.parent_molregno = md.molregno
    WHERE dm.disease_efficacy = 1
        and dm.tid is not null
    '''

    df_dti = pd.read_sql_query(sql, con=chembl_con)
    return df_dti

def get_related_tids(chembl_con):
    """
    Query target_relations for related target ids to increase the number of target ids for which there is data in the drug_mechanisms table.
    The following mappings are considered:

    +-------------------------------+-----------------------+---------------+
    |protein family                 |-[superset of]->       | single protein|
    +-------------------------------+-----------------------+---------------+
    |protein complex                |-[superset of]->       | single protein|
    +-------------------------------+-----------------------+---------------+
    |protein complex group          |-[superset of]->       | single protein|
    +-------------------------------+-----------------------+---------------+
    |single protein                 |-[equivalent to]->     | single protein|
    +-------------------------------+-----------------------+---------------+
    |chimeric protein               |-[superset of]->       | single protein|
    +-------------------------------+-----------------------+---------------+
    |protein-protein interaction    |-[superset of]->       | single protein|
    +-------------------------------+-----------------------+---------------+

    For example, for *protein family -[superset of]-> single protein* this means:  
    If there is a known relevant interaction between a compound and a protein family, interactions between the compound and single proteins of that protein family are considered to be known interactions as well.

    :param chembl_con: _description_
    :type chembl_con: _type_
    :return: _description_
    :rtype: _type_
    """


    sql = '''
    SELECT tr.tid, tr.relationship, tr.related_tid, 
        td1.pref_name as pref_name_1, td1.target_type as target_type_1, td1.organism as organism_1, 
        td2.pref_name as pref_name_2, td2.target_type as target_type_2, td2.organism as organism_2 
    FROM target_relations tr
    INNER JOIN target_dictionary td1
        ON tr.tid = td1.tid
    INNER JOIN target_dictionary td2
        ON tr.related_tid = td2.tid
    '''

    df_related_targets = pd.read_sql_query(sql, con=chembl_con)
    return df_related_targets

def get_relevant_mappings(chembl_con):
    df_related_targets = get_related_tids(chembl_con)

    protein_family_mapping = df_related_targets[(df_related_targets["target_type_1"] == "PROTEIN FAMILY") 
                    & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                    & (df_related_targets["relationship"] == "SUPERSET OF")]

    protein_complex_mapping = df_related_targets[(df_related_targets["target_type_1"] == "PROTEIN COMPLEX") 
                        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                        & (df_related_targets["relationship"] == "SUPERSET OF")]

    protein_complex_group_mapping = df_related_targets[(df_related_targets["target_type_1"] == "PROTEIN COMPLEX GROUP") 
                        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                        & (df_related_targets["relationship"] == "SUPERSET OF")]

    single_protein_mapping = df_related_targets[(df_related_targets["target_type_1"] == "SINGLE PROTEIN") 
                        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                        & (df_related_targets["relationship"] == "EQUIVALENT TO")]

    chimeric_protein_mapping = df_related_targets[(df_related_targets["target_type_1"] == "CHIMERIC PROTEIN") 
                        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                        & (df_related_targets["relationship"] == "SUPERSET OF")]

    ppi_mapping = df_related_targets[(df_related_targets["target_type_1"] == "PROTEIN-PROTEIN INTERACTION") 
                        & (df_related_targets["target_type_2"] == "SINGLE PROTEIN")
                        & (df_related_targets["relationship"] == "SUPERSET OF")]

    relevant_mappings = pd.concat([protein_family_mapping, 
                                protein_complex_mapping, 
                                protein_complex_group_mapping,
                                single_protein_mapping, 
                                chimeric_protein_mapping, 
                                ppi_mapping])
    return relevant_mappings


def get_cti_from_drug_mechanisms(chembl_con):
    df_dti = get_drug_mechanisms_interactions(chembl_con)
    relevant_mappings = get_relevant_mappings(chembl_con)

    # drug-target-interactions (DTI) and target ids (dti_tids) based on the drug_mechanisms table
    DTIs_original = set(df_dti.agg('{0[parent_molregno]}_{0[tid]}'.format, axis=1))
    dti_tids_original = set(df_dti['tid'])

    # drug-target-interactions (DTI) and target ids (dti_tids) based on mapped target ids
    df_dti_mapped_targets = df_dti.merge(relevant_mappings, on = 'tid', how = 'inner')
    DTIs_mapped = set(df_dti_mapped_targets.agg('{0[parent_molregno]}_{0[related_tid]}'.format, axis=1))
    dti_tids_mapped = set(df_dti_mapped_targets['related_tid'])

    # combined drug-target-interactions (DTI) and target ids (dti_tids) 
    # based on drug_mechanisms table and mapped target ids
    DTIs_set = DTIs_original.union(DTIs_mapped)
    dti_tids_set = dti_tids_original.union(dti_tids_mapped)

    return df_dti, df_dti_mapped_targets, dti_tids_set, DTIs_set


########### Add Compounds From the drug_mechanism Table to the Dataset ###########
def add_cti_from_drug_mechanisms(chembl_con, df_combined):
    df_dti, df_dti_mapped_targets, dti_tids_set, DTIs_set = get_cti_from_drug_mechanisms(chembl_con)

    # Add compound-target pairs from the drug_mechanism table that are not in the dataset based on the initial ChEMBL query.
    # These are compound-target pairs for which there is no associated pchembl value data.
    # Since the pairs are known interactions, they are added to the dataset despite not having a pchembl value.
    cpd_target_pairs = pd.concat([df_dti[['parent_molregno', 'tid']], 
                              df_dti_mapped_targets[['parent_molregno', 'related_tid']]
                              .rename(columns={'related_tid': 'tid'})]).drop_duplicates()

    # Set columns existing in the df_combined table.
    # None of the targets from the drug mechanism table have any mutation annotation, hence tid_mutation = tid
    cpd_target_pairs['tid_mutation'] = cpd_target_pairs['tid'].astype('str')
    cpd_target_pairs['cpd_target_pair'] = cpd_target_pairs.agg('{0[parent_molregno]}_{0[tid]}'.format, axis=1)
    cpd_target_pairs['cpd_target_pair_mutation'] = cpd_target_pairs.agg('{0[parent_molregno]}_{0[tid_mutation]}'.format, axis=1)

    # Add a new column *in_dm_table* which is set to True if the compound target pair (taking mutation annotations into account) is in the drug_mechanism table. 
    # New column: is the compound target pair (taking mutation annotations into account) in the drug_mechanism table?
    cpd_target_pairs['in_dm_table'] = True

    # Set *in_dm_table* for the initial dataset based on the ChEMBL query (df_combined).
    df_combined['in_dm_table'] = False
    df_combined.loc[(df_combined['cpd_target_pair_mutation'].isin(set(cpd_target_pairs['cpd_target_pair_mutation']))), 
                    'in_dm_table'] = True
    
    # Limit the pairs to the ones that are not yet in the dataset.  
    # Mutation annotations are taken into account. 
    # Therefore, *(cpd A, target B without mutation)* will be added if a pchembl is present for *(cpd A, target B with mutation C)* but not for *(cpd A, target B without mutation)*.
    # pairs for which there is no information based in binding or functional assays in the original ChEMBL query
    cpd_target_pairs = cpd_target_pairs[~(cpd_target_pairs['cpd_target_pair_mutation'].isin(set(df_combined['cpd_target_pair_mutation'])))].copy()
    print("#Pairs not yet present based on binding or functional assays:", len(cpd_target_pairs))

    # Query compound and target information and combine it with the new compound-target pairs table.
    sql = '''
    SELECT md.molregno as parent_molregno, 
        md.chembl_id as parent_chemblid, md.pref_name as parent_pref_name,
        md.max_phase, md.first_approval, md.usan_year, md.black_box_warning, 
        md.prodrug, md.oral, md.parenteral, md.topical
    FROM molecule_dictionary md
    '''

    df_compound_info = pd.read_sql_query(sql, con=chembl_con)
    cpd_target_pairs = cpd_target_pairs.merge(df_compound_info, on = 'parent_molregno', how = 'left')

    sql = '''
    SELECT td.tid, td.chembl_id as target_chembl_id, td.pref_name as target_pref_name, td.target_type, td.organism
    FROM target_dictionary td
    '''

    df_target_info = pd.read_sql_query(sql, con=chembl_con)
    # Fix problems with null not being recognised as None
    df_target_info.loc[df_target_info['organism'].astype(str) == 'null', 'organism'] = None
    cpd_target_pairs = cpd_target_pairs.merge(df_target_info, on = 'tid', how = 'left')

    # Combined data of existing query with new compound-target pairs.
    df_combined = pd.concat([df_combined, cpd_target_pairs]) 

    # Add a new column *keep_for_binding* which is set to True if the row should be kept 
    # if you want to limit the dataset to only data based on binding assays.   
    # Rows are kept if 
    # - there is a binding data-based pchembl value or
    # - the compound-target pair is in the drug_mechanism table
    df_combined['keep_for_binding'] = False
    df_combined.loc[((df_combined['pchembl_value_mean_B'].notnull()) | 
                    (df_combined['in_dm_table'] == True )), 'keep_for_binding'] = True
    
    return df_combined, dti_tids_set, DTIs_set