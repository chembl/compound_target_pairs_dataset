import pandas as pd

########### Add Target Class Annotations Based on ChEMBL Data ###########
def add_chembl_target_class_annotations(df_combined, chembl_con):
    # Add information about level 1 and level 2 target class annotations in ChEMBL.
    sql = '''
    SELECT DISTINCT tc.tid, 
        pc.protein_class_id, pc.pref_name, pc.short_name, pc.protein_class_desc, pc.definition
    FROM protein_classification pc
    -- join several tables to get the corresponding target id
    INNER JOIN component_class cc
        ON pc.protein_class_id = cc.protein_class_id
    INNER JOIN component_sequences cs
        ON cc.component_id = cs.component_id
    INNER JOIN target_components tc
        ON cs.component_id = tc.component_id
    '''

    df_target_classes = pd.read_sql_query(sql, con=chembl_con)
    # only interested in the target ids that are in the current dataset
    current_tids = set(df_combined['tid'])
    df_target_classes = df_target_classes[df_target_classes['tid'].isin(current_tids)]

    # Query the protein_classification table for the protein classification hierarchy 
    # and merge it with the target class information for specific tids.
    sql = '''
    WITH RECURSIVE pc_hierarchy AS (
        SELECT protein_class_id,
                parent_id,
                class_level,
                pref_name AS names
        FROM protein_classification
        WHERE parent_id IS NULL

        UNION ALL
    
        SELECT pc.protein_class_id,
            pc.parent_id,
            pc.class_level,
            -- recursively add current protein classification pref_name to string, separated by |
            pc_hierarchy.names || '|' || pc.pref_name 
        FROM protein_classification pc, pc_hierarchy
        WHERE pc.parent_id = pc_hierarchy.protein_class_id
    )
    SELECT *
    FROM pc_hierarchy
    '''

    target_class_hierarchy = pd.read_sql_query(sql, con=chembl_con)
    target_class_hierarchy[['l0', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6']] = target_class_hierarchy['names'].str.split('|', expand=True)
    target_class_hierarchy = target_class_hierarchy[target_class_hierarchy['protein_class_id'] != 0][['protein_class_id', 'l1', 'l2']]
    df_target_classes = df_target_classes.merge(target_class_hierarchy, on = 'protein_class_id', how = 'left')



    # Summarise the information for a target id with several assigned target classes of level 1 into one description. If a target id has more than one assigned target class, the target class 'Unclassified protein' is discarded.
    level = 'l1'
    between_str_join = '|'
    target_classes_level1 = df_target_classes[['tid', level]].drop_duplicates().dropna()

    # remove 'Unclassified protein' from targets with more than one target class, level 1
    nof_classes = target_classes_level1.groupby(['tid'])[level].count()
    target_classes_level1 = target_classes_level1[
        (target_classes_level1['tid'].isin(nof_classes[nof_classes == 1].index.tolist())) 
        | ((target_classes_level1['tid'].isin(nof_classes[nof_classes > 1].index.tolist())) 
        & (target_classes_level1['l1'] != 'Unclassified protein'))]

    target_classes_level1['target_class_l1'] = target_classes_level1.groupby(['tid'])[level].transform(lambda x: between_str_join.join(sorted(x)))
    target_classes_level1 = target_classes_level1[['tid', 'target_class_l1']].drop_duplicates()

    df_combined = df_combined.merge(target_classes_level1, on='tid', how = 'left')



    # Repeat the summary step for target classes of level 2.
    level = 'l2'
    target_classes_level2 = df_target_classes[['tid', level]].drop_duplicates().dropna()
    target_classes_level2['target_class_l2'] = target_classes_level2.groupby(['tid'])[level].transform(lambda x: between_str_join.join(sorted(x)))
    target_classes_level2 = target_classes_level2[['tid', 'target_class_l2']].drop_duplicates()

    df_combined = df_combined.merge(target_classes_level2, on='tid', how = 'left')


    # # TODO: output targets with more than one target class to reassign by hand if wanted
    # # Instances with targets with more than one target class assigned to them.  
    # # These could be reassigned by hand if a single target class is preferable.
    # ############### TESTING: which targets have more than one level 1 target class assigned to them? ###############
    # test = df_combined[(df_combined['target_class_l1'].notnull()) & (df_combined['target_class_l1'].str.contains('|', regex=False))][['tid', 'target_pref_name', 'target_type', 'target_class_l1', 'target_class_l2']].drop_duplicates()
    # print("#Instances with >1 level 1 target class:", len(test))
    # test
    # ############### TESTING: which targets have more than one level 2 target class assigned to them? ###############
    # test = df_combined[(df_combined['target_class_l2'].notnull()) & (df_combined['target_class_l2'].str.contains('|', regex=False))][['tid', 'target_pref_name', 'target_type', 'target_class_l1', 'target_class_l2']].drop_duplicates()
    # print("#Instances with >1 level 2 target class:", len(test))
    # test

    return df_combined, target_classes_level1, target_classes_level2

