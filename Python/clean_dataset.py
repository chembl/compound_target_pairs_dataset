import pandas as pd

def remove_irrelevant_compounds(df_combined, chembl_con):
    # TODO: REFACTOR
    # Remove Compounds Without a Smiles and Mixtures
    # Remove Compounds With a Smiles Containing a '.'
    # Double-check that rows with a SMILES containing a '.' are the parent structures, 
    # i.e., there was no error in using salt information instead of parent information.  
    # These compounds are salts or mixtures and will be removed in the next step.
    sql = '''
    SELECT DISTINCT mh.molregno as salt_molregno, mh.parent_molregno
    FROM molecule_hierarchy mh
    '''
    df_hierarchy = pd.read_sql_query(sql, con=chembl_con)

    smiles_with_dot = df_combined[df_combined['canonical_smiles'].notnull() & df_combined['canonical_smiles'].str.contains('.', regex=False)]
    smiles_with_dot = smiles_with_dot[['canonical_smiles', 'parent_molregno']].drop_duplicates()

    issue_ctr = 0
    for parent_molregno in set(smiles_with_dot['parent_molregno']):
        # the molrego should occur at least once as a parent_molregno
        if not len(df_hierarchy[df_hierarchy['parent_molregno'] == parent_molregno]) > 0:
            print(df_hierarchy[df_hierarchy['parent_molregno'] == parent_molregno])
            issue_ctr += 1
        # if it occurs as a salt_molregno, its parent_molregno should be identical, 
        # i.e., the molregno is a parent_molregno
        df_salt_molregno = df_hierarchy[df_hierarchy['salt_molregno'] == parent_molregno]
        if not df_salt_molregno['salt_molregno'].equals(df_salt_molregno['parent_molregno']):
            print(df_hierarchy[df_hierarchy['salt_molregno'] == parent_molregno])
            issue_ctr += 1

    print("#Problems:", issue_ctr)

    sql = '''
    SELECT DISTINCT mh.parent_molregno, struct.canonical_smiles
    FROM molecule_hierarchy mh
    INNER JOIN compound_structures struct
        ON mh.parent_molregno = struct.molregno
    '''
    df_parent_smiles = pd.read_sql_query(sql, con=chembl_con)

    issue_ctr = 0
    for parent_molregno in set(smiles_with_dot['parent_molregno']):
        df_smiles = df_parent_smiles[df_parent_smiles['parent_molregno'] == parent_molregno]['canonical_smiles'].item()
        smiles_w_dot = smiles_with_dot[smiles_with_dot['parent_molregno'] == parent_molregno]['canonical_smiles'].item()
        if df_smiles != smiles_w_dot:
            print(parent_molregno)
            print(df_smiles)
            print(smiles_w_dot)
            print()
            issue_ctr += 1

    print("#Problems:", issue_ctr)

    # Remove rows that contain a SMILES with a dot or that don't have a SMILES.
    len_missing_smiles = len(df_combined[df_combined['canonical_smiles'].isnull()])
    len_smiles_w_dot = len(df_combined[df_combined['parent_molregno'].isin(set(smiles_with_dot['parent_molregno']))])
    print("{:30} {}".format("#Rows w/o SMILES:", len_missing_smiles))
    print("{:30} {}".format("#Rows w SMILES with dot:", len_smiles_w_dot))
    print("{:30} {}".format("Predicted size after removal:", len(df_combined)-len_missing_smiles-len_smiles_w_dot))
    df_combined = df_combined[(df_combined['canonical_smiles'].notnull()) & ~(df_combined['parent_molregno'].isin(set(smiles_with_dot['parent_molregno'])))]
    print("{:30} {}".format("Size:", len(df_combined)))

    return df_combined


def clean_none_values(df_combined):
    # Change nan values and empty strings to None for consistency. 
    # Change all None / nan values to None
    df_combined = df_combined.where(pd.notnull(df_combined), None)
    # replace empty strings with None
    df_combined = df_combined.replace('', None).reset_index(drop=True)

    return df_combined

def set_types(df_combined, calculate_RDKit):
    # Set relevant columns to int
    df_combined = df_combined.astype({
        'first_approval': 'Int64',
        'usan_year': 'Int64',
        'first_publication_cpd_target_pair_BF': 'Int64',
        'first_publication_cpd_target_pair_w_pchembl_BF': 'Int64',
        'first_publication_cpd_target_pair_B': 'Int64',
        'first_publication_cpd_target_pair_w_pchembl_B': 'Int64',
        'first_publication_cpd': 'Int64',
        'hba': 'Int64',
        'hbd': 'Int64',
        'rtb': 'Int64',
        'num_ro5_violations': 'Int64',
        'aromatic_rings': 'Int64',
        'heavy_atoms': 'Int64',
        'hba_lipinski': 'Int64',
        'hbd_lipinski': 'Int64',
        'num_lipinski_ro5_violations': 'Int64'
    })

    if calculate_RDKit:
        df_combined = df_combined.astype({
            'num_aliphatic_carbocycles': 'Int64',
            'num_aliphatic_heterocycles': 'Int64',
            'num_aliphatic_rings': 'Int64',
            'num_aromatic_carbocycles': 'Int64',
            'num_aromatic_heterocycles': 'Int64',
            'num_aromatic_rings': 'Int64',
            'num_heteroatoms': 'Int64',
            'num_saturated_carbocycles': 'Int64',
            'num_saturated_heterocycles': 'Int64',
            'num_saturated_rings': 'Int64',
            'ring_count': 'Int64',
            'num_stereocentres': 'Int64',
            'aromatic_atoms': 'Int64',
            'aromatic_c': 'Int64',
            'aromatic_n': 'Int64',
            'aromatic_hetero': 'Int64'
        })

    return df_combined

def round_floats(df_combined, decimal_places = 4):
    # Round float columns to <decimal_places> decimal places
    for i, (col, dtype) in enumerate(df_combined.dtypes.to_dict().items()):
        if ((dtype == 'float64') or (dtype == 'Float64')) and col != 'max_phase':
            df_combined[col] = df_combined[col].round(decimals=decimal_places)
    
    return df_combined

def reorder_columns(df_combined, calculate_RDKit):
    len_columns_before = len(df_combined.columns)
    if calculate_RDKit:
        df_combined = df_combined[['parent_molregno', 'parent_chemblid', 'parent_pref_name', 
                                'max_phase', 'first_approval', 'usan_year', 'black_box_warning', 
                                'prodrug', 'oral', 'parenteral', 'topical',
                                'tid', 'mutation', 'target_chembl_id', 'target_pref_name', 'target_type', 
                                'organism', 'tid_mutation',
                                'cpd_target_pair', 'cpd_target_pair_mutation',
                                'pchembl_value_mean_BF', 'pchembl_value_max_BF', 'pchembl_value_median_BF',
                                'first_publication_cpd_target_pair_BF', 'first_publication_cpd_target_pair_w_pchembl_BF',
                                'pchembl_value_mean_B', 'pchembl_value_max_B', 'pchembl_value_median_B',
                                'first_publication_cpd_target_pair_B', 'first_publication_cpd_target_pair_w_pchembl_B', 
                                'therapeutic_target', 'DTI',
                                'first_publication_cpd', 'mw_freebase', 'alogp', 'hba', 'hbd', 'psa',
                                'rtb', 'ro3_pass', 'num_ro5_violations', 'cx_most_apka', 'cx_most_bpka',
                                'cx_logp', 'cx_logd', 'molecular_species', 'full_mwt', 'aromatic_rings',
                                'heavy_atoms', 'qed_weighted', 'mw_monoisotopic', 'full_molformula',
                                'hba_lipinski', 'hbd_lipinski', 'num_lipinski_ro5_violations', 
                                'standard_inchi', 'standard_inchi_key', 'canonical_smiles', 
                                'LE_B', 'BEI_B', 'SEI_B', 'LLE_B',
                                'LE_BF', 'BEI_BF', 'SEI_BF', 'LLE_BF',
                                'atc_level1', 'target_class_l1', 'target_class_l2', 
                                'fraction_csp3', 
                                'num_aliphatic_carbocycles', 'num_aliphatic_heterocycles', 'num_aliphatic_rings', 
                                'num_aromatic_carbocycles', 'num_aromatic_heterocycles', 'num_aromatic_rings',
                                'num_heteroatoms', 
                                'num_saturated_carbocycles', 'num_saturated_heterocycles', 'num_saturated_rings', 
                                'ring_count', 'num_stereocentres', 
                                'aromatic_atoms', 'aromatic_c', 'aromatic_n', 'aromatic_hetero', 
                                'scaffold_w_stereo', 'scaffold_wo_stereo',
                                'in_dm_table', 'keep_for_binding']]
    else:
        df_combined = df_combined[['parent_molregno', 'parent_chemblid', 'parent_pref_name', 
                                'max_phase', 'first_approval', 'usan_year', 'black_box_warning', 
                                'prodrug', 'oral', 'parenteral', 'topical',
                                'tid', 'mutation', 'target_chembl_id', 'target_pref_name', 'target_type', 
                                'organism', 'tid_mutation',
                                'cpd_target_pair', 'cpd_target_pair_mutation',
                                'pchembl_value_mean_BF', 'pchembl_value_max_BF', 'pchembl_value_median_BF',
                                'first_publication_cpd_target_pair_BF', 'first_publication_cpd_target_pair_w_pchembl_BF',
                                'pchembl_value_mean_B', 'pchembl_value_max_B', 'pchembl_value_median_B',
                                'first_publication_cpd_target_pair_B', 'first_publication_cpd_target_pair_w_pchembl_B', 
                                'therapeutic_target', 'DTI',
                                'first_publication_cpd', 'mw_freebase', 'alogp', 'hba', 'hbd', 'psa',
                                'rtb', 'ro3_pass', 'num_ro5_violations', 'cx_most_apka', 'cx_most_bpka',
                                'cx_logp', 'cx_logd', 'molecular_species', 'full_mwt', 'aromatic_rings',
                                'heavy_atoms', 'qed_weighted', 'mw_monoisotopic', 'full_molformula',
                                'hba_lipinski', 'hbd_lipinski', 'num_lipinski_ro5_violations', 
                                'standard_inchi', 'standard_inchi_key', 'canonical_smiles', 
                                'LE_B', 'BEI_B', 'SEI_B', 'LLE_B',
                                'LE_BF', 'BEI_BF', 'SEI_BF', 'LLE_BF',
                                'atc_level1', 'target_class_l1', 'target_class_l2', 
                                'in_dm_table', 'keep_for_binding']]

    len_columns_after = len(df_combined.columns)
    assert(len_columns_before == len_columns_after), "Different number of columns after reordering."

    return df_combined

def clean_dataset(df_combined, chembl_con, calculate_RDKit):
    df_combined = clean_none_values(df_combined)
    df_combined = set_types(df_combined, calculate_RDKit)
    df_combined = round_floats(df_combined, decimal_places = 4)
    df_combined = reorder_columns(df_combined, calculate_RDKit)
    return df_combined



