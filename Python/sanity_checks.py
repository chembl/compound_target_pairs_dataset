def check_for_mixed_types(df_combined):
    # Check if there are mixed types in columns with dtype=object.
    # check that there are no mixed types in object columns
    issue_ctr = 0
    for i, (col, dtype) in enumerate(df_combined.dtypes.to_dict().items()):
        if dtype == object:
            a = set(df_combined[col])
            b = set(df_combined[col].astype(str))
            x = a-b
            y = b-a
            # is there a difference in the two sets
            if(len(a-b) > 0 or len(b-a) > 0):
                # is the difference something other than None being mapped to 'None' (string)?
                if len(x.difference({None})) > 0 or len(y.difference({'None'})) > 0:
                    print("Mixed types in column ", col)
                    print(a-b, '/', b-a)
                    issue_ctr += 1

    print("#Problems:", issue_ctr)


def check_null_values(df_combined):
    # Check if any columns contain nan or null which aren't recognised as null values. 
    # Do any columns have potential issues with null types?
    issue_ctr = 0
    for i, (col, dtype) in enumerate(df_combined.dtypes.to_dict().items()):
        if 'nan' in set(df_combined[df_combined[col].notnull()][col].astype(str)):
            print("Issue with nan in column", col)
            issue_ctr += 1
        if 'null' in set(df_combined[df_combined[col].notnull()][col].astype(str)):
            print("Issue with null in column", col)
            issue_ctr += 1

    print("#Problems:", issue_ctr)


def check_pairs_without_pchembl_are_in_drug_mechanisms(df_combined):
    # Check that rows without a pchembl value based on binding+functional assays (pchembl_x_BF) are in the drug_mechanism table.  
    # Note that this is not true for the pchembl_x_B columns based on binding data. 
    # They may be in the table because there is data based on functional assays but no data based on binding assays. 
    # All pchembl_value_x_BF columns without a pchembl should be in the dm table
    issue_ctr = 0
    for pchembl in ['pchembl_value_mean_BF', 'pchembl_value_max_BF', 'pchembl_value_median_BF']:
        if not df_combined[(df_combined[pchembl].isnull())].equals(
            df_combined[(df_combined['in_dm_table'] == True) & (df_combined[pchembl].isnull())]):
            print("Problem with", pchembl)
            issue_ctr += 1
            
    print("#Problems:", issue_ctr)

def check_ligand_efficiency_metrics(df_combined):
    # Check that ligand efficiency metrics are only null when at least one of the values used to calculate them is null.
    # Ligand efficiency metrics are only null when at least one of the values used to calculate them is null.
    issue_ctr = 0
    for suffix in ['BF', 'B']:
        if not df_combined[(df_combined['LE_'+suffix].isnull())].equals(
        df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull()) 
                    | (df_combined['heavy_atoms'].isnull()) 
                    | (df_combined['heavy_atoms'] == 0)]):
            print("Problem with LE_"+suffix)
            issue_ctr += 1

        if not df_combined[(df_combined['BEI_'+suffix].isnull())].equals(
        df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull()) 
                    | (df_combined['mw_freebase'].isnull()) 
                    | (df_combined['mw_freebase'] == 0)]):
            print("Problem with BEI_"+suffix)
            issue_ctr += 1

        if not df_combined[(df_combined['SEI_'+suffix].isnull())].equals(
        df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull()) 
                    | (df_combined['psa'].isnull())
                    | (df_combined['psa'] == 0)]):
            print("Problem with SEI_"+suffix)
            issue_ctr += 1
            
        if not df_combined[(df_combined['LLE_'+suffix].isnull())].equals(
        df_combined[(df_combined['pchembl_value_mean_'+suffix].isnull()) 
                    | (df_combined['alogp'].isnull())]):
            print("Problem with LLE_"+suffix)
            issue_ctr += 1
            
    print("#Problems:", issue_ctr)

def check_compound_props(df_combined, df_cpd_props):
    # Check that compound props are only null 
    # if the property in the parent_molregno is not in the compound props table or if the value in the compound props table is null.
    # Check for issues with compound properties that are null
    issue_ctr = 0

    # missing values because the parent_molregno is not in the compound props table
    no_cpd_prop_info = len(df_combined[~df_combined['parent_molregno'].isin(set(df_cpd_props['parent_molregno']))])

    for col in df_cpd_props.columns:
        if col != 'parent_molregno':
            # missing values because the compound props query returns null (exists but is null)
            missing_values = len(df_combined[df_combined['parent_molregno'].isin(set(df_cpd_props[df_cpd_props[col].isnull()]['parent_molregno']))])
            null_values = no_cpd_prop_info+missing_values
            if null_values != len(df_combined[df_combined[col].isnull()]):
                print("Problem with column", col)
                issue_ctr += 1
            
    print("#Problems:", issue_ctr)


def check_atc_and_target_classes(df_combined, atc_levels, target_classes_level1, target_classes_level2):
    # Check that atc_level1 and target class information is only null if the parent_molregno / target id is not in the respective table.
    # issues with atc or target classes
    issue_ctr = 0

    if not df_combined[(df_combined['atc_level1'].isnull())].equals(
        df_combined[~df_combined['parent_molregno'].isin(set(atc_levels['parent_molregno']))]):
        print("Problem with atc_level1")
        issue_ctr += 1
        
    if not df_combined[(df_combined['target_class_l1'].isnull())].equals(
        df_combined[~df_combined['tid'].isin(set(target_classes_level1['tid']))]):
        print("Problem with target_class_l1")
        issue_ctr += 1
        
    if not df_combined[(df_combined['target_class_l2'].isnull())].equals(
        df_combined[~df_combined['tid'].isin(set(target_classes_level2['tid']))]):
        print("Problem with target_class_l2")
        issue_ctr += 1

    print("#Problems:", issue_ctr)

def check_rdkit_props(df_combined):
    # Check that columns set by the RDKit are only null if there is no canonical SMILES for the molecule.  
    # Scaffolds are excluded from this test because they can be None if the molecule is acyclic. 
    # issues with RDKit methods
    issue_ctr = 0

    for col in ['fraction_csp3', 'num_aliphatic_carbocycles', 'num_aliphatic_heterocycles', 'num_aliphatic_rings', 
                'num_aromatic_carbocycles', 'num_aromatic_heterocycles', 'num_aromatic_rings', 
                'num_heteroatoms', 'num_saturated_carbocycles', 'num_saturated_heterocycles', 
                'num_saturated_rings', 'ring_count', 'num_stereocentres',
                'aromatic_atoms', 'aromatic_c', 'aromatic_n', 'aromatic_hetero']:
        if len(df_combined[df_combined[col].isnull()]) != len(df_combined[df_combined['canonical_smiles'].isnull()].copy()):
            print("Problem with ", col)
            issue_ctr += 1

    print("#Problems:", issue_ctr)



def sanity_checks(df_combined, df_cpd_props, atc_levels, target_classes_level1, target_classes_level2, calculate_RDKit):
    check_for_mixed_types(df_combined)
    check_null_values(df_combined)
    check_pairs_without_pchembl_are_in_drug_mechanisms(df_combined)
    check_ligand_efficiency_metrics(df_combined)
    check_compound_props(df_combined, df_cpd_props)
    check_atc_and_target_classes(df_combined, atc_levels, target_classes_level1, target_classes_level2)
    if calculate_RDKit:
        check_rdkit_props(df_combined)
