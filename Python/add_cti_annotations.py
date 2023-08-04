########### CTI (Compound-Target Interaction) Annotations ###########
def add_cti_annotations(df_combined, drug_mechanism_pairs_set, drug_mechanism_targets_set):
    """
    Every compound-target pair is assigned a DTI (drug target interaction) annotation.  

    The assignement is based on three questions:
    - Is the compound-target pair in the drug_mechanisms table? = Is it a known relevant compound-target interaction?
    - What is the max_phase of the compound? = Is it a drug / clinical compound?
    - Is the target in the drug_mechanisms table = Is it a therapeutic target?

    The assigments are based on the following table:

    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |in drug_mechanisms table?  |max_phase? |therapeutic target?    |DTI annotation |explanation                                            |
    +===========================+===========+=======================+===============+=======================================================+
    |yes                        |4          |--                     |D_DT           |drug - drug target                                     |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |yes                        |3          |--                     |C3_DT          |clinical candidate in phase 3 - drug target            |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |yes                        |2          |--                     |C2_DT          |clinical candidate in phase 2 - drug target            |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |yes                        |1          |--                     |C1_DT          |clinical candidate in phase 1 - drug target            |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |yes                        |<1         |--                     |C0_DT          |compound in unknown clinical phase\[1\] - drug target  |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |no                         |--         |yes                    |DT             |drug target                                            |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+
    |no                         |--         |no                     |NDT            |not drug target                                        |
    +---------------------------+-----------+-----------------------+---------------+-------------------------------------------------------+

    \[1\] There are three possible annotations in ChEMBL with max_phase not between 1 and 4:
    - 0.5 = early phase 1 clinical trials  
    - -1 = clinical phase unknown for drug or clinical candidate drug, i.e., where ChEMBL cannot assign a clinical phase
    - NULL = preclinical compounds with bioactivity data

    All three are grouped together into the annotation C0_DT.

    :param df_combined: _description_
    :type df_combined: _type_
    :param dti_tids_set: _description_
    :type dti_tids_set: _type_
    :param DTIs_set: _description_
    :type DTIs_set: _type_
    """
    # Add a new column *therapeutic_target* which is set to True if target is in the drug_mechanism table
    df_combined['therapeutic_target'] = df_combined['tid'].isin(drug_mechanism_targets_set)

    # Assign the annotations based on the table.
    # Compound-target pairs from the drug mechanism table
    df_combined.loc[(df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set) & (df_combined['max_phase'] == 4)), 'DTI'] = "D_DT"
    df_combined.loc[(df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set) & (df_combined['max_phase'] == 3)), 'DTI'] = "C3_DT"
    df_combined.loc[(df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set) & (df_combined['max_phase'] == 2)), 'DTI'] = "C2_DT"
    df_combined.loc[(df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set) & (df_combined['max_phase'] == 1)), 'DTI'] = "C1_DT"
    # Compounds that are in the drug_mechanism table but don't have a known phase between 1-4:
    df_combined.loc[(df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set) & 
                    (~df_combined['max_phase'].isin([1, 2, 3, 4]))), 'DTI'] = "C0_DT"

    # Target from the drug mechanism table
    df_combined.loc[((~df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set)) 
                    & (df_combined['therapeutic_target'] == True)), 'DTI'] = "DT"

    # Other compound-target pairs
    # if target is not a therapeutic target, 'cpd_target_pair' cannot be in DTIs_set
    # (~df_combined['cpd_target_pair'].isin(DTIs_set)) is included for clarity
    df_combined.loc[((~df_combined['cpd_target_pair'].isin(drug_mechanism_pairs_set)) 
                    & (df_combined['therapeutic_target'] == False)), 'DTI'] = "NDT"

    # TODO: include?
    # ############### TESTING: before discarding NDT rows ###############
    # add_dataset_sizes(df_combined, "pre DTI")

    # Discard rows that were annotated with NDT, 
    # i.e., compound-target pairs that are not in the drug_mechanisms table 
    # and for which the target was also not in the drug_mechanisms table (not a comparator compound).
    df_combined = df_combined[(df_combined['DTI'].isin(['D_DT', 'C3_DT', 'C2_DT', 'C1_DT', 'C0_DT', 'DT']))]

    # # TODO: include?
    # ############### TESTING: after discarding NDT rows ###############
    # add_dataset_sizes(df_combined, "post DTI")
    
    return df_combined





