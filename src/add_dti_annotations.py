"""
Add DTI (Drug-Target Interaction) Annotations to the dataset.
"""

from dataset import Dataset


########### DTI (Drug-Target Interaction) Annotations ###########
def add_dti_annotations(
    dataset: Dataset,
):
    """
    Every compound-target pair is assigned a DTI (drug target interaction) annotation.  

    The assignment is based on three questions:

    - Is the compound-target pair in the drug_mechanisms table? = \
        Is it a known relevant compound-target interaction?
    - What is the max_phase of the compound? = Is it a drug / clinical compound?
    - Is the target in the drug_mechanisms table = Is it a therapeutic target?

    The assigments are based on the following table:

    +------------+----------+-----------+-----------+---------------------------------------------+
    |in DM table?|max_phase?|th. target?|DTI        |explanation                                  |
    +============+==========+===========+===========+=============================================+
    | yes        | 4        | --        | D_DT [#]_ | drug - drug target                          |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | yes        | 3        | --        | C3_DT     | clinical candidate in phase 3 - drug target |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | yes        | 2        | --        | C2_DT     | clinical candidate in phase 2 - drug target |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | yes        | 1        | --        | C1_DT     | clinical candidate in phase 1 - drug target |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | yes        | <1       | --        | C0_DT     |compound in unknown phase [#]_ - drug target |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | no         | --       | yes       | DT        | drug target                                 |
    +------------+----------+-----------+-----------+---------------------------------------------+
    | no         | --       | no        | NDT       | not drug target                             |
    +------------+----------+-----------+-----------+---------------------------------------------+

    .. [#] The annotation D_DT instead of C4_DT was chosen to be consistent \
    with the annotations in a previous version of the dataset. \\ \
    For the same reason the column is named DTI (drug-target interaction) \
    instead of CTI (compound-target interaction) \
    despite having specific annotations for clinical canidates.

    .. [#] C0_DT groups together all compounds with a max_phase not between 1 and 4. 
    
    Since ChEMBL32 there are three possible annotations in ChEMBL \
    with a max_phase value not between 1 and 4:

    - 0.5 = early phase 1 clinical trials  
    - -1 = clinical phase unknown for drug or clinical candidate drug, \
            i.e., where ChEMBL cannot assign a clinical phase
    - NULL = preclinical compounds with bioactivity data

    All three are grouped together into the annotation C0_DT.

    Compound-target pairs that were annotated with NDT, \
    i.e., compound-target pairs that are not in the drug_mechanisms table \
    and for which the target was also not in the drug_mechanisms table \
    (not a comparator compound), are discarded.

    :param dataset: Dataset with all relevant information:
        - Pandas DataFrame with compound-target pairs \
            based on activities AND drug_mechanism table
        - set of compound-target pairs in the drug_mechanism table
        - set of targets in the drug_mechanism table
    :type dataset: Dataset
    """
    # Add a new column *therapeutic_target* which is set to True
    # if the target is in the drug_mechanism table
    dataset.df_result["therapeutic_target"] = dataset.df_result["tid"].isin(
        dataset.drug_mechanism_targets_set
    )

    # Assign the annotations based on the table.
    # Compound-target pairs from the drug mechanism table
    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)
            & (dataset.df_result["max_phase"] == 4)
        ),
        "DTI",
    ] = "D_DT"

    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)
            & (dataset.df_result["max_phase"] == 3)
        ),
        "DTI",
    ] = "C3_DT"

    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)
            & (dataset.df_result["max_phase"] == 2)
        ),
        "DTI",
    ] = "C2_DT"

    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)
            & (dataset.df_result["max_phase"] == 1)
        ),
        "DTI",
    ] = "C1_DT"

    # Compounds that are in the drug_mechanism table but don't have a known phase between 1-4:
    dataset.df_result.loc[
        (
            dataset.df_result["cpd_target_pair"].isin(dataset.drug_mechanism_pairs_set)
            & (~dataset.df_result["max_phase"].isin([1, 2, 3, 4]))
        ),
        "DTI",
    ] = "C0_DT"

    # Target is in the drug mechanism table
    dataset.df_result.loc[
        (
            (
                ~dataset.df_result["cpd_target_pair"].isin(
                    dataset.drug_mechanism_pairs_set
                )
            )
            & (dataset.df_result["therapeutic_target"])
        ),
        "DTI",
    ] = "DT"

    # Other compound-target pairs
    # if target is not a therapeutic target, 'cpd_target_pair' cannot be in DTIs_set
    # (~dataset.df_result['cpd_target_pair'].isin(DTIs_set)) is included for clarity
    dataset.df_result.loc[
        (
            (
                ~dataset.df_result["cpd_target_pair"].isin(
                    dataset.drug_mechanism_pairs_set
                )
            )
            & ~(dataset.df_result["therapeutic_target"])
        ),
        "DTI",
    ] = "NDT"

    # Discard NDT rows
    dataset.df_result = dataset.df_result[
        (
            dataset.df_result["DTI"].isin(
                ["D_DT", "C3_DT", "C2_DT", "C1_DT", "C0_DT", "DT"]
            )
        )
    ]
