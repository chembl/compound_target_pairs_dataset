import numpy as np
import pandas as pd
import sqlite3

########### Add Compound Properties Based on ChEMBL Data ###########
def add_first_publication_date(df_combined: pd.DataFrame, chembl_con: sqlite3.Connection, limit_to_literature: bool) -> pd.DataFrame:        
    """
    Query and calculate the first publication of a compound based on ChEMBL data (column name: first_publication_cpd).
    If limit_to_literature is True, this corresponds to the first appearance of the compound in the literature according to ChEMBL. 
    Otherwise this is the first appearance in any source in ChEMBL. 

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Base first_publication_cpd on literature sources only if True. 
    :type limit_to_literature: bool
    :return: Pandas DataFrame with added first_publication_cpd.
    :rtype: pd.DataFrame
    """
    # information about salts is aggregated in the parent
    sql = '''
    SELECT DISTINCT docs.year, mh.parent_molregno
    FROM docs
    LEFT JOIN compound_records cr
        ON docs.doc_id = cr.doc_id
    INNER JOIN molecule_hierarchy mh 
        ON cr.molregno = mh.molregno   -- cr.molregno = salt_molregno
    WHERE docs.year is not null
    '''
    if limit_to_literature:
        sql += '''    and docs.src_id = 1'''
    df_docs = pd.read_sql_query(sql, con=chembl_con)
    
    df_docs['first_publication_cpd'] = df_docs.groupby('parent_molregno')['year'].transform('min')
    df_docs = df_docs[['parent_molregno', 'first_publication_cpd']].drop_duplicates()

    df_combined = df_combined.merge(df_docs, on = 'parent_molregno', how = 'left')

    return df_combined


def add_chembl_properties_and_structures(df_combined: pd.DataFrame, chembl_con: sqlite3.Connection) -> (pd.DataFrame, pd.DataFrame):
    """
    Add compound properties from the compound_properties table (e.g., alogp, #hydrogen bond acceptors / donors, etc.).
    Add InChI, InChI key and canonical smiles. 

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: - Pandas DataFrame with added compound properties and structures. \\
        - Pandas DataFrame with compound properties and structures for all compound ids in ChEMBL.
    :rtype: (pd.DataFrame, pd.DataFrame)
    """
    sql = '''
    SELECT DISTINCT mh.parent_molregno, 
        cp.mw_freebase, cp.alogp, cp.hba, cp.hbd, cp.psa, cp.rtb, cp.ro3_pass, cp.num_ro5_violations, 
        cp.cx_most_apka, cp.cx_most_bpka, cp.cx_logp, cp.cx_logd, cp.molecular_species, cp.full_mwt, 
        cp.aromatic_rings, cp.heavy_atoms, cp.qed_weighted, cp.mw_monoisotopic, cp.full_molformula, 
        cp.hba_lipinski, cp.hbd_lipinski, cp.num_lipinski_ro5_violations, 
        struct.standard_inchi, struct.standard_inchi_key, struct.canonical_smiles
    FROM compound_properties cp
    INNER JOIN molecule_hierarchy mh
        ON cp.molregno = mh.parent_molregno
    INNER JOIN compound_structures struct
        ON mh.parent_molregno = struct.molregno
    '''

    df_cpd_props = pd.read_sql_query(sql, con=chembl_con)

    df_combined = df_combined.merge(df_cpd_props, on = 'parent_molregno', how = 'left')

    # # TODO: include?
    # ############### TESTING: compound props ###############
    # add_dataset_sizes(df_combined, "cpd props")

    return df_combined, df_cpd_props


def add_ligand_efficiency_metrics(df_combined: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the ligand efficiency metrics for the compounds 
    based on the mean pchembl values for a compound-target pair and the following ligand efficiency (LE) formulas: 

    .. math::
        LE &= \\frac{\\Delta G}{HA}
            \\qquad \\qquad \\text{where } \\Delta G = - RT \\ln(K_d)
            \\text{, } - RT\\ln(K_i)
            \\text{,  or} - RT\\ln(IC_{50}) 

        LE &= \\frac{2.303 \\cdot 298 \\cdot 0.00199 \\cdot pchembl \\_ value} {heavy \\_ atoms}

        BEI &= \\frac{pchembl \\_ mean \cdot 1000}{mw \\_ freebase}

        SEI &= \\frac{pchembl \\_ mean \cdot 100}{PSA}

        LLE &= pchembl \\_ mean - ALOGP
    
    Since LE metrics are based on pchembl values, they are calculated twice.
    Once for the pchembl values based on binding + functional assays (BF) 
    and once for the pchembl values based on binding assays only (B).
    
    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :return: Pandas DataFrame with added ligand efficiency metrics
    :rtype: pd.DataFrame
    """
    for suffix in ['BF', 'B']:
        df_combined['LE_'+suffix] = df_combined['pchembl_value_mean_'+suffix]/df_combined['heavy_atoms']*(2.303*298*0.00199)
        # replace infinity values with None as they are not useful
        df_combined['LE_'+suffix] = df_combined['LE_'+suffix].replace(np.inf, None)
        
        df_combined['BEI_'+suffix] = df_combined['pchembl_value_mean_'+suffix]*1000/df_combined["mw_freebase"]
        df_combined['BEI_'+suffix] = df_combined['BEI_'+suffix].replace(np.inf, None)
        
        df_combined['SEI_'+suffix] = df_combined['pchembl_value_mean_'+suffix]*100/df_combined["psa"]
        df_combined['SEI_'+suffix] = df_combined['SEI_'+suffix].replace(np.inf, None)
        
        df_combined['LLE_'+suffix] = df_combined['pchembl_value_mean_'+suffix]-df_combined["alogp"]
        
        df_combined = df_combined.astype({
        'LE_'+suffix: 'float64',
        'BEI_'+suffix: 'float64',
        'SEI_'+suffix: 'float64',
        'LLE_'+suffix: 'float64'
        })

    return df_combined


def add_atc_classification(df_combined: pd.DataFrame, chembl_con: sqlite3.Connection) -> (pd.DataFrame, pd.DataFrame):
    """
    Query and add ATC classifications (level 1) from the atc_classification and molecule_atc_classification tables.
    ATC level annotations for the same parent_molregno are combined into one description 
    that concatenates all descriptions sorted alphabetically into one string with ' | ' as a separator.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :return: - Pandas DataFrame with added ATC classifications \\
        - Pandas DataFrame with ATC annotations in ChEMBL
    :rtype: (pd.DataFrame, pd.DataFrame)
    """
    sql = '''
    SELECT DISTINCT mh.parent_molregno, atc.level1, level1_description
    FROM atc_classification atc
    INNER JOIN molecule_atc_classification matc
        ON atc.level5 = matc.level5
    INNER JOIN molecule_hierarchy mh
        ON matc.molregno = mh.molregno
    '''

    atc_levels = pd.read_sql_query(sql, con=chembl_con)
    atc_levels["l1_full"] = atc_levels["level1"] + "_" + atc_levels["level1_description"]

    # Combine ATC level annotations
    between_str_join = ' | '
    atc_levels['atc_level1'] = atc_levels.groupby(['parent_molregno'])['l1_full'].transform(lambda x: between_str_join.join(sorted(x)))
    atc_levels = atc_levels[['parent_molregno', 'atc_level1']].drop_duplicates()

    df_combined = df_combined.merge(atc_levels, on='parent_molregno', how = 'left')

    return df_combined, atc_levels



def add_all_chembl_compound_properties(df_combined: pd.DataFrame, chembl_con: sqlite3.Connection, limit_to_literature: bool) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    """
    Add ChEMBL-based compound properties to the given compound-target pairs, specifically:

    - the first publication date of a compound (first_publication_cpd)
    - ChEMBL compound properties
    - InChI, InChI key and canonical smiles
    - ligand efficiency metrics
    - ATC classifications

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Base first_publication_cpd on literature sources only if True. Base it on all available sources otherwise.
    :type limit_to_literature: bool
    :return: - Pandas DataFrame with added compound properties \\
        - Pandas DataFrame with compound properties and structures for all compound ids in ChEMBL \\
        - Pandas DataFrame with ATC annotations in ChEMBL
    :rtype: (pd.DataFrame, pd.DataFrame, pd.DataFrame)
    """
    df_combined = add_first_publication_date(df_combined, chembl_con, limit_to_literature)
    df_combined, df_cpd_props = add_chembl_properties_and_structures(df_combined, chembl_con)
    df_combined = add_ligand_efficiency_metrics(df_combined)
    df_combined, atc_levels = add_atc_classification(df_combined, chembl_con)
    return df_combined, df_cpd_props, atc_levels

