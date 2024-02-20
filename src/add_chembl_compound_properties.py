import sqlite3

import pandas as pd

from dataset import Dataset


########### Add Compound Properties Based on ChEMBL Data ###########
def add_first_publication_date(
    dataset: Dataset, chembl_con: sqlite3.Connection, limit_to_literature: bool
):
    """
    Query and calculate the first publication of a compound
    based on ChEMBL data (column name: first_publication_cpd).
    If limit_to_literature is True, this corresponds to the first appearance
    of the compound in the literature according to ChEMBL.
    Otherwise this is the first appearance in any source in ChEMBL.

    :param dataset: Dataset with compound-target pairs.
        Will be updated to include first_publication_cpd
    :type dataset: Dataset
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Base first_publication_cpd on literature sources only if True.
    :type limit_to_literature: bool
    """
    # information about salts is aggregated in the parent
    sql = """
    SELECT DISTINCT docs.year, mh.parent_molregno
    FROM docs
    LEFT JOIN compound_records cr
        ON docs.doc_id = cr.doc_id
    INNER JOIN molecule_hierarchy mh 
        ON cr.molregno = mh.molregno   -- cr.molregno = salt_molregno
    WHERE docs.year is not null
    """
    if limit_to_literature:
        sql += """    and docs.src_id = 1"""
    df_docs = pd.read_sql_query(sql, con=chembl_con)

    df_docs["first_publication_cpd"] = df_docs.groupby("parent_molregno")[
        "year"
    ].transform("min")
    df_docs = df_docs[["parent_molregno", "first_publication_cpd"]].drop_duplicates()

    dataset.df_result = dataset.df_result.merge(
        df_docs, on="parent_molregno", how="left"
    )


def add_chembl_properties_and_structures(
    dataset: Dataset, chembl_con: sqlite3.Connection
):
    """
    Add compound properties from the compound_properties table
    (e.g., alogp, #hydrogen bond acceptors / donors, etc.).
    Add InChI, InChI key and canonical smiles.

    :param dataset: Dataset with compound-target pairs.
        Will be updated to include compound properties and structures.
        dataset.df_cpd_props will be set to
        compound properties and structures for all compound ids in ChEMBL.
    :type dataset: Dataset
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    """
    sql = """
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
    """

    df_cpd_props = pd.read_sql_query(sql, con=chembl_con)
    dataset.df_cpd_props = df_cpd_props

    dataset.df_result = dataset.df_result.merge(
        df_cpd_props, on="parent_molregno", how="left"
    )


def add_ligand_efficiency_metrics(dataset: Dataset):
    """
    Calculate the ligand efficiency metrics for the compounds
    based on the mean pchembl values for a compound-target pair and
    the following ligand efficiency (LE) formulas:

    .. math::
        LE &= \\frac{\\Delta G}{HA}
            \\qquad \\qquad \\text{where } \\Delta G = - RT \\ln(K_d)
            \\text{, } - RT\\ln(K_i)
            \\text{,  or} - RT\\ln(IC_{50})

        LE &= \\frac{2.303 \\cdot 298 \\cdot 0.00199 \\cdot pchembl \\_ value} {heavy \\_ atoms}

        BEI &= \\frac{pchembl \\_ mean \\cdot 1000}{mw \\_ freebase}

        SEI &= \\frac{pchembl \\_ mean \\cdot 100}{PSA}

        LLE &= pchembl \\_ mean - ALOGP

    Since LE metrics are based on pchembl values, they are calculated twice.
    Once for the pchembl values based on binding + functional assays (BF)
    and once for the pchembl values based on binding assays only (B).

    :param dataset: Dataset with compound-target pairs.
        Will be updated to include ligand efficiency metrics.
    :type dataset: Dataset
    """
    for suffix in ["BF", "B"]:
        dataset.df_result.loc[dataset.df_result["heavy_atoms"] != 0, f"LE_{suffix}"] = (
            dataset.df_result[f"pchembl_value_mean_{suffix}"]
            / dataset.df_result["heavy_atoms"]
            * (2.303 * 298 * 0.00199)
        )

        dataset.df_result.loc[
            dataset.df_result["mw_freebase"] != 0, f"BEI_{suffix}"
        ] = (
            dataset.df_result[f"pchembl_value_mean_{suffix}"]
            * 1000
            / dataset.df_result["mw_freebase"]
        )

        dataset.df_result.loc[dataset.df_result["psa"] != 0, f"SEI_{suffix}"] = (
            dataset.df_result[f"pchembl_value_mean_{suffix}"]
            * 100
            / dataset.df_result["psa"]
        )

        dataset.df_result[f"LLE_{suffix}"] = (
            dataset.df_result[f"pchembl_value_mean_{suffix}"]
            - dataset.df_result["alogp"]
        )

        dataset.df_result = dataset.df_result.astype(
            {
                f"LE_{suffix}": "float64",
                f"BEI_{suffix}": "float64",
                f"SEI_{suffix}": "float64",
                f"LLE_{suffix}": "float64",
            }
        )


def add_atc_classification(dataset: Dataset, chembl_con: sqlite3.Connection):
    """
    Query and add ATC classifications (level 1) from the atc_classification and
    molecule_atc_classification tables.
    ATC level annotations for the same parent_molregno are combined into one description
    that concatenates all descriptions sorted alphabetically
    into one string with ' | ' as a separator.

    :param dataset: Dataset with compound-target pairs.
        Will be updated to include ATC classifications.
        dataset.atc_levels will be set to ATC annotations in ChEMBL.
    :type dataset: Dataset
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    """
    sql = """
    SELECT DISTINCT mh.parent_molregno, atc.level1, atc.level1_description
    FROM atc_classification atc
    INNER JOIN molecule_atc_classification matc
        ON atc.level5 = matc.level5
    INNER JOIN molecule_hierarchy mh
        ON matc.molregno = mh.molregno
    """

    atc_levels = pd.read_sql_query(sql, con=chembl_con)
    atc_levels["l1_full"] = (
        atc_levels["level1"] + "_" + atc_levels["level1_description"]
    )

    # Combine ATC level annotations
    between_str_join = " | "
    atc_levels["atc_level1"] = atc_levels.groupby(["parent_molregno"])[
        "l1_full"
    ].transform(lambda x: between_str_join.join(sorted(x)))
    atc_levels = atc_levels[["parent_molregno", "atc_level1"]].drop_duplicates()

    dataset.atc_levels = atc_levels

    dataset.df_result = dataset.df_result.merge(
        atc_levels, on="parent_molregno", how="left"
    )


def add_all_chembl_compound_properties(
    dataset: Dataset, chembl_con: sqlite3.Connection, limit_to_literature: bool
):
    """
    Add ChEMBL-based compound properties to the given compound-target pairs, specifically:

    - the first publication date of a compound (first_publication_cpd)
    - ChEMBL compound properties
    - InChI, InChI key and canonical smiles
    - ligand efficiency metrics
    - ATC classifications

    :param dataset: Dataset with compound-target pairs.
        Will be updated to include compound properties.
    :type dataset: Dataset
    :param chembl_con: Sqlite3 connection to ChEMBL database.
    :type chembl_con: sqlite3.Connection
    :param limit_to_literature: Base first_publication_cpd on literature sources only if True.
        Base it on all available sources otherwise.
    :type limit_to_literature: bool
    """
    add_first_publication_date(dataset, chembl_con, limit_to_literature)

    add_chembl_properties_and_structures(dataset, chembl_con)

    add_ligand_efficiency_metrics(dataset)

    add_atc_classification(dataset, chembl_con)
