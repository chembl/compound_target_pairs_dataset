import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from tqdm import tqdm


def add_built_in_descriptors(df_combined: pd.DataFrame) -> pd.DataFrame:
    """
    Add RDKit built-in compound descriptors.

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :return: Pandas DataFrame with added built-in RDKit compound descriptors
    :rtype: pd.DataFrame
    """
    # add a column with RDKit molecules, used to calculate the descriptors
    PandasTools.AddMoleculeColumnToFrame(
        df_combined, "canonical_smiles", "mol", includeFingerprints=False
    )

    df_combined.loc[:, "fraction_csp3"] = df_combined["mol"].apply(
        Descriptors.FractionCSP3
    )
    df_combined.loc[:, "ring_count"] = df_combined["mol"].apply(Descriptors.RingCount)
    df_combined.loc[:, "num_aliphatic_rings"] = df_combined["mol"].apply(
        Descriptors.NumAliphaticRings
    )
    df_combined.loc[:, "num_aliphatic_carbocycles"] = df_combined["mol"].apply(
        Descriptors.NumAliphaticCarbocycles
    )
    df_combined.loc[:, "num_aliphatic_heterocycles"] = df_combined["mol"].apply(
        Descriptors.NumAliphaticHeterocycles
    )
    df_combined.loc[:, "num_aromatic_rings"] = df_combined["mol"].apply(
        Descriptors.NumAromaticRings
    )
    df_combined.loc[:, "num_aromatic_carbocycles"] = df_combined["mol"].apply(
        Descriptors.NumAromaticCarbocycles
    )
    df_combined.loc[:, "num_aromatic_heterocycles"] = df_combined["mol"].apply(
        Descriptors.NumAromaticHeterocycles
    )
    df_combined.loc[:, "num_saturated_rings"] = df_combined["mol"].apply(
        Descriptors.NumSaturatedRings
    )
    df_combined.loc[:, "num_saturated_carbocycles"] = df_combined["mol"].apply(
        Descriptors.NumSaturatedCarbocycles
    )
    df_combined.loc[:, "num_saturated_heterocycles"] = df_combined["mol"].apply(
        Descriptors.NumSaturatedHeterocycles
    )
    df_combined.loc[:, "num_stereocentres"] = df_combined["mol"].apply(
        Chem.rdMolDescriptors.CalcNumAtomStereoCenters
    )
    df_combined.loc[:, "num_heteroatoms"] = df_combined["mol"].apply(
        Descriptors.NumHeteroatoms
    )

    # add scaffolds
    PandasTools.AddMurckoToFrame(df_combined, "mol", "scaffold_w_stereo")
    # remove stereo information of the molecule to add scaffolds without stereo information
    df_combined["mol"].apply(Chem.RemoveStereochemistry)
    PandasTools.AddMurckoToFrame(df_combined, "mol", "scaffold_wo_stereo")

    # drop the column with RDKit molecules
    df_combined = df_combined.drop(["mol"], axis=1)

    return df_combined


def calculate_aromatic_atoms(
    smiles_set: set[str],
) -> tuple[dict[str, int], dict[str, int], dict[str, int], dict[str, int]]:
    """
    Get dictionaries with number of aromatic atoms for each smiles.

    :param smiles_set: Set of smiles to calculate the number of aromatic atoms for
    :type smiles_set: set[str]
    :return: Dictionaries with:

        - SMILES -> # aromatics atoms
        - SMILES -> # aromatic carbon atoms
        - SMILES -> # aromatic nitrogen atoms
        - SMILES -> # aromatic hetero atoms
    :rtype: (dict[str, int], dict[str, int], dict[str, int], dict[str, int])
    """
    aromatic_atoms_dict = {}
    aromatic_c_dict = {}
    aromatic_n_dict = {}
    aromatic_hetero_dict = {}

    for smiles in tqdm(smiles_set):
        mol = Chem.MolFromSmiles(smiles)
        aromatic_atoms_dict[smiles] = sum(
            mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())
        )
        aromatic_c_dict[smiles] = sum(
            (
                mol.GetAtomWithIdx(i).GetIsAromatic()
                & (mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            )
            for i in range(mol.GetNumAtoms())
        )
        aromatic_n_dict[smiles] = sum(
            (
                mol.GetAtomWithIdx(i).GetIsAromatic()
                & (mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
            )
            for i in range(mol.GetNumAtoms())
        )
        aromatic_hetero_dict[smiles] = sum(
            (
                mol.GetAtomWithIdx(i).GetIsAromatic()
                & (mol.GetAtomWithIdx(i).GetAtomicNum() != 6)
                & (mol.GetAtomWithIdx(i).GetAtomicNum() != 1)
            )
            for i in range(mol.GetNumAtoms())
        )

    return aromatic_atoms_dict, aromatic_c_dict, aromatic_n_dict, aromatic_hetero_dict


def add_aromaticity_descriptors(df_combined: pd.DataFrame) -> pd.DataFrame:
    """
    Add number of aromatic atoms in a compounds, specifically:

    - total # aromatics atoms (aromatic_atoms)
    - # aromatic carbon atoms (aromatic_c)
    - # aromatic nitrogen atoms (aromatic_n)
    - # aromatic hetero atoms (aromatic_hetero)

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :return: Pandas DataFrame with added counts of aromatic atoms
    :rtype: pd.DataFrame
    """
    # use df_combined_w_smiles to exclude null values
    smiles_set = set(df_combined["canonical_smiles"])
    aromatic_atoms_dict, aromatic_c_dict, aromatic_n_dict, aromatic_hetero_dict = (
        calculate_aromatic_atoms(smiles_set)
    )

    df_combined["aromatic_atoms"] = df_combined["canonical_smiles"].map(
        aromatic_atoms_dict
    )
    df_combined["aromatic_c"] = df_combined["canonical_smiles"].map(aromatic_c_dict)
    df_combined["aromatic_n"] = df_combined["canonical_smiles"].map(aromatic_n_dict)
    df_combined["aromatic_hetero"] = df_combined["canonical_smiles"].map(
        aromatic_hetero_dict
    )

    return df_combined


def add_rdkit_compound_descriptors(df_combined: pd.DataFrame) -> pd.DataFrame:
    """
    Add RDKit-based compound descriptors (built-in and numbers of aromatic atoms).

    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    :return: Pandas DataFrame with added built-in RDKit compound descriptors
        and numbers of aromatic atoms
    :rtype: pd.DataFrame
    """
    df_combined = add_built_in_descriptors(df_combined)
    df_combined = add_aromaticity_descriptors(df_combined)

    return df_combined
