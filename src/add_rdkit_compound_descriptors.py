"""
Add RDKit-based compound properties to the dataset.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from tqdm import tqdm

from dataset import Dataset
import sanity_checks


def add_built_in_descriptors(dataset: Dataset):
    """
    Add RDKit built-in compound descriptors.

    :param dataset: Dataset with compound-target pairs.
        Will be updated to only include built-in RDKit compound descriptors.
    :type dataset: Dataset
    :param df_combined: Pandas DataFrame with compound-target pairs
    :type df_combined: pd.DataFrame
    """
    # add a column with RDKit molecules, used to calculate the descriptors
    PandasTools.AddMoleculeColumnToFrame(
        dataset.df_result, "canonical_smiles", "mol", includeFingerprints=False
    )

    dataset.df_result.loc[:, "fraction_csp3"] = dataset.df_result["mol"].apply(
        Descriptors.FractionCSP3
    )
    dataset.df_result.loc[:, "ring_count"] = dataset.df_result["mol"].apply(
        Descriptors.RingCount
    )
    dataset.df_result.loc[:, "num_aliphatic_rings"] = dataset.df_result["mol"].apply(
        Descriptors.NumAliphaticRings
    )
    dataset.df_result.loc[:, "num_aliphatic_carbocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumAliphaticCarbocycles)
    dataset.df_result.loc[:, "num_aliphatic_heterocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumAliphaticHeterocycles)
    dataset.df_result.loc[:, "num_aromatic_rings"] = dataset.df_result["mol"].apply(
        Descriptors.NumAromaticRings
    )
    dataset.df_result.loc[:, "num_aromatic_carbocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumAromaticCarbocycles)
    dataset.df_result.loc[:, "num_aromatic_heterocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumAromaticHeterocycles)
    dataset.df_result.loc[:, "num_saturated_rings"] = dataset.df_result["mol"].apply(
        Descriptors.NumSaturatedRings
    )
    dataset.df_result.loc[:, "num_saturated_carbocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumSaturatedCarbocycles)
    dataset.df_result.loc[:, "num_saturated_heterocycles"] = dataset.df_result[
        "mol"
    ].apply(Descriptors.NumSaturatedHeterocycles)
    dataset.df_result.loc[:, "num_stereocentres"] = dataset.df_result["mol"].apply(
        Chem.rdMolDescriptors.CalcNumAtomStereoCenters
    )
    dataset.df_result.loc[:, "num_heteroatoms"] = dataset.df_result["mol"].apply(
        Descriptors.NumHeteroatoms
    )

    # add scaffolds
    PandasTools.AddMurckoToFrame(dataset.df_result, "mol", "scaffold_w_stereo")
    # remove stereo information of the molecule to add scaffolds without stereo information
    dataset.df_result["mol"].apply(Chem.RemoveStereochemistry)
    PandasTools.AddMurckoToFrame(dataset.df_result, "mol", "scaffold_wo_stereo")

    # drop the column with RDKit molecules
    dataset.df_result = dataset.df_result.drop(["mol"], axis=1)


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


def add_aromaticity_descriptors(dataset: Dataset):
    """
    Add number of aromatic atoms in a compounds, specifically:

    - total # aromatics atoms (aromatic_atoms)
    - # aromatic carbon atoms (aromatic_c)
    - # aromatic nitrogen atoms (aromatic_n)
    - # aromatic hetero atoms (aromatic_hetero)

    :param dataset: Dataset with compound-target pairs.
        Will be updated to only include counts of aromatic atoms
    :type dataset: Dataset
    """
    # use df_combined_w_smiles to exclude null values
    smiles_set = set(dataset.df_result["canonical_smiles"])
    aromatic_atoms_dict, aromatic_c_dict, aromatic_n_dict, aromatic_hetero_dict = (
        calculate_aromatic_atoms(smiles_set)
    )

    dataset.df_result["aromatic_atoms"] = dataset.df_result["canonical_smiles"].map(
        aromatic_atoms_dict
    )
    dataset.df_result["aromatic_c"] = dataset.df_result["canonical_smiles"].map(
        aromatic_c_dict
    )
    dataset.df_result["aromatic_n"] = dataset.df_result["canonical_smiles"].map(
        aromatic_n_dict
    )
    dataset.df_result["aromatic_hetero"] = dataset.df_result["canonical_smiles"].map(
        aromatic_hetero_dict
    )


def add_rdkit_compound_descriptors(dataset: Dataset):
    """
    Add RDKit-based compound descriptors (built-in and numbers of aromatic atoms).

    :param dataset: Dataset with compound-target pairs.
        Will be updated to only include
        built-in RDKit compound descriptors
        and numbers of aromatic atoms.
    :type dataset: Dataset
    """
    add_built_in_descriptors(dataset)
    add_aromaticity_descriptors(dataset)
    sanity_checks.check_rdkit_props(dataset.df_result)
