import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import PandasTools
from tqdm import tqdm

def add_built_in_descriptors(df_combined):
    # Built-in Compound Descriptors
    # Add relevant compound descriptors using built-in RDKit methods. 
    # Split table into two sections; PandasTools has difficulties working with null values
    df_combined_w_smiles = df_combined[df_combined['canonical_smiles'].notnull()].copy()
    df_combined_wo_smiles = df_combined[df_combined['canonical_smiles'].isnull()].copy()

    # add a column with RDKit molecules, used to calculate the descriptors
    PandasTools.AddMoleculeColumnToFrame(df_combined_w_smiles, 'canonical_smiles', 'mol', includeFingerprints=False)

    df_combined_w_smiles.loc[:,'fraction_csp3'] = df_combined_w_smiles['mol'].apply(Descriptors.FractionCSP3)
    df_combined_w_smiles.loc[:,'num_aliphatic_carbocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAliphaticCarbocycles)
    df_combined_w_smiles.loc[:,'num_aliphatic_heterocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAliphaticHeterocycles)
    df_combined_w_smiles.loc[:,'num_aliphatic_rings'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAliphaticRings)
    df_combined_w_smiles.loc[:,'num_aromatic_carbocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAromaticCarbocycles)
    df_combined_w_smiles.loc[:,'num_aromatic_heterocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAromaticHeterocycles)
    df_combined_w_smiles.loc[:,'num_aromatic_rings'] = df_combined_w_smiles['mol'].apply(Descriptors.NumAromaticRings)
    df_combined_w_smiles.loc[:,'num_heteroatoms'] = df_combined_w_smiles['mol'].apply(Descriptors.NumHeteroatoms)
    df_combined_w_smiles.loc[:,'num_saturated_carbocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumSaturatedCarbocycles)
    df_combined_w_smiles.loc[:,'num_saturated_heterocycles'] = df_combined_w_smiles['mol'].apply(Descriptors.NumSaturatedHeterocycles)
    df_combined_w_smiles.loc[:,'num_saturated_rings'] = df_combined_w_smiles['mol'].apply(Descriptors.NumSaturatedRings)
    df_combined_w_smiles.loc[:,'ring_count'] = df_combined_w_smiles['mol'].apply(Descriptors.RingCount)
    df_combined_w_smiles.loc[:,'num_stereocentres'] = df_combined_w_smiles['mol'].apply(Chem.rdMolDescriptors.CalcNumAtomStereoCenters)

    # add scaffolds
    PandasTools.AddMurckoToFrame(df_combined_w_smiles, 'mol', 'scaffold_w_stereo')
    # remove stereo information of the molecule to add scaffolds without stereo information
    df_combined_w_smiles['mol'].apply(Chem.RemoveStereochemistry)
    PandasTools.AddMurckoToFrame(df_combined_w_smiles, 'mol', 'scaffold_wo_stereo')

    # drop the column with RDKit molecules
    df_combined_w_smiles = df_combined_w_smiles.drop(['mol'] , axis=1)

    # combined both sections of the table
    df_combined = pd.concat([df_combined_w_smiles, 
                             df_combined_wo_smiles]).reset_index(drop=True)
    
    return df_combined, df_combined_w_smiles

def calculate_aromatic_atoms(smiles_set):
    aromatic_atoms_dict = dict()
    aromatic_c_dict = dict()
    aromatic_n_dict = dict()
    aromatic_hetero_dict = dict()
    
    for smiles in tqdm(smiles_set):
        mol = Chem.MolFromSmiles(smiles)
        aromatic_atoms_dict[smiles] = sum(mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms()))
        aromatic_c_dict[smiles] = sum((mol.GetAtomWithIdx(i).GetIsAromatic() & (mol.GetAtomWithIdx(i).GetAtomicNum() == 6)) for i in range(mol.GetNumAtoms()))
        aromatic_n_dict[smiles] = sum((mol.GetAtomWithIdx(i).GetIsAromatic() & (mol.GetAtomWithIdx(i).GetAtomicNum() == 7)) for i in range(mol.GetNumAtoms()))
        aromatic_hetero_dict[smiles] = sum((mol.GetAtomWithIdx(i).GetIsAromatic() & (mol.GetAtomWithIdx(i).GetAtomicNum() != 6) & (mol.GetAtomWithIdx(i).GetAtomicNum() != 1)) for i in range(mol.GetNumAtoms()))
        
    return aromatic_atoms_dict, aromatic_c_dict, aromatic_n_dict, aromatic_hetero_dict

def add_aromaticity_descriptors(df_combined, df_combined_w_smiles):
    # use df_combined_w_smiles to exclude null values
    smiles_set = set(df_combined_w_smiles["canonical_smiles"])
    aromatic_atoms_dict, aromatic_c_dict, aromatic_n_dict, aromatic_hetero_dict = calculate_aromatic_atoms(smiles_set)

    df_combined['aromatic_atoms'] = df_combined['canonical_smiles'].map(aromatic_atoms_dict)
    df_combined['aromatic_c'] = df_combined['canonical_smiles'].map(aromatic_c_dict)
    df_combined['aromatic_n'] = df_combined['canonical_smiles'].map(aromatic_n_dict)
    df_combined['aromatic_hetero'] = df_combined['canonical_smiles'].map(aromatic_hetero_dict)

    return df_combined


def add_rdkit_compound_descriptors(df_combined):
    df_combined, df_combined_w_smiles = add_built_in_descriptors(df_combined)
    df_combined = add_aromaticity_descriptors(df_combined, df_combined_w_smiles)

    return df_combined



    