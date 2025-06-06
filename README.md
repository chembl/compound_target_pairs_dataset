<h1 align="center">
    Dataset of Interacting Compound-Target Pairs in ChEMBL
</h1>

<p align="center">
    <a href="https://zenodo.org/doi/10.5281/zenodo.10723114"><img src="https://zenodo.org/badge/550876229.svg" alt="Code DOI"></a>
    <a href="https://doi.org/10.5281/zenodo.10721939"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10721939.svg" alt="Dataset DOI"></a>
</p>

# Introduction
This code extracts a dataset of compound-target pairs from the open-source bioactivity database [ChEMBL](https://www.ebi.ac.uk/chembl/) [Zdrazil2023]. 

The compound-target pairs are known to interact because 

- they have at least one corresponding measured activity value in ChEMBL or 
- they are part of a set of manually curated known interactions in ChEMBL.

Furthermore, the dataset contains a number of compound and target annotations to enable future analyses. 

Previously, a similar dataset has been curated manually and has been used to investigate target-based differences in drug-like properties and ligand efficiencies [Leeson2021]. 
This code can generate an extended version of the previous dataset for every ChEMBL version from ChEMBL 26 onwards.  

[Zdrazil2023]: Zdrazil et al., "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods",
    Nucleic Acids Research, gkad1004, 2023, https://doi.org/10.1093/nar/gkad1004

[Leeson2021]: Leeson et al., "Target-Based Evaluation of “Drug-Like” Properties and Ligand Efficiencies", 
    Journal of Medicinal Chemistry, 64(11), 7210-7230, 2021, https://doi.org/10.1021/acs.jmedchem.1c00416

# Dataset
The dataset for different ChEMBL versions from ChEMBL 26 onwards is available [here](https://ftp.ebi.ac.uk/pub/databases/chembl/Drug_Target_dataset/).

# Quick Start
## Dependencies
Install the required dependencies with
```
pip install .
```

Note: Using Pandas version 2.2 will lead to warnings regarding the RDKit PandasTools when running the code. 
However, the final dataset is not impacted. 


## Generating the Dataset
The default version of the dataset (the full dataset as a CSV file based on the newest ChEMBL version) can be generated by calling 
```
python main.py -o <output_path>
```

An overview of the available arguments to modify the output is available by calling 

```
python main.py --help
```

## Documentation
The full documentation is available [here](https://chembl.github.io/compound_target_pairs_dataset/).

The corresponding article is available [here](https://doi.org/10.1038/s41597-024-03582-9).
