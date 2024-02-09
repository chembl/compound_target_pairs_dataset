Introduction
============
This code extract a dataset of compound-target pairs from the open-source bioactivity database `ChEMBL`_ [Zdrazil2023]_. 

The compound-target pairs are known to interact because 

- they have at least one corresponding measured activity values in ChEMBL or 
- they are part of a set of manually curated known interactions in ChEMBL.

Furthermore, the dataset contains a number of compounds and target annotations to enable future analyses. 

Previously, a similar dataset has been curated manually and has been used to investigate target-based differences in drug-like properties and ligand efficiencies [Leeson2021]_. 
This code can generate an extended version of the previous dataset for every ChEMBL version from ChEMBL 26 onwards.  

.. _ChEMBL: https://www.ebi.ac.uk/chembl/

.. [Zdrazil2023] Zdrazil et al., "The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods",
    Nucleic Acids Research, gkad1004, 2023, https://doi.org/10.1093/nar/gkad1004
.. [Leeson2021] Leeson et al., "Target-Based Evaluation of “Drug-Like” Properties and Ligand Efficiencies", 
    Journal of Medicinal Chemistry, 64(11), 7210-7230, 2021, https://doi.org/10.1021/acs.jmedchem.1c00416


Dataset Documentation
*********************
If you are interested in understanding the fields in the resulting dataset, see :doc:`columns_docs`

User Guide
**********
If you are interested in using the code, see :doc:`user_guide`

Code Documentation
******************
If you are interested in understanding the code, see :doc:`modules`

