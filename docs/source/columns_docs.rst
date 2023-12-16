Columns in the Final Dataset
============================
This page provides explanations for all columns available in the final dataset. 

More information on ChEMBL-based columns can be found in the respective `ChEMBL schema documentation`_. 
The information on this page mostly corresponds to the `ChEMBL 32 schema documentation`_.

.. _ChEMBL 32 schema documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_32/schema_documentation.html
.. _ChEMBL schema documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/


Initial Query
*************
.. csv-table:: 
   :file: tables/InitialQuery.csv
   :widths: 10, 10, 10, 10, 30
   :header-rows: 1

.. [#] There have been changes to the max_phase field in ChEMBL with `version 32`_. See `Maximum Phase in ChEMBL`_.

.. _version 32: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_32/chembl_32_release_notes.txt



Aggregated Values
*****************
Aggregated per compound-target pair using parent_molregno and tid_mutation.

.. csv-table:: 
   :file: tables/AggregatedValues.csv
   :widths: 20, 10, 10, 60
   :header-rows: 1

Naming Convention: B vs. BF
---------------------------
These values are aggregated based on different subsets of the full dataset. 
The corresponding columns in the final dataset have a suffix that corresponds to the assay types the value is based on: 

- _BF: based on binding + functional assays
- _B: based on binding assays



DTI (Drug-Target Interaction) Annotations
*****************************************
Based on cpd_target_pair, does not include mutation information.

.. csv-table:: 
   :file: tables/DTIAnnotations.csv
   :widths: 10, 10, 10, 10, 40
   :header-rows: 1


Mechanism to Assign DTI
-----------------------
.. csv-table:: 
   :file: tables/DTIMechanism.csv
   :widths: 10, 10, 10, 10, 60
   :header-rows: 1

.. [#] Is the compound-target pair in the drug_mechanisms table? = Is it a known relevant compound-target interaction?

.. [#] What is the max_phase of the compound? = Is it a drug / clinical compound?

.. [#] Is the target in the drug_mechanisms table? = Is it a therapeutic target?


.. [#] There have been changes to the max_phase field in ChEMBL with `version 32`_. 
   C0_DT groups together all compounds with a max_phase not between 1 and 4. See `Maximum Phase in ChEMBL`_


Maximum Phase in ChEMBL
-----------------------
Before ChEMBL 32, compounds with a max_phase not between 1 and 4 were assigned a max_phase of 0. 

| From ChEMBL 32 onwards, compounds with a max_phase not between 1 and 4 can have three possible values: 
|   - 0.5 = early phase 1 clinical trials
|   - -1 = clinical phase unknown for drug or clinical candidate drug, i.e., where ChEMBL cannot assign a clinical phase
|   - NULL = preclinical compounds with bioactivity data



Compound and Target Properties Based on ChEMBL Data
***************************************************
.. csv-table:: 
   :file: tables/CompoundProps.csv
   :widths: 10, 10, 10, 10, 60
   :header-rows: 1


Ligand Efficiency Metrics 
*************************
Calculated based on pchembl_value_mean. 

Since LE metrics are based on pchembl values, they are calculated twice.
Once for the pchembl values based on binding + functional assays (suffix _BF) 
and once for the pchembl values based on binding assays only (suffix _B).

.. csv-table:: 
   :file: tables/LEMetrics.csv
   :widths: 10, 10, 10, 20
   :header-rows: 1

.. [#] .. math:: LE &= \frac{\Delta G}{HA} \text{where } \Delta G = - RT \ln(K_d) \text{, } - RT\ln(K_i) \text{,  or} - RT\ln(IC_{50}) \\ LE &= \frac{2.303 \cdot 298 \cdot 0.00199 \cdot \text{pchembl_value}} {\text{heavy_atoms}}

.. [#] .. math:: BEI  = \frac{pchembl\text{_}mean \cdot 1000}{mw\text{_}freebase}

.. [#] .. math:: SEI = \frac{pchembl\text{_}mean \cdot 100}{PSA}

.. [#] .. math:: LLE = pchembl\text{_}mean - ALogP



RDKit-Based Compound Descriptors
********************************
Most RDKit-based compound descriptors are calculated using built-in RDKit methods from `Descriptors`_ and `rdMolDescriptors`_.
A few are calculated with bespoke RDKit-based methods. 

.. _Descriptors: https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
.. _rdMolDescriptors: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html

.. csv-table:: 
   :file: tables/RDKitProps.csv
   :widths: 10, 10, 10, 10, 60
   :header-rows: 1



Annotations for Filtering
*************************
Columns are only available for the full dataset to facilitate the filtering into subsets.


Helper Columns
--------------
.. csv-table:: 
   :file: tables/Filtering.csv
   :widths: 10, 10, 10, 70
   :header-rows: 1


Filtering Columns
-----------------
.. csv-table:: 
   :file: tables/Filtering2.csv
   :widths: 10, 10, 10, 10, 10, 70
   :header-rows: 1

.. [#] Comparator compounds in this context are all compounds with a pchembl_value_mean_BF / _B. 
   I.e., this includes compounds with a DTI of D_DT or C<p>_DT. 


