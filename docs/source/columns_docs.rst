Columns in the Final Dataset
============================
This page provides explanations for all columns available in the final dataset. 

More information on ChEMBL-based columns can be found in the respective `ChEMBL schema documentation`_. 
The information on this page mostly corresponds to the `ChEMBL 32 schema documentation`_.

.. _ChEMBL 32 schema documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_32/schema_documentation.html
.. _ChEMBL schema documentation: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/


Initial Query
*************
Pchembl Values
---------------
The pchembl_value is later aggregated into mean, max and median per compound-target pair and dropped.

.. csv-table:: 
   :file: tables/InitialQuery1.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

Compound Information
--------------------
.. csv-table:: 
   :file: tables/InitialQuery2.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

.. [#] There have been changes to the max_phase field in ChEMBL with `version 32`_. See `Maximum Phase in ChEMBL`_.

.. _version 32: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_32/chembl_32_release_notes.txt


Target Information
------------------
.. csv-table:: 
   :file: tables/InitialQuery3.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

Helper Columns
--------------
These columns are combination of other columns, used for easier processing of the dataset. 

.. csv-table:: 
   :file: tables/InitialQuery4.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1



Aggregated Values
*****************
Aggregated per compound-target pair using parent_molregno and tid_mutation.

.. csv-table:: 
   :file: tables/AggregatedValues.csv
   :widths: 30, 10, 15, 45
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
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1


Mechanism to Assign DTI
-----------------------
.. csv-table:: 
   :file: tables/DTIMechanism.csv
   :widths: 15, 15, 15, 10, 45
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
First publication
-----------------
In contrast to the aggregated time-related fields, 
this field takes all of ChEMBL and not just the time-related data within the dataset into account. 

.. csv-table:: 
   :file: tables/CompoundProps1.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

Compound Properties
-------------------
.. csv-table:: 
   :file: tables/CompoundProps2.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

Compound Structures
-------------------
.. csv-table:: 
   :file: tables/CompoundProps3.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

ATC and Target Class
--------------------
.. csv-table:: 
   :file: tables/CompoundProps4.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1


Ligand Efficiency Metrics 
*************************
Calculated based on pchembl_value_mean. 

Since LE metrics are based on pchembl values, they are calculated twice.
Once for the pchembl values based on binding and functional assays (suffix _BF) 
and once for the pchembl values based on binding assays only (suffix _B).

.. csv-table:: 
   :file: tables/LEMetrics.csv
   :widths: 20, 20, 20, 40
   :header-rows: 1


Equations
---------
.. math::
   :nowrap:

   \begin{flalign*}
   LE &= \frac{2.303 \cdot 298 \cdot 0.00199 \cdot pchembl\_value} {heavy\_atoms} \\
   BEI  &= \frac{pchembl\_mean \cdot 1000}{mw\_freebase} \\
   SEI &= \frac{pchembl\_mean \cdot 100}{PSA} \\
   LLE &= pchembl\_mean - ALogP \\
   \end{flalign*}



RDKit-Based Compound Descriptors
********************************
Built-in Methods
----------------
These compound descriptors are calculated using built-in RDKit methods from `Descriptors`_ and `rdMolDescriptors`_.

.. _Descriptors: https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html
.. _rdMolDescriptors: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html

.. csv-table:: 
   :file: tables/RDKitProps1.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1

Bespoke Methods
---------------
These compound descriptors are calculated using custom RDKit-based methods.

.. csv-table:: 
   :file: tables/RDKitProps2.csv
   :widths: 20, 10, 15, 20, 35
   :header-rows: 1



Annotations for Filtering
*************************
Columns are only available for the full dataset to facilitate the filtering into subsets.


Helper Columns
--------------
pair_mutation_in_dm_table and pair_in_dm_table are similar fields. 
They differ in whether mutation information is taken into account, 
reflecting that mutation information is only sometimes taken into account 
when calculating fields and adding rows to the dataset. 

   - pair_mutation_in_dm_table: 
      Is the compound-target pair in the drug_mechanism table 
      when taking mutation information into account? 
      Mutation information IS taken into account when adding pairs to the dataset 
      because they appear in the drug_mechanism table. 
      (cpd A, target B without mutation) will be added to the set of existing 
      compound-target pairs with pchembl values 
      if there is a pair with a pchembl value for (cpd A, target B with mutation C) 
      but there is no pair with a pchembl value for (cpd A, target B without mutation).
      It is used to determine keep_for_binding which in turn is used 
      to determine the B subset of data based on binding assays. 
   - pair_in_dm_table: 
      Is the compound-target pair in the drug_mechanism table
      when ignoring mutation information? 
      Mutation information is NOT taken into account when assigning DTI values. 

.. csv-table:: 
   :file: tables/Filtering.csv
   :widths: 20, 10, 15, 55
   :header-rows: 1


Filtering Columns
-----------------
.. csv-table:: 
   :file: tables/Filtering2.csv
   :widths: 20, 10, 15, 15, 15, 25
   :header-rows: 1

.. [#] Comparator compounds in this context are all compounds with a pchembl_value_mean_BF / _B. 
   I.e., this includes compounds with a DTI of D_DT or C<p>_DT. 


