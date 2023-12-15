Columns Documentation
=====================

Initial Query
*************
.. csv-table:: 
   :file: tables/InitialQuery.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1

Aggregated Values
*****************
Aggregated per compound-target pair using parent_molregno and tid_mutation:

- _BF: based on binding + functional assays

- _B: based on binding assays

.. csv-table:: 
   :file: tables/AggregatedValues.csv
   :widths: 10, 10, 10, 10
   :header-rows: 1

DTI (Drug-Target Interaction) Annotations
*****************************************
.. csv-table:: 
   :file: tables/DTIAnnotations.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1

Mechanism to assign DTI
***********************
.. csv-table:: 
   :file: tables/DTIMechanism.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1

- Is the compound-target pair in the drug_mechanisms table?	= Is it a known relevant compound-target interaction?

- What is the max_phase of the compound? 				        = Is it a drug / clinical compound?

- Is the target in the drug_mechanisms table? 			    = Is it a therapeutic target?


There are three possible annotations in ChEMBL with max_phase not between 1 and 4:

- 0.5 = early phase 1 clinical trials

- -1 = clinical phase unknown for drug or clinical candidate drug, i.e., where ChEMBL cannot assign a clinical phase

- NULL = preclinical compounds with bioactivity data

All three are grouped together into the annotation C0_DT.



Compound and Target Properties Based on ChEMBL Data
***************************************************
.. csv-table:: 
   :file: tables/CompoundProps.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1


Ligand Efficiency Metrics 
*************************
.. csv-table:: 
   :file: tables/LEMetrics.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1

- _BF: based on pchembl_value_mean_BF (based on binding + functional assays)

- _B: based on pchembl_value_mean_B (based on binding assays)


RDKit-Based Compound Descriptors
********************************
.. csv-table:: 
   :file: tables/RDKitProps.csv
   :widths: 10, 10, 10, 10, 10
   :header-rows: 1


Annotations for Filtering
*************************
available for full dataset to facilitate filtering into subsets

.. csv-table:: 
   :file: tables/Filtering.csv
   :widths: 10, 10, 10, 10
   :header-rows: 1

* comparator compounds must have a pchembl value but don't have a required max_phase, i.e., drugs and clinical candidates are counted as comparator compounds

