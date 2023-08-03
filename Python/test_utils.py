def calculate_dataset_sizes(data):
    now_mols = len(set(data["parent_molregno"]))
    now_targets = len(set(data["tid"]))
    now_targets_mutation = len(set(data["tid_mutation"]))
    now_pairs = len(set(data['cpd_target_pair']))
    now_pairs_mutation = len(set(data['cpd_target_pair_mutation']))
    
    if 'DTI' in data.columns:
        # drugs = compounds of a compound-target pair with a known interaction  
        data_drugs = data[data["DTI"] == "D_DT"]
    else: 
        data_drugs = data[data["max_phase"] == 4]
        
    now_drugs = len(set(data_drugs["parent_molregno"]))
    now_drug_targets = len(set(data_drugs["tid"]))
    now_drug_targets_mutation = len(set(data_drugs["tid_mutation"]))
    now_drug_pairs = len(set(data_drugs['cpd_target_pair']))
    now_drug_pairs_mutation = len(set(data_drugs['cpd_target_pair_mutation']))

    return [now_mols, now_drugs, 
            now_targets, now_drug_targets,
            now_targets_mutation, now_drug_targets_mutation,
            now_pairs, now_drug_pairs,
            now_pairs_mutation, now_drug_pairs_mutation]


def add_dataset_sizes(data, label, 
                      all_lengths, all_lengths_pchembl):
    data_test = data.copy()
    all_lengths.append([label] + calculate_dataset_sizes(data_test))
    
    # restrict to data with any pchembl value (any data with a pchembl, even if it is based on only functional data)
    # these statistics are purely based on removing compound-target pairs without pchembl information
    # i.e., the subset of the dataset is determined by the given data parameter and not recalculated
    data_pchembl = data_test.dropna(subset=[x for x in data_test.columns if x.startswith('pchembl_value')], how = 'all')
    all_lengths_pchembl.append([label] + calculate_dataset_sizes(data_pchembl))