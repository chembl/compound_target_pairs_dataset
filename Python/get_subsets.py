import os
import pandas as pd

# function to calculate and return the different subsets of interest
def get_data_subsets(data, min_nof_cpds, desc):
    if desc == 'B':
        drop_desc = 'BF'
    else:
        drop_desc = 'B'
    data = data.drop(columns=['pchembl_value_mean_'+drop_desc, 
                              'pchembl_value_max_'+drop_desc, 
                              'pchembl_value_median_'+drop_desc, 
                              'first_publication_cpd_target_pair_'+drop_desc, 
                              'first_publication_cpd_target_pair_w_pchembl_'+drop_desc, 
                              'LE_'+drop_desc,
                              'BEI_'+drop_desc,
                              'SEI_'+drop_desc,
                              'LLE_'+drop_desc]).drop_duplicates()
    
    # Restrict the dataset to targets with at least *min_nof_cpds* compounds with a pchembl value.
    comparator_counts = data[data['pchembl_value_mean_'+desc].notnull()].groupby(['tid_mutation'])['parent_molregno'].count()
    targets_w_enough_cpds = comparator_counts[comparator_counts >= min_nof_cpds].index.tolist()
    df_enough_cpds = data.query('tid_mutation in @targets_w_enough_cpds')
    
    # Restrict the dataset further to targets with at least one compound-target pair labelled as 'D_DT', 'C3_DT', 'C2_DT', 'C1_DT' or 'C0_DT', 
    # i.e., compound-target pairs with a known interactions.
    c_dt_d_dt_targets = set(df_enough_cpds[df_enough_cpds['DTI'].isin(['D_DT', 'C3_DT', 'C2_DT', 'C1_DT', 'C0_DT'])].tid_mutation.to_list())
    df_c_dt_d_dt = df_enough_cpds.query('tid_mutation in @c_dt_d_dt_targets')
    
    # Restrict the dataset further to targets with at least one compound-target pair labelled as 'D_DT', 
    # i.e., known drug-target interactions. 
    d_dt_targets = set(df_enough_cpds[df_enough_cpds['DTI'] == 'D_DT'].tid_mutation.to_list())
    df_d_dt = df_enough_cpds.query('tid_mutation in @d_dt_targets')
    
    return data, df_enough_cpds, df_c_dt_d_dt, df_d_dt

def write_output(df, filename, write_to_csv, write_to_excel, delimiter):
    """
    Write dataframe df to outputfile named filename.

    :return: Returns False if writing to excel was unsuccessful, True otherwise.
    :rtype: bool
    """
    if write_to_csv:
        df.to_csv(filename+".csv", sep = delimiter, index = False)
    if write_to_excel:
        try:
            with pd.ExcelWriter(filename + ".xlsx",engine='xlsxwriter') as writer: 
                writer.book.use_zip64()
                df.to_excel(writer, index = False)
        except ValueError as e: # full dataset may be too large to write to excel
            # remove empty file in case of error to avoid confusion
            if os.path.exists(filename + ".xlsx"):
                os.remove(filename + ".xlsx")
            print(e)
            return False
    return True



def write_BF_to_file(df_combined, chembl_version, min_nof_cpds_BF, path_results, write_BF, write_to_csv, write_to_excel, delimiter):
    # consider binding and functional assays
    # assay description = binding+functional
    desc = 'BF'
    df_combined_BF = df_combined.copy()
    df_combined_BF, df_combined_BF_enough_cpds, df_combined_BF_c_dt_d_dt, df_combined_BF_d_dt = get_data_subsets(df_combined_BF, min_nof_cpds_BF, desc)

    if write_BF:
        # note that this is almost identical to the full dataset which will be saved later on
        # however, the binding-related columns are dropped
        name_BF = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_BF")
        write_BF_success = write_output(df_combined_BF, name_BF, write_to_csv, write_to_excel, delimiter)

        name_BF_100 = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_BF_"+ str(min_nof_cpds_BF))
        write_BF_success &= write_output(df_combined_BF_enough_cpds, name_BF_100, write_to_csv, write_to_excel, delimiter)
        
        name_BF_100_c_dt_d_dt = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_BF_"+ str(min_nof_cpds_BF) + "_c_dt_d_dt")
        write_BF_success &= write_output(df_combined_BF_c_dt_d_dt, name_BF_100_c_dt_d_dt, write_to_csv, write_to_excel, delimiter)

        name_BF_100_d_dt = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_BF_"+ str(min_nof_cpds_BF) + "_d_dt")
        write_BF_success &= write_output(df_combined_BF_d_dt, name_BF_100_d_dt, write_to_csv, write_to_excel, delimiter)

        return df_combined_BF, df_combined_BF_enough_cpds, df_combined_BF_c_dt_d_dt, df_combined_BF_d_dt, \
            name_BF,  name_BF_100, name_BF_100_c_dt_d_dt, name_BF_100_d_dt, \
            write_BF_success
    else:
        return df_combined_BF, df_combined_BF_enough_cpds, df_combined_BF_c_dt_d_dt, df_combined_BF_d_dt, \
            "",  "", "", "", False

    # # TODO: include?
    # ############### TESTING: binding and functional assays ###############
    # add_dataset_sizes(df_combined_BF, "all assays")
    # add_dataset_sizes(df_combined_BF_enough_cpds, "all, >= 100")
    # add_dataset_sizes(df_combined_BF_c_dt_d_dt, "all, >= 100, c_dt and d_dt")
    # add_dataset_sizes(df_combined_BF_d_dt, "all, >= 100, d_dt")

    


def write_B_to_file(df_combined, chembl_version, min_nof_cpds_B, path_results, write_B, write_to_csv, write_to_excel, delimiter):
    # consider only binding assays
    # assay description = binding
    desc = 'B'
    df_combined_B = df_combined[df_combined['keep_for_binding'] == True].copy()
    df_combined_B, df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt = get_data_subsets(df_combined_B, min_nof_cpds_B, desc)

    if write_B:
        name_B = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_B")
        write_B_success = write_output(df_combined_B, name_B, write_to_csv, write_to_excel, delimiter)

        name_B_100 = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_B_"+ str(min_nof_cpds_B))
        write_B_success &= write_output(df_combined_B_enough_cpds, name_B_100, write_to_csv, write_to_excel, delimiter)

        name_B_100_c_dt_d_dt = os.path.join(path_results, path_results+"ChEMBL"+chembl_version+"_CTI_B_"+ str(min_nof_cpds_B) + "_c_dt_d_dt")
        write_B_success &= write_output(df_combined_B_c_dt_d_dt, name_B_100_c_dt_d_dt, write_to_csv, write_to_excel, delimiter)

        name_B_100_d_dt = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_B_"+ str(min_nof_cpds_B) + "_d_dt")
        write_B_success &= write_output(df_combined_B_d_dt, name_B_100_d_dt, write_to_csv, write_to_excel, delimiter)
    
        return df_combined_B, df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt, \
            name_B,  name_B_100, name_B_100_c_dt_d_dt, name_B_100_d_dt, \
            write_B_success
    else:
        return df_combined_B, df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt, \
            "",  "", "", "", False

    # # TODO: include?
    # ############### TESTING: binding assays ###############
    # add_dataset_sizes(df_combined_B, "binding")
    # add_dataset_sizes(df_combined_B_enough_cpds, "b, >= 100")
    # add_dataset_sizes(df_combined_B_c_dt_d_dt, "b, >= 100, c_dt and d_dt")
    # add_dataset_sizes(df_combined_B_d_dt, "b, >= 100, d_dt")

    


def write_full_dataset_to_file(df_combined, path_results, chembl_version, write_to_csv, write_to_excel, delimiter, write_full_dataset, 
                               df_combined_BF_enough_cpds, df_combined_BF_c_dt_d_dt, df_combined_BF_d_dt, min_nof_cpds_BF,
                               df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt, min_nof_cpds_B):
    issue_ctr = 0
    for df, name in zip([df_combined_BF_enough_cpds, 
                        df_combined_BF_c_dt_d_dt, 
                        df_combined_BF_d_dt, 
                        df_combined_B_enough_cpds, 
                        df_combined_B_c_dt_d_dt, 
                        df_combined_B_d_dt
                        ], 
                        ['BF_'+str(min_nof_cpds_BF), 
                        'BF_'+str(min_nof_cpds_BF)+'_c_dt_d_dt', 
                        'BF_'+str(min_nof_cpds_BF)+'_d_dt', 
                        'B_'+str(min_nof_cpds_B), 
                        'B_'+str(min_nof_cpds_B)+'_c_dt_d_dt', 
                        'B_'+str(min_nof_cpds_B)+'_d_dt']):
        df_combined[name] = False
        df_combined.loc[(df_combined.index.isin(df.index)), name] = True
        # check that filtering works
        if not df_combined[df_combined[name]==True][df.columns].equals(df):
            print("Problem with", name)
            issue_ctr += 1
            
    print("Number of problems:", issue_ctr)

    name_all = os.path.join(path_results, "ChEMBL"+chembl_version+"_CTI_all")
    if write_full_dataset: 
        write_full_success = write_output(df_combined, name_all, write_to_csv, write_to_excel, delimiter)
    
    return name_all, write_full_success
