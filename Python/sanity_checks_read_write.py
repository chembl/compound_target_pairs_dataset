import pandas as pd


def test_equality(current_file, read_file_name, assay_type, file_type_list, calculate_RDKit):
    curr_file_copy = current_file.copy().reset_index(drop=True)

    for file_type in file_type_list:
        if file_type == 'csv':
            try:
                read_file = pd.read_csv(read_file_name+".csv", sep=";",
                                        dtype={'mutation': 'str',
                                               'tid_mutation': 'str',
                                               'atc_level1': 'str',
                                               'target_class_l2': 'str',
                                               'ro3_pass': 'str',
                                               'molecular_species': 'str',
                                               'full_molformula': 'str',
                                               'standard_inchi': 'str',
                                               'standard_inchi_key': 'str',
                                               'canonical_smiles': 'str',
                                               'scaffold_w_stereo': 'str',
                                               'scaffold_wo_stereo': 'str',
                                               })
            except FileNotFoundError:
                print(read_file_name+".csv not found")
                continue
        elif file_type == 'xlsx':
            try:
                read_file = pd.read_excel(read_file_name+".xlsx")
            except FileNotFoundError:
                print(read_file_name+".xlsx not found")
                continue

        if assay_type == 'BF' or assay_type == 'all':
            read_file = read_file.astype({'first_publication_cpd_target_pair_BF': 'Int64',
                                          'first_publication_cpd_target_pair_w_pchembl_BF': 'Int64',
                                          })
        if assay_type == 'B' or assay_type == 'all':
            read_file = read_file.astype({'first_publication_cpd_target_pair_B': 'Int64',
                                          'first_publication_cpd_target_pair_w_pchembl_B': 'Int64',
                                          })
        read_file = read_file.astype({'first_approval': 'Int64',
                                      'usan_year': 'Int64',
                                      'first_publication_cpd': 'Int64',
                                      'hba': 'Int64',
                                      'hbd': 'Int64',
                                      'rtb': 'Int64',
                                      'num_ro5_violations': 'Int64',
                                      'aromatic_rings': 'Int64',
                                      'heavy_atoms': 'Int64',
                                      'hba_lipinski': 'Int64',
                                      'hbd_lipinski': 'Int64',
                                      'num_lipinski_ro5_violations': 'Int64',
                                      })
        if calculate_RDKit:
            read_file = read_file.astype({'num_aliphatic_carbocycles': 'Int64',
                                          'num_aliphatic_heterocycles': 'Int64',
                                          'num_aliphatic_rings': 'Int64',
                                          'num_aromatic_carbocycles': 'Int64',
                                          'num_aromatic_heterocycles': 'Int64',
                                          'num_aromatic_rings': 'Int64',
                                          'num_heteroatoms': 'Int64',
                                          'num_saturated_carbocycles': 'Int64',
                                          'num_saturated_heterocycles': 'Int64',
                                          'num_saturated_rings': 'Int64',
                                          'ring_count': 'Int64',
                                          'num_stereocentres': 'Int64',
                                          'aromatic_atoms': 'Int64',
                                          'aromatic_c': 'Int64',
                                          'aromatic_n': 'Int64',
                                          'aromatic_hetero': 'Int64',
                                          })

        print(read_file_name)
        print("{:5} file is ok: {}".format(
            file_type, read_file.equals(curr_file_copy)))
    print("----------")


def check_read_write(write_to_csv, write_to_excel, calculate_RDKit,
                     write_BF, df_combined_BF, df_combined_BF_enough_cpds, df_combined_BF_c_dt_d_dt, df_combined_BF_d_dt,
                     name_BF, name_BF_100, name_BF_100_c_dt_d_dt, name_BF_100_d_dt, write_BF_success,
                     write_B, df_combined_B, df_combined_B_enough_cpds, df_combined_B_c_dt_d_dt, df_combined_B_d_dt,
                     name_B, name_B_100, name_B_100_c_dt_d_dt, name_B_100_d_dt, write_B_success,
                     write_full_dataset, df_combined, name_all, write_full_success):
    # Check that output files can be written and read and are identical to the original dataframes.
    file_type_list = []
    if write_to_csv:
        file_type_list.append('csv')
    if write_to_excel:
        file_type_list.append('xlsx')

    # Some output was written
    if len(file_type_list) > 0:
        # binding + functional
        if write_BF:
            print("Check BF subset")
            test_equality(df_combined_BF, name_BF, 'BF',
                          file_type_list if write_BF_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_BF_enough_cpds, name_BF_100, 'BF',
                          file_type_list if write_BF_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_BF_c_dt_d_dt, name_BF_100_c_dt_d_dt, 'BF',
                          file_type_list if write_BF_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_BF_d_dt, name_BF_100_d_dt, 'BF',
                          file_type_list if write_BF_success else file_type_list[:-1], calculate_RDKit)

        # binding only
        if write_B:
            print("Check B subset")
            test_equality(df_combined_B, name_B, 'B',
                          file_type_list if write_B_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_B_enough_cpds, name_B_100, 'B',
                          file_type_list if write_B_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_B_c_dt_d_dt, name_B_100_c_dt_d_dt, 'B',
                          file_type_list if write_B_success else file_type_list[:-1], calculate_RDKit)
            test_equality(df_combined_B_d_dt, name_B_100_d_dt, 'B',
                          file_type_list if write_B_success else file_type_list[:-1], calculate_RDKit)

        # full dataset
        if write_full_dataset:
            print("Check full dataset")
            test_equality(df_combined, name_all, 'all',
                          file_type_list if write_full_success else file_type_list[:-1], calculate_RDKit)





