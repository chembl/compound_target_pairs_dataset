def get_stats(data, column):
    print(column)
    print(f"{'Full dataset:' : <40} {len(set(data[column]))}")
    print(f"{'Comparators' : <40} {len(set(data[data['DTI'].isin(['DT'])][column]))}")
    print(f"{'Drugs' : <40} {len(set(data[data['DTI'] == 'D_DT'][column]))}")
    print(f"{'Clinical Candidates' : <40} {len(set(data[data['DTI'].isin(['C0_DT', 'C1_DT', 'C2_DT', 'C3_DT'])][column]))}")
    print(f"{'Clinical Candidates, max_phase = 3' : <40} {len(set(data[data['DTI'] == 'C3_DT'][column]))}")
    print(f"{'Clinical Candidates, max_phase = 2' : <40} {len(set(data[data['DTI'] == 'C2_DT'][column]))}")
    print(f"{'Clinical Candidates, max_phase = 1' : <40} {len(set(data[data['DTI'] == 'C1_DT'][column]))}")
    print(f"{'Clinical Candidates, max_phase < 1' : <40} {len(set(data[data['DTI'] == 'C0_DT'][column]))}")
    print()

def print_stats(df):
    get_stats(df, "parent_molregno")
    get_stats(df, "tid")
    get_stats(df, "tid_mutation")
    get_stats(df, "cpd_target_pair")
    get_stats(df, "cpd_target_pair_mutation")
