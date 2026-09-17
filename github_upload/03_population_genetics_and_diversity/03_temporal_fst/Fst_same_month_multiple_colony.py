import pandas as pd
import numpy as np
from itertools import combinations


# Define the month order (if not already defined)
month_order = ['May', 'June', 'July', 'August', 'September', 'October', 'November', 'January', 'February']

# --- Fst Calculation Functions (same as before) ---
def calculate_nucleotide_diversity(df):
    if df.empty:
        return np.nan
    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)
    df['freq'] = df.apply(
        lambda row: float(
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ) if row['coverage'] > 0 else np.nan,
        axis=1
    )
    df['heterozygosity'] = (1 - df['freq']) * df['site_cov']
    return df['heterozygosity'].sum() / len(df)

def calculate_heterozygosity_for_combined(df):
    if df.empty:
        return df
    df['site_cov'] = df['coverage'] / (df['coverage'] - 1)
    df['freq'] = df.apply(
        lambda row: float(
            (row['A'] / row['coverage'])**2 +
            (row['T'] / row['coverage'])**2 +
            (row['C'] / row['coverage'])**2 +
            (row['G'] / row['coverage'])**2
        ) if row['coverage'] > 0 else np.nan,
        axis=1
    )
    df['heterozygosity'] = (1 - df['freq']) * df['site_cov']
    return df

def calculate_fst(df1, df2, pos_column='pos'):
    if df1.empty or df2.empty:
        return np.nan, np.nan

    theta1 = calculate_nucleotide_diversity(df1)
    theta2 = calculate_nucleotide_diversity(df2)
    hs = (theta1 + theta2) / 2

    combined_df = pd.merge(df1, df2, on=pos_column, how='outer', suffixes=('_1', '_2'))
    cols_to_sum = ['coverage', 'A', 'C', 'T', 'G']
    for col in cols_to_sum:
        combined_df[col] = combined_df[f'{col}_1'].fillna(0) + combined_df[f'{col}_2'].fillna(0)
        if f'{col}_1' in combined_df.columns:
            combined_df.drop(columns=[f'{col}_1', f'{col}_2'], inplace=True)
    combined_df = calculate_heterozygosity_for_combined(combined_df)
    ht = combined_df['heterozygosity'].sum() / len(combined_df) if len(combined_df) > 0 else np.nan

    fst_method1 = (ht - hs) / ht if ht > 0 else np.nan
    fst_method2 = (ht - hs) / ht if ht > 0 else np.nan # User's second method is the same formula

    return fst_method1, fst_method2

# Load the DataFrame
try:
    multiple_colony_shared_month_df = pd.read_csv('/home/robinch/projects/BEE-WITCH/22_anvio/multiple_colony_shared_month.csv')
    print("Loaded multiple_colony_shared_month.csv")
except FileNotFoundError:
    print("Error: multiple_colony_shared_month.csv not found. Please ensure the file is in the correct directory.")
    exit()

# Dictionary to store all Fst results
all_mag_id_fst_results = {}

# Iterate through each unique MAG_id
for mag_id in multiple_colony_shared_month_df['MAG_id'].unique():
    print(f"\nProcessing MAG_id: {mag_id}")
    mag_id_df = multiple_colony_shared_month_df[multiple_colony_shared_month_df['MAG_id'] == mag_id].copy()
    mag_id_results = []

    # Iterate through each unique month
    unique_months = sorted(mag_id_df['Month'].unique(), key=lambda m: month_order.index(m) if m in month_order else len(month_order))
    for month in unique_months:
        month_df = mag_id_df[mag_id_df['Month'] == month].copy()
        unique_colonies = month_df['Colony'].unique()

        # Perform pairwise Fst between each Colony combination
        for colony1, colony2 in combinations(unique_colonies, 2):
            df_colony1 = month_df[month_df['Colony'] == colony1].copy()
            df_colony2 = month_df[month_df['Colony'] == colony2].copy()

            # Ensure coverage and allele counts are numeric
            cols_to_numeric = ['coverage', 'A', 'C', 'T', 'G']
            for col in cols_to_numeric:
                if df_colony1[col].dtype == 'object':
                    df_colony1[col] = pd.to_numeric(df_colony1[col], errors='coerce')
                if df_colony2[col].dtype == 'object':
                    df_colony2[col] = pd.to_numeric(df_colony2[col], errors='coerce')

            fst1, fst2 = calculate_fst(df_colony1, df_colony2)
            mag_id_results.append({
                'MAG_id': mag_id,
                'Month': month,
                'Colony1': colony1,
                'Colony2': colony2,
                'Fst_Method1': fst1,
                'Fst_Method2': fst2
            })

    if mag_id_results:
        results_df = pd.DataFrame(mag_id_results)

        # Calculate average Fst per Month
        average_fst_per_month = results_df.groupby('Month')[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id} - Average Fst per Month:")
        print(average_fst_per_month)

        # Calculate average Fst per Colony comparison
        results_df['Colony_Comparison'] = results_df.apply(lambda row: tuple(sorted((row['Colony1'], row['Colony2']))), axis=1)
        average_fst_per_colony_pair = results_df.groupby('Colony_Comparison')[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id} - Average Fst per Colony Comparison:")
        print(average_fst_per_colony_pair)

        # Calculate average Fst across all comparisons
        average_fst_overall = results_df[['Fst_Method1', 'Fst_Method2']].mean()
        print(f"\n{mag_id} - Average Fst across all comparisons:")
        print(average_fst_overall)

        all_mag_id_fst_results[mag_id] = {
            'pairwise_fst': results_df.to_dict(orient='records'),
            'average_fst_per_month': average_fst_per_month.to_dict('index'),
            'average_fst_per_colony_pair': average_fst_per_colony_pair.to_dict('index'),
            'average_fst_overall': average_fst_overall.to_dict()
        }
    else:
        print(f"\nNo valid pairwise comparisons found for {mag_id}.")

# You can further process or save the all_mag_id_fst_results dictionary if needed

# 1. Save Pairwise Fst Results
pairwise_fst_list = []
for mag_id, results in all_mag_id_fst_results.items():
    for record in results['pairwise_fst']:
        pairwise_fst_list.append(record)

if pairwise_fst_list:
    pairwise_fst_df = pd.DataFrame(pairwise_fst_list)
    pairwise_fst_df.to_csv('pairwise_fst_results.csv', index=False)
    print("\nPairwise Fst results saved to pairwise_fst_results.csv")
else:
    print("\nNo pairwise Fst results to save.")

# 2. Save Average Fst per Month
average_fst_month_list = []
for mag_id, results in all_mag_id_fst_results.items():
    for month, fst_values in results['average_fst_per_month'].items():
        average_fst_month_list.append({
            'MAG_id': mag_id,
            'Month': month,
            'Average_Fst_Method1': fst_values.get('Fst_Method1'),
            'Average_Fst_Method2': fst_values.get('Fst_Method2')
        })

if average_fst_month_list:
    average_fst_month_df = pd.DataFrame(average_fst_month_list)
    average_fst_month_df.to_csv('average_fst_per_month.csv', index=False)
    print("Average Fst per Month saved to average_fst_per_month.csv")
else:
    print("No average Fst per month results to save.")

# 3. Save Average Fst per Colony Pair
average_fst_colony_pair_list = []
for mag_id, results in all_mag_id_fst_results.items():
    for colony_pair_tuple, fst_values in results['average_fst_per_colony_pair'].items():
        colony1, colony2 = colony_pair_tuple
        average_fst_colony_pair_list.append({
            'MAG_id': mag_id,
            'Colony1': colony1,
            'Colony2': colony2,
            'Average_Fst_Method1': fst_values.get('Fst_Method1'),
            'Average_Fst_Method2': fst_values.get('Fst_Method2')
        })

if average_fst_colony_pair_list:
    average_fst_colony_pair_df = pd.DataFrame(average_fst_colony_pair_list)
    average_fst_colony_pair_df.to_csv('average_fst_per_colony_pair.csv', index=False)
    print("Average Fst per Colony Pair saved to average_fst_per_colony_pair.csv")
else:
    print("No average Fst per colony pair results to save.")

# 4. Save Overall Average Fst per MAG_id
overall_average_fst_list = []
for mag_id, results in all_mag_id_fst_results.items():
    overall_fst_values = results['average_fst_overall']
    overall_average_fst_list.append({
        'MAG_id': mag_id,
        'Average_Fst_Method1': overall_fst_values.get('Fst_Method1'),
        'Average_Fst_Method2': overall_fst_values.get('Fst_Method2')
    })

if overall_average_fst_list:
    overall_average_fst_df = pd.DataFrame(overall_average_fst_list)
    overall_average_fst_df.to_csv('overall_average_fst_per_mag_id.csv', index=False)
    print("Overall Average Fst per MAG_id saved to overall_average_fst_per_mag_id.csv")
else:
    print("No overall average Fst results to save.")