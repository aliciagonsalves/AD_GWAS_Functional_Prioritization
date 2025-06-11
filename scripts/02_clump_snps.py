# scripts/02_clump_snps.py

import pandas as pd
import os
import json

def perform_clumping(df: pd.DataFrame, window: int) -> pd.DataFrame:
    """
    Performs distance-based clumping on a dataframe of significant SNPs.
    
    This function iterates through SNPs, sorted by p-value, and selects 'lead'
    SNPs. Any subsequent SNP within a defined genomic window of a lead SNP on the
    same chromosome is 'clumped' with it and removed from consideration. This
    reduces correlated signals from Linkage Disequilibrium (LD) to a set of 
    independent loci.

    Args:
        df (pd.DataFrame): DataFrame of SNPs, must be sorted by p-value ascending.
        window (int): The clumping window size in base pairs.

    Returns:
        pd.DataFrame: A DataFrame containing only the independent lead SNPs.
    """
    # The input DataFrame is already sorted by p-value from the previous script.
    # Iterate through it to select lead SNPs.
    
    lead_snps = []
    snps_to_discard = set()

    # Iterate through each SNP, from most to least significant
    for index, lead_snp in df.iterrows():
        # If this SNP was already 'clumped' with a more significant lead SNP, skip it.
        if index in snps_to_discard:
            continue
        
        # If not clumped, this is a new lead SNP. Add it to our list.
        lead_snps.append(lead_snp)
        
        # Now, find all other SNPs in the original dataframe that are within the
        # clumping window of this new lead SNP and mark them to be discarded.
        clump_chr = lead_snp['CHR']
        clump_min_bp = lead_snp['BP'] - window
        clump_max_bp = lead_snp['BP'] + window
        
        # Identify indices of SNPs to discard
        clumped_indices = df[
            (df['CHR'] == clump_chr) &
            (df['BP'] >= clump_min_bp) &
            (df['BP'] <= clump_max_bp)
        ].index
        
        # Add these indices to our discard set
        snps_to_discard.update(clumped_indices)
            
    return pd.DataFrame(lead_snps)

def main():
    """
    Main function to execute the clumping workflow.
    """
    print("--- Running Script 02: Identifying Independent Loci via Clumping ---")

    # --- 1. Configuration ---
    clumping_window_kb = 500  # Window size in kilobases
    clumping_window = clumping_window_kb * 1000 # Convert to base pairs

    # Define input and output file paths
    tables_dir = 'results/tables'
    risk_file = os.path.join(tables_dir, 'significant_risk_snps.tsv')
    protective_file = os.path.join(tables_dir, 'significant_protective_snps.tsv')
    output_path_lead_risk = os.path.join(tables_dir, 'lead_risk_snps.tsv')
    output_path_lead_protective = os.path.join(tables_dir, 'lead_protective_snps.tsv')
    summary_stats_file = os.path.join(tables_dir, 'summary_stats.json')

    # --- 2. Load Data with Error Handling ---
    try:
        df_risk = pd.read_csv(risk_file, sep='\t')
        df_protective = pd.read_csv(protective_file, sep='\t')
        print(f"Loaded {len(df_risk)} risk SNPs and {len(df_protective)} protective SNPs.")
    except FileNotFoundError as e:
        print(f"ERROR: Could not find input file from Script 01. Please run it first. Details: {e}")
        return

    # --- 3. Perform Clumping ---
    print(f"Performing clumping with a {clumping_window_kb}kb window...")
    lead_risk_snps = perform_clumping(df_risk, clumping_window)
    lead_protective_snps = perform_clumping(df_protective, clumping_window)

    # --- 4. Report and Save Results ---
    print("\nClumping Complete:")
    print(f"  > Risk Loci: Reduced {len(df_risk)} SNPs to {len(lead_risk_snps)} independent lead SNPs.")
    print(f"  > Protective Loci: Reduced {len(df_protective)} SNPs to {len(lead_protective_snps)} independent lead SNPs.")
    
    lead_risk_snps.to_csv(output_path_lead_risk, sep='\t', index=False)
    lead_protective_snps.to_csv(output_path_lead_protective, sep='\t', index=False)
    print(f"Saved lead risk SNPs to: {output_path_lead_risk}")
    print(f"Saved lead protective SNPs to: {output_path_lead_protective}")

    # --- 5. Update Summary Statistics ---
    try:
        with open(summary_stats_file, 'r') as f:
            summary_data = json.load(f)

        summary_data['independent_lead_risk'] = len(lead_risk_snps)
        summary_data['independent_lead_protective'] = len(lead_protective_snps)

        with open(summary_stats_file, 'w') as f:
            json.dump(summary_data, f, indent=4)
        print(f"Updated summary statistics in: {summary_stats_file}")
    except FileNotFoundError:
        print(f"WARNING: {summary_stats_file} not found. Skipping update. Please run Script 01 first.")

    print("\n--- Script 02 Complete ---")

if __name__ == '__main__':
    main()