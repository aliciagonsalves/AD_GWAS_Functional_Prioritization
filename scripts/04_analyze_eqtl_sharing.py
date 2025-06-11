# scripts/04_analyze_eqtl_sharing.py

import pandas as pd
import os

def analyze_eqtl_sharing():
    """
    Analyzes the consolidated eQTL results to quantify how many brain tissues
    each lead SNP has a regulatory effect in. This helps identify SNPs with
    widespread vs. tissue-specific effects, a key metric for prioritization.
    This script generates a data table that will be used for visualization later.
    """
    print("--- Running Script 04: Analyzing eQTL Sharing Across Brain Tissues ---")

    # --- 1. Configuration ---
    tables_dir = 'results/tables'
    eqtl_file = os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv')
    output_path = os.path.join(tables_dir, 'eqtl_sharing_by_snp.tsv')

    # --- 2. Load Data ---
    try:
        df_eqtl = pd.read_csv(eqtl_file, sep='\t')
        print(f"Loaded {len(df_eqtl)} eQTL associations from: {eqtl_file}")
    except FileNotFoundError as e:
        print(f"ERROR: Could not find input file from Script 03. Please run it first. Details: {e}")
        return

    # --- 3. Perform Analysis ---
    # For each unique SNP, count the number of unique tissues it appears in.
    print("Calculating tissue count for each lead SNP...")
    if df_eqtl.empty:
        print("Input eQTL file is empty. Cannot perform analysis.")
        # Create an empty output file to allow pipeline to continue
        pd.DataFrame(columns=['SNP', 'tissue_count', 'BETA', 'variant_type']).to_csv(output_path, sep='\t', index=False)
        return
        
    tissue_counts = df_eqtl.groupby('SNP')['tissue'].nunique().reset_index()
    tissue_counts.rename(columns={'tissue': 'tissue_count'}, inplace=True)

    # Sort by the number of tissues, from most to least widespread
    tissue_counts = tissue_counts.sort_values(by='tissue_count', ascending=False)
    
    # Add the BETA value back in to determine risk/protective status
    # Drop the duplicates to get a single BETA value for each unique SNP
    snp_info = df_eqtl[['SNP', 'BETA']].drop_duplicates(subset=['SNP'])
    
    # Merge the counts with the SNP info
    final_df = pd.merge(tissue_counts, snp_info, on='SNP')
    
    # Add a human-readable variant type column
    final_df['variant_type'] = final_df['BETA'].apply(lambda x: 'Risk' if x > 0 else 'Protective')

    # --- 4. Save Results ---
    final_df.to_csv(output_path, sep='\t', index=False)
    print(f"\nAnalysis complete. Found tissue sharing counts for {len(final_df)} SNPs.")
    print(f"Results saved to: {output_path}")
    print("This file will be used to generate the eQTL sharing bar chart in the final visualization script.")
    
    print("\n--- Script 04 Complete ---")


if __name__ == '__main__':
    analyze_eqtl_sharing()