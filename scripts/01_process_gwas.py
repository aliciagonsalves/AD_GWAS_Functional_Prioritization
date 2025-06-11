# scripts/01_process_gwas.py

import pandas as pd
import os
import json

def process_gwas_data():
    """
    Loads and processes raw GWAS summary statistics to extract genome-wide significant
    risk and protective variants.
    """
    print("--- Running Script 01: Processing GWAS Summary Statistics ---")

    # --- 1. Configuration ---
    # Define input file and the columns required for the analysis.
    gwas_file = 'data/gwas/harmonised.qc.tsv'
    required_cols = [
        'hm_rsid', 'hm_chrom', 'hm_pos', 
        'p_value', 'hm_beta', 'hm_effect_allele', 'hm_other_allele'
    ]
    
    # Define output directory and file paths
    output_dir = 'results/tables'
    summary_stats_file = os.path.join(output_dir, 'summary_stats.json')
    output_path_risk = os.path.join(output_dir, 'significant_risk_snps.tsv')
    output_path_protective = os.path.join(output_dir, 'significant_protective_snps.tsv')
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # --- 2. Load and Filter Data ---
    print(f"Loading GWAS data from: {gwas_file}")
    try:
        df = pd.read_csv(gwas_file, sep='\t', usecols=required_cols)
    except FileNotFoundError:
        print(f"ERROR: Input GWAS file not found at '{gwas_file}'. Please check the path.")
        return

    print(f"Loaded {len(df)} variants.")

    # Filter for genome-wide significant SNPs (P < 5e-8)
    significant_df = df[df['p_value'] < 5e-8].copy()
    print(f"Found {len(significant_df)} genome-wide significant SNPs (P < 5e-8).")

    # For consistency, rename columns to a simpler format
    significant_df.rename(columns={
        'hm_rsid': 'SNP',
        'hm_chrom': 'CHR',
        'hm_pos': 'BP',
        'p_value': 'P',
        'hm_beta': 'BETA',
        'hm_effect_allele': 'A1',
        'hm_other_allele': 'A2'
    }, inplace=True)

    # --- 3. Separate Risk and Protective Variants ---
    # Protective variants have a negative BETA (resilience factors)
    protective_snps = significant_df[significant_df['BETA'] < 0].sort_values(by='P')

    # Risk variants have a positive BETA
    risk_snps = significant_df[significant_df['BETA'] > 0].sort_values(by='P')

    print(f"Separated into {len(risk_snps)} risk SNPs and {len(protective_snps)} protective SNPs.")

    # --- 4. Save Processed Data ---
    risk_snps.to_csv(output_path_risk, sep='\t', index=False)
    protective_snps.to_csv(output_path_protective, sep='\t', index=False)
    print(f"Saved risk SNPs to: {output_path_risk}")
    print(f"Saved protective SNPs to: {output_path_protective}")
    
    # --- 5. Create and Save Initial Summary Statistics ---
    summary_data = {
        'gwas_significant_risk': len(risk_snps),
        'gwas_significant_protective': len(protective_snps)
    }
    with open(summary_stats_file, 'w') as f:
        json.dump(summary_data, f, indent=4)
    print(f"Saved initial counts to: {summary_stats_file}")
    
    print("\n--- Script 01 Complete ---")

if __name__ == '__main__':
    process_gwas_data()