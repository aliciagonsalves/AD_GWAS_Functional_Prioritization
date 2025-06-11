# scripts/06_hic_analysis.py

import pandas as pd
import os
import glob
import re
import json
from utils import parse_gencode_to_promoters

def load_hic_data(hic_dir: str) -> pd.DataFrame:
    """Loads and consolidates all Hi-C loop files from a directory."""
    print("Loading and consolidating Hi-C loop data...")
    hic_files = glob.glob(os.path.join(hic_dir, '*.bedpe.gz'))
    if not hic_files:
        raise FileNotFoundError(f"No Hi-C .bedpe.gz files found in '{hic_dir}'")

    hic_cols = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    df_hic_list = [pd.read_csv(f, sep='\t', header=None, usecols=range(6), names=hic_cols, comment='#') for f in hic_files]
    df_hic = pd.concat(df_hic_list, ignore_index=True).drop_duplicates()

    # Clean and standardize Hi-C data for merging
    for col in ['start1', 'end1', 'start2', 'end2']:
        df_hic[col] = pd.to_numeric(df_hic[col], errors='coerce')
    df_hic.dropna(inplace=True)
    df_hic['chr1'] = df_hic['chr1'].astype(str).str.replace('chr', '')
    df_hic['chr2'] = df_hic['chr2'].astype(str).str.replace('chr', '')
    print(f"Loaded a total of {len(df_hic):,} unique Hi-C loops from {len(hic_files)} files.")
    return df_hic

def main():
    """Main function to execute the Hi-C integration workflow."""
    print("\n--- Running Script 06: Integrating Chromatin Interaction (Hi-C) Data ---")

    # --- 1. Configuration ---
    tables_dir = 'results/tables'
    hic_dir = 'data/hic'
    gencode_dir = 'data/annotation'
    
    ld_block_file = os.path.join(tables_dir, 'all_ld_blocks.tsv')
    eqtl_file = os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv')
    gff_file_path = os.path.join(gencode_dir, 'gencode.v48.annotation.gff3')
    output_path = os.path.join(tables_dir, 'prioritized_candidates.tsv')
    summary_stats_file = os.path.join(tables_dir, 'summary_stats.json')

    # --- 2. Load All Input Data ---
    try:
        print("Loading previously generated results...")
        df_ld_blocks = pd.read_csv(ld_block_file, sep='\t')
        df_eqtl = pd.read_csv(eqtl_file, sep='\t')
        print(f"Loaded {len(df_ld_blocks)} LD blocks and {len(df_eqtl)} eQTL associations.")
        
        df_hic = load_hic_data(hic_dir)
        df_promoters = parse_gencode_to_promoters(gff_file_path, promoter_window=2000)
    except FileNotFoundError as e:
        print(f"ERROR: An input file was not found. Please ensure all previous scripts have run successfully. Details: {e}")
        return

    # --- 3. Find Supported Interactions ---
    print("\nSearching for Hi-C loops connecting LD blocks to eGene promoters...")
    df_merged = pd.merge(df_eqtl, df_ld_blocks, on='SNP')
    df_merged = pd.merge(df_merged, df_promoters, on='eGene')
    print(f"Testing {len(df_merged)} potential SNP-eGene-Promoter combinations...")

    hic_supported_results = []
    for index, row in df_merged.iterrows():
        ld_chr, ld_start, ld_end = str(row['CHR_y']), row['ld_start'], row['ld_end']
        promoter_chr, p_start, p_end = str(row['promoter_chr']), row['promoter_start'], row['promoter_end']

        # Condition 1: LD block in anchor 1, Promoter in anchor 2
        cond1 = df_hic[(df_hic['chr1'] == ld_chr) & (df_hic['start1'] < ld_end) & (df_hic['end1'] > ld_start) & (df_hic['chr2'] == promoter_chr) & (df_hic['start2'] < p_end) & (df_hic['end2'] > p_start)]
        # Condition 2: LD block in anchor 2, Promoter in anchor 1
        cond2 = df_hic[(df_hic['chr2'] == ld_chr) & (df_hic['start2'] < ld_end) & (df_hic['end2'] > ld_start) & (df_hic['chr1'] == promoter_chr) & (df_hic['start1'] < p_end) & (df_hic['end1'] > p_start)]

        if not cond1.empty or not cond2.empty:
            hic_supported_results.append(row.to_dict())

    # --- 4. Consolidate Results and Update Summary Stats ---
    if not hic_supported_results:
        print("\nNo Hi-C supported interactions were found for the tested LD blocks and eGene pairs.")
        hic_risk, hic_prot = 0, 0
    else:
        df_final = pd.DataFrame(hic_supported_results).drop_duplicates(subset=['SNP', 'eGene', 'tissue'])
        final_cols = ['SNP', 'CHR_x', 'BP', 'P', 'BETA', 'A1', 'A2', 'eGene', 'eQTL_P_Value', 'tissue']
        df_final = df_final[final_cols].rename(columns={'CHR_x': 'CHR'})
        df_final = df_final.sort_values(by=['P', 'eQTL_P_Value'])
        
        df_final.to_csv(output_path, sep='\t', index=False)
        print(f"\n--- Analysis Complete ---")
        print(f"Found {len(df_final)} high-confidence SNP-eGene pairs supported by Hi-C.")
        print(f"Results saved to: {output_path}")

        hic_risk = df_final[df_final['BETA'] > 0]['SNP'].nunique()
        hic_prot = df_final[df_final['BETA'] < 0]['SNP'].nunique()

    # Update the summary_stats.json file
    try:
        with open(summary_stats_file, 'r') as f:
            summary_data = json.load(f)
        
        summary_data['hic_supported_risk'] = hic_risk
        summary_data['hic_supported_protective'] = hic_prot

        with open(summary_stats_file, 'w') as f:
            json.dump(summary_data, f, indent=4)
        print(f"Updated summary statistics in: {summary_stats_file}")
    except FileNotFoundError:
        print(f"WARNING: {summary_stats_file} not found. Cannot update Hi-C counts.")

    print("\n--- Script 06 Complete ---")

if __name__ == '__main__':
    main()