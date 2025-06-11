# scripts/03_annotate_eqtls.py

import pandas as pd
import os
import json
import glob

def annotate_snps_with_eqtls():
    """
    Annotates lead GWAS SNPs with eQTL data from GTEx brain tissues.
    
    This function discovers GTEx brain data files automatically, loads the lead SNPs
    from the previous step, and identifies which lead SNPs are significant eQTLs
    by matching their genomic coordinates. This provides the first layer of
    functional evidence for our prioritized loci.
    """
    print("--- Running Script 03: Annotating Lead SNPs with Brain eQTLs ---")

    # --- 1. Configuration ---
    tables_dir = 'results/tables'
    gtex_dir = 'data/eqtl'
    
    lead_risk_file = os.path.join(tables_dir, 'lead_risk_snps.tsv')
    lead_protective_file = os.path.join(tables_dir, 'lead_protective_snps.tsv')
    output_path = os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv')
    summary_stats_file = os.path.join(tables_dir, 'summary_stats.json')

    # --- 2. Automatically Discover GTEx Brain Tissue Files ---
    search_pattern = os.path.join(gtex_dir, 'Brain_*.signif_variant_gene_pairs.txt.gz')
    gtex_file_paths = glob.glob(search_pattern)

    if not gtex_file_paths:
        print(f"ERROR: No GTEx brain tissue files found in '{gtex_dir}' matching the pattern. Halting.")
        return

    gtex_files_dict = {}
    for path in gtex_file_paths:
        # Extracts a clean tissue name (e.g., 'Cortex') from the long filename
        filename = os.path.basename(path)
        key_name = filename.replace('Brain_', '').split('.v8.')[0]
        gtex_files_dict[key_name] = path
    print(f"Automatically discovered and prepared {len(gtex_files_dict)} GTEx brain tissue files.")

    # --- 3. Load Lead SNP Data ---
    try:
        df_lead_risk = pd.read_csv(lead_risk_file, sep='\t')
        df_lead_protective = pd.read_csv(lead_protective_file, sep='\t')
        all_lead_snps = pd.concat([df_lead_risk, df_lead_protective], ignore_index=True)
        print(f"Loaded a total of {len(all_lead_snps)} unique lead SNPs to check.")
    except FileNotFoundError as e:
        print(f"ERROR: Could not find lead SNP files from Script 02. Please run it first. Details: {e}")
        return

    # --- 4. Find eQTLs by Iterating Through Tissues ---
    all_eqtl_results = []
    print("\nProcessing GTEx brain tissues...")
    for tissue, file_path in gtex_files_dict.items():
        print(f"  - Processing {tissue}...")
        try:
            df_eqtl = pd.read_csv(
                file_path, sep='\t', 
                usecols=['variant_id', 'gene_id', 'pval_nominal'],
                dtype={'variant_id': str, 'gene_id': str}
            )
            
            # Parse CHR and BP from GTEx 'variant_id' for a robust coordinate-based merge
            split_id = df_eqtl['variant_id'].str.split('_', expand=True)
            df_eqtl['CHR'] = pd.to_numeric(split_id[0].str.replace('chr', ''), errors='coerce')
            df_eqtl['BP'] = pd.to_numeric(split_id[1], errors='coerce')
            df_eqtl.dropna(subset=['CHR', 'BP'], inplace=True)
            df_eqtl[['CHR', 'BP']] = df_eqtl[['CHR', 'BP']].astype(int)

            merged_df = pd.merge(all_lead_snps, df_eqtl, on=['CHR', 'BP'], how='inner')
            
            if not merged_df.empty:
                merged_df['tissue'] = tissue
                all_eqtl_results.append(merged_df)
                print(f"    > Found {len(merged_df)} eQTL associations for lead SNPs.")
            
        except FileNotFoundError:
            print(f"    > WARNING: GTEx file not found: {file_path}. Skipping.")
        except Exception as e:
            print(f"    > An error occurred processing file {file_path}: {e}")

    # --- 5. Consolidate, Save, and Summarize Results ---
    if not all_eqtl_results:
        print("\nNo eQTLs were found for any lead SNPs in the specified brain tissues.")
        print("\n--- Script 03 Complete ---")
        return

    final_eqtl_df = pd.concat(all_eqtl_results, ignore_index=True)
    
    # Clean up Ensembl gene IDs (e.g., 'ENSG0000012345.6' -> 'ENSG0000012345')
    final_eqtl_df['eGene'] = final_eqtl_df['gene_id'].str.split('.').str[0]
    final_eqtl_df.rename(columns={'pval_nominal': 'eQTL_P_Value'}, inplace=True)
    
    output_cols = ['SNP', 'CHR', 'BP', 'P', 'BETA', 'A1', 'A2', 'eGene', 'eQTL_P_Value', 'tissue']
    final_eqtl_df = final_eqtl_df[output_cols].sort_values(by=['P', 'eQTL_P_Value'])
    final_eqtl_df.to_csv(output_path, sep='\t', index=False)
    
    unique_eqtl_snps_count = final_eqtl_df['SNP'].nunique()
    print(f"\nSuccess! Found a total of {len(final_eqtl_df):,} eQTL associations.")
    print(f"A unique set of {unique_eqtl_snps_count} lead SNPs are eQTLs in at least one brain tissue.")
    print(f"Results saved to: {output_path}")

    # --- 6. Update Summary Statistics ---
    try:
        with open(summary_stats_file, 'r') as f:
            summary_data = json.load(f)

        risk_eqtl_count = final_eqtl_df[final_eqtl_df['BETA'] > 0]['SNP'].nunique()
        protective_eqtl_count = final_eqtl_df[final_eqtl_df['BETA'] < 0]['SNP'].nunique()
        summary_data['eqtl_associated_risk'] = risk_eqtl_count
        summary_data['eqtl_associated_protective'] = protective_eqtl_count

        with open(summary_stats_file, 'w') as f:
            json.dump(summary_data, f, indent=4)
        print(f"Updated summary statistics in: {summary_stats_file}")
    except FileNotFoundError:
        print(f"WARNING: {summary_stats_file} not found. Skipping update. Please run Script 01 first.")

    print("\n--- Script 03 Complete ---")

if __name__ == '__main__':
    annotate_snps_with_eqtls()