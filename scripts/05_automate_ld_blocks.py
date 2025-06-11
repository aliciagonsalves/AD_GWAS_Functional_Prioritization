# scripts/05_automate_ld_blocks.py

import pandas as pd
import requests
import time
import os
from io import StringIO

def define_ld_blocks_from_api():
    """
    Automates the process of defining LD blocks for all lead eQTL SNPs.

    This script reads the list of eQTL SNPs, queries the LDlink API for each one
    to find its proxies in high LD (R² >= 0.8), and then calculates the
    genomic start and end coordinates of the resulting LD block.
    """
    print("--- Running Script 05: Automating LD Block Definition via API ---")

    # --- 1. Configuration ---
    # IMPORTANT: Paste your personal LDlink API token here.
    # To get a token, register at: https://ldlink.nci.nih.gov/?tab=apiaccess
    LDLINK_API_TOKEN = "YOUR_API_TOKEN_HERE" 

    # Define input and output files
    tables_dir = 'results/tables'
    eqtl_file = os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv')
    output_path = os.path.join(tables_dir, 'all_ld_blocks.tsv')

    # --- 2. Initial Checks ---
    if "YOUR_API_TOKEN_HERE" in LDLINK_API_TOKEN:
        print("ERROR: Please replace 'YOUR_API_TOKEN_HERE' with your actual LDlink API token in the script.")
        return

    try:
        df_eqtl = pd.read_csv(eqtl_file, sep='\t')
        unique_snps = df_eqtl['SNP'].unique()
        print(f"Found {len(unique_snps)} unique eQTL SNPs to process.")
    except FileNotFoundError:
        print(f"ERROR: Input eQTL file not found at '{eqtl_file}'. Please run Script 03 first.")
        return

    # --- 3. Automated API Calls to LDlink ---
    base_url = "https://ldlink.nih.gov/LDlinkRest/ldproxy"
    ld_block_results = []
    
    for i, snp_id in enumerate(unique_snps):
        print(f"Processing SNP {i+1}/{len(unique_snps)}: {snp_id}...")
        
        # Construct the API request URL, specifying the EUR population and GRCh38 build
        url = f"{base_url}?var={snp_id}&pop=EUR&r2_d=r2&genome_build=grch38&token={LDLINK_API_TOKEN}"
        
        try:
            # Make the API call with a timeout to prevent hanging indefinitely
            response = requests.get(url)
            response.raise_for_status() # Raise an error for bad responses (e.g., 404, 500)
            raw_response_text = response.text

            # Handle the specific case where the variant is not in the reference panel
            if "variant is not in 1000G reference panel" in raw_response_text.lower():
                print(f"  > WARNING: {snp_id} not found in 1000G reference panel. Skipping.")
                time.sleep(1)
                continue

            # Read the response text (which should be a TSV) into a pandas DataFrame
            df_proxy = pd.read_csv(StringIO(raw_response_text), sep='\t')
            
            # Handle other API-returned errors that come in table format
            if 'error' in df_proxy.columns:
                print(f"  > API Error for {snp_id}: {df_proxy['error'].iloc[0]}")
                continue

            # The API returns CHR and BP combined in a 'Coord' column, which we must parse
            if 'Coord' in df_proxy.columns:
                coord_split = df_proxy['Coord'].str.split(':', expand=True)
                df_proxy['Chr'] = coord_split[0].str.replace('chr', '')
                df_proxy['Position (GRCh38)'] = pd.to_numeric(coord_split[1], errors='coerce')
                df_proxy.dropna(subset=['Position (GRCh38)'], inplace=True)
            else:
                print(f"  > UNEXPECTED RESPONSE FORMAT for {snp_id}. 'Coord' column not found. Skipping.")
                continue

            # Filter for proxies in high LD (R² >= 0.8)
            df_filtered = df_proxy[df_proxy['R2'] >= 0.8]
            
            if not df_filtered.empty:
                # Calculate the min and max position to define the LD block
                chrom = df_filtered['Chr'].iloc[0]
                min_pos = int(df_filtered['Position (GRCh38)'].min())
                max_pos = int(df_filtered['Position (GRCh38)'].max())
                
                ld_block_results.append({'SNP': snp_id, 'CHR': chrom, 'ld_start': min_pos, 'ld_end': max_pos})
            else:
                print(f"  > No proxies found with R2 >= 0.8 for {snp_id}")

        except requests.exceptions.Timeout:
            print(f"  > ERROR: The request for {snp_id} timed out. Skipping.")
        except Exception as e:
            print(f"  > An unexpected error occurred for {snp_id}: {e}")
        
        # Pause for 1 second between API calls to be respectful to the server
        time.sleep(1)

    # --- 4. Save the Complete Results ---
    if ld_block_results:
        df_all_ld_blocks = pd.DataFrame(ld_block_results)
        df_all_ld_blocks.to_csv(output_path, sep='\t', index=False)
        
        print(f"\nSuccess! Defined LD blocks for {len(df_all_ld_blocks)} out of {len(unique_snps)} SNPs.")
        print(f"Complete results saved to: {output_path}")
    else:
        print("\nCould not define any LD blocks.")
        
    print("\n--- Script 05 Complete ---")

if __name__ == '__main__':
    define_ld_blocks_from_api()