# scripts/07_fetch_hpa_data.py

import pandas as pd
import requests
import time
import os
from io import StringIO

def fetch_hpa_data():
    """
    Fetches rich annotation data for prioritized genes from the Human Protein Atlas (HPA) API.
    """
    print("--- Running Script 07: Enriching Candidate Genes with HPA Data ---")

    # --- 1. Configuration ---
    results_file = 'results/tables/prioritized_candidates.tsv'
    output_path = os.path.join('results/tables', 'candidates_enriched.tsv')

    # Define the specific data columns to download from HPA
    selected_columns = [
        'g', 'gd', 'pc', 'scml', 'ecbrain', 'brain_RNA__tau',
        'Brain_sn_RNA_astrocyte', 'Brain_sn_RNA_oligodendrocyte',
        'Brain_sn_RNA_central_nervous_system_macrophage', # Microglia
        'Brain_sn_RNA_hippocampal_CA1-3', 'Brain_sn_RNA_hippocampal_dentate_gyrus',
        'Brain_sn_RNA_deep-layer_intratelencephalic'
    ]
    columns_string = ",".join(selected_columns)
    base_url = "https://www.proteinatlas.org/api/search_download.php"
    all_gene_data = []

    # --- 2. Load Input Genes ---
    try:
        df_results = pd.read_csv(results_file, sep='\t')
        unique_genes_list = df_results['eGene'].unique()
        print(f"Found {len(unique_genes_list)} unique eGenes from final results to query.")
    except FileNotFoundError:
        print(f"ERROR: Input file not found: {results_file}. Please run Script 06 first.")
        return

    # --- 3. Query HPA API One Gene at a Time ---
    for i, gene_id in enumerate(unique_genes_list):
        print(f"Processing gene {i+1}/{len(unique_genes_list)}: {gene_id}")
        api_url = f"{base_url}?search={gene_id}&format=json&columns={columns_string}&compress=no"
        
        try:
            response = requests.get(api_url, timeout=60)
            response.raise_for_status()
            data = response.json()
            
            if data:
                gene_data = data[0]
                gene_data['eGene'] = gene_id
                all_gene_data.append(gene_data)
            else:
                print(f"  > WARNING: No data returned from HPA for {gene_id}")

        except requests.exceptions.RequestException as e:
            print(f"  > An error occurred during the API request for {gene_id}: {e}")
        
        time.sleep(0.5) # Be polite to the server

    # --- 4. Create, Clean, and Merge Final Enriched Table ---
    if not all_gene_data:
        print("\nCould not retrieve any data from HPA. Saving original results.")
        df_results.to_csv(output_path, sep='\t', index=False) # Save original file to this name
        return

    df_hpa = pd.DataFrame(all_gene_data)
    
    df_hpa.rename(columns={
        'Gene': 'Gene_Symbol', 'Gene description': 'Gene_Description',
        'Protein class': 'Protein_Class', 'Subcellular main location': 'Subcellular_Location',
        'Brain expression cluster': 'Brain_Expression_Cluster', 'TAU score - Brain': 'Brain_Tau_Score'
    }, inplace=True)
    
    df_final_enriched = pd.merge(df_results, df_hpa, on='eGene', how='left')
    
    df_final_enriched.to_csv(output_path, sep='\t', index=False)
    
    print(f"\nSuccess! Saved final, enriched results table to {output_path}")
    print("This table now contains rich functional and brain-specific data from HPA.")
    print("\n--- Script 07 Complete ---")

if __name__ == '__main__':
    fetch_hpa_data()