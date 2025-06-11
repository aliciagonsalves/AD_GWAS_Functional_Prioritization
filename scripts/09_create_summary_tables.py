# scripts/09_create_summary_tables.py

import pandas as pd
import os
import glob

def create_gene_summary_table(enriched_file, output_path):
    """
    Creates a clean, human-readable summary table of the final prioritized genes
    and their most important biological annotations.
    """
    print("--- Creating Final Prioritized Gene Summary Table ---")
    
    try:
        df = pd.read_csv(enriched_file, sep='\t')
    except FileNotFoundError:
        print(f"ERROR: Enriched results file not found at {enriched_file}. Cannot create summary table.")
        return

    # Define the most impactful columns for the final summary
    summary_cols = {
        'Gene_Symbol': 'Gene', 'SNP': 'Lead SNP', 'BETA': 'Effect (BETA)',
        'Gene_Description': 'Description', 'Brain_Expression_Cluster': 'Brain Co-Expression Cluster',
        'Protein_Class': 'Protein Class', 'Subcellular_Location': 'Subcellular Location'
    }
    
    # Ensure all expected columns exist before trying to select them
    existing_cols = {k: v for k, v in summary_cols.items() if k in df.columns}
    
    # Keep only one entry per unique gene, picking the one with the best eQTL p-value
    df_summary = df.sort_values(by='eQTL_P_Value').drop_duplicates(subset='eGene', keep='first')
    df_summary = df_summary[list(existing_cols.keys())].rename(columns=existing_cols)
    df_summary.dropna(subset=['Gene'], inplace=True) # Remove rows without a gene symbol

    # Function to clean and truncate the list-like strings from HPA
    def clean_and_truncate(cell_text, max_items=2):
        if isinstance(cell_text, str) and cell_text.startswith('['):
            cleaned_text = cell_text.strip("[]' ")
            items = [item.strip().strip("'") for item in cleaned_text.split(',')]
            if len(items) > max_items:
                return ', '.join(items[:max_items]) + ', ...'
            else:
                return ', '.join(items)
        return cell_text

    # Apply cleaning to specific columns
    for col in ['Protein Class', 'Subcellular Location']:
        if col in df_summary.columns:
            df_summary[col] = df_summary[col].apply(clean_and_truncate)

    df_summary['Effect (BETA)'] = pd.to_numeric(df_summary['Effect (BETA)'], errors='coerce').round(3)
    
    # Save and print the table
    df_summary.to_csv(output_path, sep='\t', index=False)
    print(f"  > Success! Clean summary table saved to: {output_path}")
    print("\n--- Top Prioritized Genes Summary ---")
    try:
        print(df_summary.to_markdown(index=False))
    except ImportError:
        print("NOTE: Install 'tabulate' for Markdown output (`pip install tabulate`)")

def create_prs_summary_table(prs_dir, output_path):
    """
    Finds all PRS weight files and creates a summary table of their contents.
    """
    print("\n--- Creating PRS Definition Summary Table ---")
    
    search_pattern = os.path.join(prs_dir, 'prs_*.tsv')
    prs_files = glob.glob(search_pattern)

    if not prs_files:
        print(f"  > No PRS definition files found in '{prs_dir}'. Skipping.")
        return
        
    summary_list = []
    for file_path in prs_files:
        try:
            df = pd.read_csv(file_path, sep='\t')
            # Parse filename to get PRS tier and type
            filename = os.path.basename(file_path)
            parts = filename.replace('prs_', '').replace('_weights.tsv', '').split('_')
            variant_type = parts[-1].capitalize()
            prs_tier = " ".join(parts[:-1]).replace('_', ' ').capitalize()
            
            summary_list.append({
                'PRS Tier': prs_tier,
                'Variant Type': variant_type,
                'Number of SNPs': len(df)
            })
        except Exception:
            print(f"  > WARNING: Could not process file {file_path}.")

    if not summary_list:
        print("  > No valid PRS files could be processed.")
        return
        
    df_summary = pd.DataFrame(summary_list).sort_values(by=['PRS Tier', 'Variant Type'])
    df_summary.to_csv(output_path, sep='\t', index=False)
    print(f"  > Success! PRS summary table saved to: {output_path}")

    print("\n--- Summary of Defined Polygenic Risk Score Tiers ---")
    try:
        print(df_summary.to_markdown(index=False))
    except ImportError:
        print("NOTE: Install 'tabulate' for Markdown output (`pip install tabulate`)")

def main():
    """Main function to generate all summary tables."""
    tables_dir = 'results/tables'
    prs_dir = 'results/prs_definitions'
    
    # Ensure output directory exists
    os.makedirs(tables_dir, exist_ok=True)
    
    # Define paths
    enriched_file = os.path.join(tables_dir, 'candidates_enriched.tsv')
    gene_summary_output = os.path.join(tables_dir, 'summary_of_prioritized_genes.tsv')
    prs_summary_output = os.path.join(tables_dir, 'summary_of_prs_definitions.tsv')
    
    # Generate the tables
    create_gene_summary_table(enriched_file, gene_summary_output)
    create_prs_summary_table(prs_dir, prs_summary_output)
    
    print("\n--- Script 09 Complete ---")

if __name__ == '__main__':
    main()