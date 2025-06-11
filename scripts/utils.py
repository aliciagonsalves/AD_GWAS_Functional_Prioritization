# scripts/utils.py

import pandas as pd
import re
import os

def parse_gencode_to_promoters(gff_file_path: str, promoter_window: int = 2000) -> pd.DataFrame:
    """
    Parses a GENCODE GFF3 file to define promoter regions for all genes.

    A promoter is defined as a window (+/- promoter_window) around the
    Transcription Start Site (TSS) of a gene.

    Args:
        gff_file_path (str): The path to the gencode.vXX.annotation.gff3 file.
        promoter_window (int): The number of base pairs upstream and downstream of the TSS.

    Returns:
        pd.DataFrame: A DataFrame with promoter coordinates for each gene.
    """
    print("  > (Util) Parsing GENCODE annotations to define promoter regions...")
    try:
        gff_cols = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df_gff = pd.read_csv(gff_file_path, sep='\t', comment='#', header=None, names=gff_cols, compression='infer')
    except FileNotFoundError:
        print(f"    > ERROR: GENCODE file not found at {gff_file_path}")
        return pd.DataFrame()

    df_genes = df_gff[df_gff['type'] == 'gene'].copy()
    
    df_genes['TSS'] = df_genes.apply(lambda row: row['start'] if row['strand'] == '+' else row['end'], axis=1)

    def get_gene_id(attributes_str: str) -> str:
        """Extracts the Ensembl gene ID from the GFF3 attributes string."""
        match = re.search(r'gene_id=([^;]+)', attributes_str)
        return match.group(1).split('.')[0] if match else None

    df_genes['eGene'] = df_genes['attributes'].apply(get_gene_id)
    df_genes.dropna(subset=['eGene'], inplace=True)

    df_promoters = pd.DataFrame({
        'promoter_chr': df_genes['seqid'].str.replace('chr', ''),
        'promoter_start': df_genes['TSS'] - promoter_window,
        'promoter_end': df_genes['TSS'] + promoter_window,
        'eGene': df_genes['eGene']
    })
    print(f"    > (Util) Defined promoter regions for {len(df_promoters):,} genes.")
    return df_promoters
