# scripts/10_generate_plots.py

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import json
import os
import networkx as nx
import plotly.graph_objects as go
import numpy as np

def generate_all_plots():
    """
    Main function to generate all final figures for the project.
    It reads the processed data tables and creates a suite of
    static (PNG) and interactive (HTML) plots.
    """
    print("--- Running Script 10: Generating Final Suite of All Visualizations ---")

    # --- 1. Style Definitions and Data Loading ---
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif']
    COLOR_RISK = "#B84406"        
    COLOR_PROTECTIVE = "#028086"  
    
    output_dir = 'results/figures'
    tables_dir = 'results/tables'
    os.makedirs(output_dir, exist_ok=True)

    print("Loading all necessary data files...")
    try:
        with open(os.path.join(tables_dir, 'summary_stats.json'), 'r') as f:
            summary = json.load(f)
        df_sharing = pd.read_csv(os.path.join(tables_dir, 'eqtl_sharing_by_snp.tsv'), sep='\t')
        df_enriched = pd.read_csv(os.path.join(tables_dir, 'candidates_enriched.tsv'), sep='\t')
        df_ld = pd.read_csv(os.path.join(tables_dir, 'all_ld_blocks.tsv'), sep='\t')
        gff_cols = ['seqid', 'type', 'start', 'end', 'attributes']
        df_genes = pd.read_csv('data/annotation/gencode.v48.annotation.gff3', sep='\t', comment='#', header=None, usecols=[0,2,3,4,8], names=gff_cols, compression='infer')
        df_genes['eGene'] = df_genes['attributes'].str.extract(r'gene_id=([^;]+)')[0].str.split('.').str[0]
    except FileNotFoundError as e:
        print(f"ERROR: A required data file is missing. Please run the full analysis pipeline first. Details: {e}")
        return

    # --- 2. Generate Plot 1: Two-Panel Funnel Plot ---
    print("\nGenerating 1/7: Two-Panel Funnel Plot...")
    stages = ['GWAS Significant', 'Independent Lead', 'eQTL-associated', 'Hi-C Supported']
    risk_keys = ['gwas_significant_risk', 'independent_lead_risk', 'eqtl_associated_risk', 'hic_supported_risk']
    prot_keys = ['gwas_significant_protective', 'independent_lead_protective', 'eqtl_associated_protective', 'hic_supported_protective']
    risk_counts = [summary.get(k, 0) for k in risk_keys]
    protective_counts = [summary.get(k, 0) for k in prot_keys]

    fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    fig1.suptitle('Project Funnel: Prioritizing AD-Associated SNPs', fontsize=18, weight='bold')

    ax1.plot(protective_counts, stages, color=COLOR_PROTECTIVE, marker='o', linestyle='--', markersize=8)
    ax1.set_title('Protective (Resilience) Loci', fontsize=14, color=COLOR_PROTECTIVE)
    ax1.set_xscale('log'); ax1.grid(axis='x', linestyle='--', alpha=0.6)
    for i, val in enumerate(protective_counts):
        ax1.text(val * 1.1, stages[i], f'{val:,}', va='center', ha='left', fontsize=11)

    ax2.plot(risk_counts, stages, color=COLOR_RISK, marker='o', linestyle='--', markersize=8)
    ax2.set_title('Risk Loci', fontsize=14, color=COLOR_RISK)
    ax2.set_xlabel('Number of Unique SNPs (Log Scale)', fontsize=12)
    ax2.grid(axis='x', linestyle='--', alpha=0.6)
    for i, val in enumerate(risk_counts):
        ax2.text(val * 0.9, stages[i], f'{val:,}', va='center', ha='right', fontsize=11)
    
    ax1.invert_yaxis(); ax2.invert_yaxis()
    max_val = max(max(risk_counts), max(protective_counts))
    ax1.set_xlim(left=0.8, right=max_val * 3)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(output_dir, '01_funnel_plot.png'), dpi=300)
    print("  > Funnel plot saved.")
    plt.close(fig1)


    # --- 3. Generate Plot 2: eQTL Sharing Bar Chart ---
    print("Generating 2/7: eQTL Sharing Bar Chart...")
    plt.figure(figsize=(12, 10))
    palette = {'Risk': COLOR_RISK, 'Protective': COLOR_PROTECTIVE}
    sns.barplot(x='tissue_count', y='SNP', data=df_sharing.head(25), hue='variant_type', palette=palette, dodge=False)
    plt.title('Top 25 SNPs by Widespread eQTL Activity', fontsize=16, weight='bold')
    plt.xlabel('Number of Brain Tissues with eQTL Effect', fontsize=12)
    plt.ylabel('Lead SNP', fontsize=12)
    plt.legend(title='Variant Type'); plt.grid(axis='x', linestyle='--', alpha=0.7); plt.tight_layout()
    plt.savefig(os.path.join(output_dir, '02_eqtl_sharing_ranking.png'), dpi=300)
    print("  > eQTL sharing plot saved.")
    plt.close()

    # --- 4. Generate Plot 3: eQTL Activity Heatmap ---
    print("Generating 3/7: eQTL Activity Heatmap...")
    df_eqtl = pd.read_csv(os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv'), sep='\t')
    df_eqtl['neg_log10_p'] = -np.log10(df_eqtl['eQTL_P_Value'])
    heatmap_df = df_eqtl.pivot_table(index='SNP', columns='tissue', values='neg_log10_p').fillna(0)

    # Filter and sort the data for a cleaner plot
    tissue_counts = df_eqtl.groupby('SNP')['tissue'].nunique()
    snps_to_keep = tissue_counts[tissue_counts >= 4].index
    # Ensure we only try to select SNPs that exist in the heatmap index
    snps_to_keep_existing = [snp for snp in snps_to_keep if snp in heatmap_df.index]
    if not snps_to_keep_existing:
        print("  > No SNPs passed the tissue count filter. Skipping eQTL heatmap.")
    else:
        heatmap_df = heatmap_df.loc[snps_to_keep_existing]
        heatmap_df = heatmap_df.loc[heatmap_df.sum(axis=1).sort_values(ascending=False).index]
        heatmap_df.columns = heatmap_df.columns.str.replace('_', ' ').str.title()

        # --- Create the Heatmap Visualization ---
        plt.figure(figsize=(15, 12))
        ax = sns.heatmap(
            heatmap_df, 
            cmap="mako",
            linewidths=.5,
            linecolor='lightgray',
            annot=False
        )
        
        plt.title('eQTL Activity of Top AD Loci Across Brain Tissues', fontsize=18, weight='bold')
        plt.xlabel('GTEx Brain Tissue', fontsize=12)
        plt.ylabel('Lead SNP', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        cbar = ax.collections[0].colorbar
        cbar.set_label('eQTL Association Strength (-log10 P)', fontsize=12, rotation=270, labelpad=20)
        
        snp_info = df_eqtl[['SNP', 'BETA']].drop_duplicates().set_index('SNP')
        for label in ax.get_yticklabels():
            snp_id = label.get_text()
            if snp_id in snp_info.index:
                beta = snp_info.loc[snp_id, 'BETA']
                label.set_color(COLOR_PROTECTIVE if beta < 0 else COLOR_RISK)
                label.set_weight('bold')

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, '03_eQTL_activity_heatmap.png'), dpi=300)
        print("  > eQTL activity heatmap saved.")

    plt.close()

    # --- 5. Generate Plot 4: HPA Cell-Type Heatmap ---
    print("Generating 4/7: HPA Cell-Type Heatmap...")

    # Define the column names we expect from the HPA enrichment script
    cell_type_cols = [
        'Single nuclei brain RNA - astrocyte [nTPM]',
        'Single nuclei brain RNA - central nervous system macrophage [nTPM]',
        'Single nuclei brain RNA - deep-layer intratelencephalic [nTPM]',
        'Single nuclei brain RNA - hippocampal CA1-3 [nTPM]',
        'Single nuclei brain RNA - hippocampal dentate gyrus [nTPM]',
        'Single nuclei brain RNA - oligodendrocyte [nTPM]'
    ]

    # Robustly check if ALL required columns are actually in the results file
    if not all(col in df_enriched.columns for col in cell_type_cols):
        print("  > WARNING: One or more required cell-type columns were not found. Skipping HPA heatmap.")
    else:
        # --- 1. Prepare Data for Heatmap ---
        # Select only the necessary columns for this plot
        heatmap_data = df_enriched[['Gene_Symbol', 'BETA'] + cell_type_cols].copy()
        heatmap_data.dropna(subset=['Gene_Symbol'], inplace=True)
        heatmap_data = heatmap_data.drop_duplicates(subset=['Gene_Symbol']).set_index('Gene_Symbol')

        # Clean up the column names to be shorter for the plot axes
        heatmap_data.columns = [
            'BETA', 'Astrocytes', 'Microglia', 'Cortical Neurons', 
            'Hippocampal Neurons (CA1-3)', 'Hippocampal Neurons (DG)', 'Oligodendrocytes'
        ]

        # Separate the BETA values for coloring the gene labels later
        betas = heatmap_data['BETA']
        heatmap_data = heatmap_data.drop(columns=['BETA'])

        # Fill any missing expression values with 0 and log-transform ONCE for better color scaling
        heatmap_data = heatmap_data.fillna(0).apply(np.log1p)

        # Sort rows (genes) by their total expression across all cell types for a cleaner look
        heatmap_data = heatmap_data.loc[heatmap_data.sum(axis=1).sort_values(ascending=False).index]

        # --- 2. Create the Heatmap Visualization ---
        plt.figure(figsize=(12, 10))
        ax = sns.heatmap(
            heatmap_data, 
            cmap="mako_r",  # A muted purple-green-teal colormap
            linewidths=.5,
            linecolor='lightgray',
            annot=True, fmt=".1f" # Annotate with values, rounded to 1 decimal place
        )
        
        plt.title('Brain Cell-Type Expression of Prioritized AD Genes', fontsize=18, weight='bold')
        plt.xlabel('Brain Cell Type', fontsize=12)
        plt.ylabel('Prioritized Gene', fontsize=12)
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)

        # Color the Y-axis labels (Gene Symbols) by risk/protective status
        for label in ax.get_yticklabels():
            gene_symbol = label.get_text()
            if gene_symbol in betas.index:
                beta_val = betas.loc[gene_symbol]
                label.set_color(COLOR_PROTECTIVE if beta_val < 0 else COLOR_RISK)
                label.set_weight('bold')

        # Adjust color bar label
        cbar = ax.collections[0].colorbar
        cbar.set_label('Gene Expression Level (log(1 + nTPM))', fontsize=12, rotation=270, labelpad=20)
        
        plt.tight_layout()
        
        # Save the final figure
        plt.savefig(os.path.join(output_dir, '04_hpa_cell_type_heatmap.png'), dpi=300)
        print("  > HPA cell-type heatmap saved.")
        plt.close()

# --- 6. Generate Plot 5: Gene Co-Expression Network ---
    print("Generating 5/7: Gene Co-Expression Network Plot...")
    # Use the enriched data, focusing on genes with a defined cluster
    network_df = df_enriched[['Gene_Symbol', 'BETA', 'Brain_Expression_Cluster']].dropna(subset=['Brain_Expression_Cluster', 'Gene_Symbol']).drop_duplicates(subset=['Gene_Symbol'])
    
    # Create a map of gene to its risk/protective status for coloring
    risk_map = {row.Gene_Symbol: 'Risk' if row.BETA > 0 else 'Protective' for i, row in network_df.iterrows()}

    # Create the graph object
    G = nx.Graph()
    
    # Add all genes from our list as nodes
    for gene in network_df['Gene_Symbol']:
        if gene in risk_map: # Ensure gene is in our map
            G.add_node(gene, type=risk_map[gene])

    # Add edges (connections) between genes that are in the same expression cluster
    clusters = network_df.groupby('Brain_Expression_Cluster')['Gene_Symbol'].apply(list)
    for cluster_name, genes_in_cluster in clusters.items():
        if len(genes_in_cluster) > 1 and cluster_name != 'NA':
            from itertools import combinations
            # Create an edge between every pair of genes in the cluster
            for gene1, gene2 in combinations(genes_in_cluster, 2):
                if G.has_node(gene1) and G.has_node(gene2):
                    G.add_edge(gene1, gene2)

    # For visual clarity, remove any genes that have no connections
    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)

    if len(G.nodes()) > 0:
        plt.figure(figsize=(14, 14))
        # Use a deterministic layout algorithm for consistent appearance
        pos = nx.kamada_kawai_layout(G)
        
        node_colors = [COLOR_RISK if G.nodes[n]['type'] == 'Risk' else COLOR_PROTECTIVE for n in G.nodes()]
        
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=4000, alpha=0.9)
        nx.draw_networkx_edges(G, pos, width=2.0, alpha=0.5, edge_color='grey')
        
        # Draw labels with a "halo" for better readability
        import matplotlib.patheffects as path_effects
        text_items = nx.draw_networkx_labels(G, pos, font_size=11, font_family='sans-serif', font_weight='bold')
        for _,t in text_items.items():
            t.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
        
        plt.title('Co-Expression Networks of Prioritized AD Genes', fontsize=20, weight='bold')
        plt.margins(0.1)
        plt.box(False)
        from matplotlib.lines import Line2D
        legend_elements = [Line2D([0], [0], marker='o', color='w', label='Risk-Associated Gene', markerfacecolor=COLOR_RISK, markersize=15),
                           Line2D([0], [0], marker='o', color='w', label='Protective-Associated Gene', markerfacecolor=COLOR_PROTECTIVE, markersize=15)]
        plt.legend(handles=legend_elements, loc='best', fontsize=12)
        plt.savefig(os.path.join(output_dir, '05_gene_network_plot.png'), dpi=300, bbox_inches='tight')
        print("  > Gene network plot saved.")
        plt.close()
    else:
        print("  > No networks with 2 or more genes found to plot.")


    # --- 7. Generate Plots 6 & 7: Interactive Locus Plots with Plotly ---
    print("Generating 6/7 & 7/7: Interactive Locus Plots...")
    
    # Define the latest Plotly JS CDN path
    PLOTLY_CDN_PATH = 'https://cdn.plot.ly/plotly-latest.min.js'

    def create_interactive_locus_plot(snp_info, ld_block, gene_info, variant_type_color):
        """Creates an interactive locus plot using Plotly."""
        fig = go.Figure()
        
        # Define plot range with some padding
        plot_start = min(ld_block['ld_start'], gene_info['start']) - 20000
        plot_end = max(ld_block['ld_end'], gene_info['end']) + 20000
        
        # Draw Chromosome Axis
        fig.add_shape(type="line", x0=plot_start, y0=0, x1=plot_end, y1=0, line=dict(color="black", width=2))

        # Draw LD Block
        fig.add_trace(go.Scatter(
            x=[ld_block['ld_start'], ld_block['ld_end']], y=[0, 0], mode='lines',
            line=dict(color=variant_type_color, width=20),
            text=f"<b>AD Locus (LD Block)</b><br>Lead SNP: {snp_info['SNP']}<br>Location: chr{ld_block['CHR']}:{ld_block['ld_start']}-{ld_block['ld_end']}",
            hoverinfo='text', name='LD Block'
        ))

        # Draw Gene Body
        fig.add_trace(go.Scatter(
            x=[gene_info['start'], gene_info['end']], y=[0, 0], mode='lines',
            line=dict(color="#6c757d", width=20),
            text=f"<b>Target eGene</b><br>{snp_info['Gene_Symbol']}<br>Location: {gene_info['seqid']}:{gene_info['start']}-{gene_info['end']}",
            hoverinfo='text', name='eGene'
        ))

        # Draw Hi-C Arc
        arc_path = f"M {ld_block['ld_start']},{0} C {(ld_block['ld_start']+plot_start)/2},{0.8} {(ld_block['ld_end']+plot_end)/2},{0.8} {gene_info['start']},{0}"
        fig.add_shape(type="path", path=arc_path, line=dict(color="grey", width=2, dash="dash"))
        
        variant_type = "Risk" if snp_info['BETA'] > 0 else "Protective"
        fig.update_layout(
            title=f"<b>Interactive View of Top {variant_type} Locus: {snp_info['SNP']} -> {snp_info['Gene_Symbol']}</b>",
            xaxis_title=f"Position on Chromosome {int(snp_info['CHR'])}",
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.5, 1]),
            showlegend=False, plot_bgcolor='white', font_family='sans-serif'
        )
        return fig

    # Generate and save plot for Top Protective Hit
    if not df_enriched[df_enriched['BETA'] < 0].empty:
        top_prot_hit = df_enriched[df_enriched['BETA'] < 0].sort_values(by='P').iloc[0]
        prot_ld_block = df_ld[df_ld['SNP'] == top_prot_hit['SNP']].iloc[0]
        prot_gene_info = df_genes[df_genes['eGene'] == top_prot_hit['eGene']].iloc[0]
        fig_prot = create_interactive_locus_plot(top_prot_hit, prot_ld_block, prot_gene_info, COLOR_PROTECTIVE)
        fig_prot.write_html(
            os.path.join(output_dir, "06_interactive_protective_locus.html"),
            include_plotlyjs=PLOTLY_CDN_PATH
        )
        print("  > Interactive protective locus plot saved.")
    else:
        print("  > No protective SNPs with Hi-C support found to plot.")

    # Generate and save plot for Top Risk Hit
    if not df_enriched[df_enriched['BETA'] > 0].empty:
        top_risk_hit = df_enriched[df_enriched['BETA'] > 0].sort_values(by='P').iloc[0]
        risk_ld_block = df_ld[df_ld['SNP'] == top_risk_hit['SNP']].iloc[0]
        risk_gene_info = df_genes[df_genes['eGene'] == top_risk_hit['eGene']].iloc[0]
        fig_risk = create_interactive_locus_plot(top_risk_hit, risk_ld_block, risk_gene_info, COLOR_RISK)
        fig_risk.write_html(
            os.path.join(output_dir, "07_interactive_risk_locus.html"),
            include_plotlyjs=PLOTLY_CDN_PATH
        )
        print("  > Interactive risk locus plot saved.")
    else:
        print("  > No risk SNPs with Hi-C support found to plot.")

    print("\n--- All visualizations generated successfully! ---")
    print("--- Project Complete ---")

if __name__ == '__main__':
    generate_all_plots()