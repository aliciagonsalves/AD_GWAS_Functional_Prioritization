# scripts/08_define_prs_weights.py

import pandas as pd
import os

def define_prs_weights():
    """
    Packages the project's findings into tiered Polygenic Risk Score (PRS)
    definition files. Each file contains the SNPs and their effect sizes (BETAs),
    which serve as weights for building predictive models.
    """
    print("--- Running Script 08: Defining Polygenic Risk Score Weights ---")

    # --- 1. Configuration ---
    tables_dir = 'results/tables'
    output_dir = 'results/prs_definitions'
    os.makedirs(output_dir, exist_ok=True)

    # Input files from previous analysis steps
    lead_risk_file = os.path.join(tables_dir, 'lead_risk_snps.tsv')
    lead_protective_file = os.path.join(tables_dir, 'lead_protective_snps.tsv')
    eqtl_file = os.path.join(tables_dir, 'lead_snps_with_brain_eqtls.tsv')

    # --- 2. Define PRS #1: GWAS-Tier ---
    # This score is based on all independent lead SNPs identified after clumping.
    print("Defining GWAS-Tier PRS weights...")
    try:
        df_risk = pd.read_csv(lead_risk_file, sep='\t')
        df_prot = pd.read_csv(lead_protective_file, sep='\t')

        # Save only the SNP and its effect size (BETA)
        df_risk[['SNP', 'BETA']].to_csv(os.path.join(output_dir, 'prs_gwas_tier_risk.tsv'), sep='\t', index=False)
        df_prot[['SNP', 'BETA']].to_csv(os.path.join(output_dir, 'prs_gwas_tier_protective.tsv'), sep='\t', index=False)
        print("  > GWAS-Tier files created.")
    except FileNotFoundError as e:
        print(f"  > WARNING: Could not find input file from Script 02, skipping GWAS-Tier. Details: {e}")

    # --- 3. Define PRS #2: eQTL-Informed ---
    # This score is based only on the subset of SNPs with known function as brain eQTLs.
    print("Defining eQTL-Informed PRS weights...")
    try:
        df_eqtl = pd.read_csv(eqtl_file, sep='\t')

        # Filter for risk SNPs (BETA > 0)
        eqtl_risk_weights = df_eqtl[df_eqtl['BETA'] > 0][['SNP', 'BETA']].drop_duplicates()
        eqtl_risk_weights.to_csv(os.path.join(output_dir, 'prs_eqtl_informed_risk.tsv'), sep='\t', index=False)

        # Filter for protective SNPs (BETA < 0)
        eqtl_prot_weights = df_eqtl[df_eqtl['BETA'] < 0][['SNP', 'BETA']].drop_duplicates()
        eqtl_prot_weights.to_csv(os.path.join(output_dir, 'prs_eqtl_informed_protective.tsv'), sep='\t', index=False)
        print("  > eQTL-Informed files created.")
    except FileNotFoundError as e:
        print(f"  > WARNING: Could not find input file from Script 03, skipping eQTL-Informed Tier. Details: {e}")

    # --- 4. Note on Future PRS Tiers ---
    print("\nConceptual PRS Tiers for future work have also been established:")
    print("  > A 'Hi-C-Informed PRS' can be built using 'prioritized_candidates.tsv'.")
    print("  > A 'Cell-Type-Specific PRS' can be built using 'candidates_enriched.tsv'.")

    print("\n--- Script 08 Complete ---")

if __name__ == '__main__':
    define_prs_weights()