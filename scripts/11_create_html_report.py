# scripts/11_create_html_report.py

import os

def create_master_report():
    """
    Combines the individual interactive Plotly figures into a single,
    standalone HTML report for easy viewing and presentation.
    """
    print("--- Running Script 11: Creating Final HTML Report ---")

    # --- 1. Configuration ---
    figures_dir = 'results/figures'
    # NOTE: These filenames match the output from the final visualization script
    protective_plot_file = os.path.join(figures_dir, "06_interactive_protective_locus.html")
    risk_plot_file = os.path.join(figures_dir, "07_interactive_risk_locus.html")
    master_report_file = os.path.join(figures_dir, "master_project_report.html")

    # --- 2. Read Individual Plot HTML ---
    try:
        print("Reading individual plot files...")
        with open(protective_plot_file, 'r', encoding='utf-8') as f:
            prot_html = f.read()
        with open(risk_plot_file, 'r', encoding='utf-8') as f:
            risk_html = f.read()
    except FileNotFoundError as e:
        print(f"ERROR: Could not find one of the input plot files. Please ensure Script 10 has run successfully. Details: {e}")
        return

    # --- 3. Create Master HTML Content ---
    print("Assembling master report...")
    # This string uses an f-string to embed the HTML content of each plot into a master page with professional styling.
    master_html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>AD GWAS Prioritization Project Report</title>
        <style>
            body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif; margin: 40px; line-height: 1.6; color: #333; background-color: #f8f9fa; }}
            h1, h2 {{ color: #2c3e50; border-bottom: 2px solid #dee2e6; padding-bottom: 10px; }}
            .plot-container {{ box-shadow: 0 4px 8px 0 rgba(0,0,0,0.1); border-radius: 8px; padding: 20px; margin-bottom: 40px; background-color: #ffffff; }}
        </style>
    </head>
    <body>

        <h1>Project Report: Top Prioritized Alzheimer's Loci</h1>
        <p>This report contains interactive visualizations for the top protective and top risk loci identified through an integrative genomics pipeline. Hover over plot elements for details.</p>

        <div class="plot-container">
            <h2>Top Protective Locus</h2>
            {prot_html}
        </div>

        <div class="plot-container">
            <h2>Top Risk Locus</h2>
            {risk_html}
        </div>

    </body>
    </html>
    """

    # --- 4. Write Master Report File ---
    with open(master_report_file, 'w', encoding='utf-8') as f:
        f.write(master_html_content)

    print(f"\nSuccess! A master report has been created.")
    print(f"Open this file in your web browser: {master_report_file}")
    
    print("\n--- Script 11 Complete ---")
    print("\n--- ENTIRE PROJECT PIPELINE COMPLETE ---")

if __name__ == '__main__':
    create_master_report()