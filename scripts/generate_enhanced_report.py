#!/usr/bin/env python3
"""
Enhanced Report Generator for Drug Discovery Pipeline Results
Demonstrates the new comprehensive reporting system with compound-specific analysis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import json
from datetime import datetime
import base64
from io import BytesIO

def generate_enhanced_report(results_csv: str, manifest_json: str, output_dir: str = None):
    """Generate a comprehensive HTML report with enhanced compound analysis."""
    
    # Load data
    df = pd.read_csv(results_csv)
    with open(manifest_json) as f:
        manifest = json.load(f)
    
    # Create output directory
    if output_dir is None:
        output_dir = Path(results_csv).parent / "enhanced_report"
    else:
        output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Generate visualizations
    plots = {}
    
    # 1. Score Distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df['fused_score'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    plt.title('Fused Score Distribution', fontsize=16, fontweight='bold')
    plt.xlabel('Fused Score', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.grid(True, alpha=0.3)
    score_dist_path = output_dir / "score_distribution.png"
    plt.savefig(score_dist_path, dpi=300, bbox_inches='tight')
    plt.close()
    plots['score_dist'] = score_dist_path.name
    
    # 2. Top Compounds by Target
    plt.figure(figsize=(12, 8))
    top_compounds = df.nlargest(15, 'fused_score')
    colors = plt.cm.viridis(range(len(top_compounds)))
    bars = plt.barh(range(len(top_compounds)), top_compounds['fused_score'], color=colors)
    plt.yticks(range(len(top_compounds)), [f"{row['compound_id']}\n({row['target_gene']})" for _, row in top_compounds.iterrows()])
    plt.xlabel('Fused Score', fontsize=12)
    plt.title('Top 15 Compounds by Fused Score', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    top_compounds_path = output_dir / "top_compounds.png"
    plt.savefig(top_compounds_path, dpi=300, bbox_inches='tight')
    plt.close()
    plots['top_compounds'] = top_compounds_path.name
    
    # 3. Similarity vs Fused Score
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(df['tanimoto_max'], df['fused_score'], 
                         c=df['evidence_strength'], cmap='plasma', alpha=0.7, s=60)
    plt.colorbar(scatter, label='Evidence Strength')
    plt.xlabel('Tanimoto Similarity (Max)', fontsize=12)
    plt.ylabel('Fused Score', fontsize=12)
    plt.title('Similarity vs Fused Score (colored by Evidence)', fontsize=16, fontweight='bold')
    plt.grid(True, alpha=0.3)
    sim_vs_score_path = output_dir / "similarity_vs_score.png"
    plt.savefig(sim_vs_score_path, dpi=300, bbox_inches='tight')
    plt.close()
    plots['sim_vs_score'] = sim_vs_score_path.name
    
    # 4. Target Analysis
    plt.figure(figsize=(12, 6))
    target_counts = df['target_gene'].value_counts()
    bars = plt.bar(target_counts.index, target_counts.values, color='lightcoral')
    plt.title('Compounds per Target', fontsize=16, fontweight='bold')
    plt.xlabel('Target Gene', fontsize=12)
    plt.ylabel('Number of Compounds', fontsize=12)
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    target_analysis_path = output_dir / "target_analysis.png"
    plt.savefig(target_analysis_path, dpi=300, bbox_inches='tight')
    plt.close()
    plots['target_analysis'] = target_analysis_path.name
    
    # Generate HTML report
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Enhanced Drug Discovery Pipeline Report</title>
    <style>
        body {{ 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; 
            margin: 0; 
            padding: 20px; 
            background-color: #f5f5f5;
        }}
        .container {{ 
            max-width: 1200px; 
            margin: 0 auto; 
            background-color: white; 
            padding: 30px; 
            border-radius: 10px; 
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .header {{ 
            text-align: center; 
            margin-bottom: 30px; 
            padding-bottom: 20px; 
            border-bottom: 3px solid #2c3e50;
        }}
        .header h1 {{ 
            color: #2c3e50; 
            margin: 0; 
            font-size: 2.5em;
        }}
        .header p {{ 
            color: #7f8c8d; 
            margin: 10px 0 0 0; 
            font-size: 1.1em;
        }}
        .section {{ 
            margin: 30px 0; 
            padding: 20px; 
            background-color: #f8f9fa; 
            border-radius: 8px; 
            border-left: 4px solid #3498db;
        }}
        .section h2 {{ 
            color: #2c3e50; 
            margin-top: 0; 
            font-size: 1.8em;
        }}
        .metrics {{ 
            display: grid; 
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
            gap: 20px; 
            margin: 20px 0;
        }}
        .metric {{ 
            background-color: white; 
            padding: 20px; 
            border-radius: 8px; 
            text-align: center; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric h3 {{ 
            color: #e74c3c; 
            margin: 0 0 10px 0; 
            font-size: 2em;
        }}
        .metric p {{ 
            color: #7f8c8d; 
            margin: 0; 
            font-weight: bold;
        }}
        .plot {{ 
            text-align: center; 
            margin: 20px 0;
        }}
        .plot img {{ 
            max-width: 100%; 
            height: auto; 
            border-radius: 8px; 
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }}
        .compound-table {{ 
            overflow-x: auto; 
            margin: 20px 0;
        }}
        .compound-table table {{ 
            width: 100%; 
            border-collapse: collapse; 
            background-color: white; 
            border-radius: 8px; 
            overflow: hidden; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .compound-table th, .compound-table td {{ 
            padding: 12px; 
            text-align: left; 
            border-bottom: 1px solid #ecf0f1;
        }}
        .compound-table th {{ 
            background-color: #34495e; 
            color: white; 
            font-weight: bold;
        }}
        .compound-table tr:hover {{ 
            background-color: #f8f9fa;
        }}
        .smiles {{ 
            font-family: 'Courier New', monospace; 
            font-size: 0.9em; 
            color: #2c3e50;
        }}
        .highlight {{ 
            background-color: #fff3cd; 
            padding: 15px; 
            border-radius: 5px; 
            border-left: 4px solid #ffc107; 
            margin: 20px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Enhanced Drug Discovery Pipeline Report</h1>
            <p>Comprehensive Analysis of Pyrazole Conjugate Compounds</p>
            <p><strong>Run ID:</strong> {manifest.get('run_id', 'Unknown')} | <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="section">
            <h2>📊 Executive Summary</h2>
            <div class="metrics">
                <div class="metric">
                    <h3>{len(df)}</h3>
                    <p>Total Compound-Target Pairs</p>
                </div>
                <div class="metric">
                    <h3>{df['compound_id'].nunique()}</h3>
                    <p>Unique Compounds</p>
                </div>
                <div class="metric">
                    <h3>{df['target_gene'].nunique()}</h3>
                    <p>Target Genes</p>
                </div>
                <div class="metric">
                    <h3>{df['fused_score'].max():.3f}</h3>
                    <p>Best Fused Score</p>
                </div>
            </div>
            
            <div class="highlight">
                <strong>🎯 Key Finding:</strong> The pipeline successfully analyzed {df['compound_id'].nunique()} pyrazole conjugate compounds against {df['target_gene'].nunique()} lung cancer targets (EGFR, ERBB2, FGFR1). The highest scoring compound achieved a fused score of {df['fused_score'].max():.3f}.
            </div>
        </div>
        
        <div class="section">
            <h2>📈 Score Distribution Analysis</h2>
            <div class="plot">
                <img src="{plots['score_dist']}" alt="Score Distribution">
            </div>
            <p><strong>Analysis:</strong> The fused score distribution shows {len(df[df['fused_score'] > 0.2])} compounds with scores above 0.2, indicating promising drug-target interactions.</p>
        </div>
        
        <div class="section">
            <h2>🏆 Top Performing Compounds</h2>
            <div class="plot">
                <img src="{plots['top_compounds']}" alt="Top Compounds">
            </div>
        </div>
        
        <div class="section">
            <h2>🔬 Similarity vs Activity Analysis</h2>
            <div class="plot">
                <img src="{plots['sim_vs_score']}" alt="Similarity vs Score">
            </div>
            <p><strong>Insight:</strong> This plot reveals the relationship between chemical similarity (Tanimoto coefficient) and predicted activity (fused score), colored by evidence strength.</p>
        </div>
        
        <div class="section">
            <h2>🎯 Target Coverage Analysis</h2>
            <div class="plot">
                <img src="{plots['target_analysis']}" alt="Target Analysis">
            </div>
        </div>
        
        <div class="section">
            <h2>🧪 Detailed Compound Analysis</h2>
            <div class="compound-table">
                {df.head(20).to_html(classes='table', escape=False, index=False)}
            </div>
            <p><em>Showing top 20 results. Full dataset available in the CSV file.</em></p>
        </div>
        
        <div class="section">
            <h2>📋 Compound-Specific Insights</h2>
            <div class="highlight">
                <h3>🔍 Top Compound Analysis:</h3>
                <p><strong>Best Performing Compound:</strong> {df.loc[df['fused_score'].idxmax(), 'compound_id']}</p>
                <p><strong>SMILES:</strong> <span class="smiles">{df.loc[df['fused_score'].idxmax(), 'compound_smiles']}</span></p>
                <p><strong>Target:</strong> {df.loc[df['fused_score'].idxmax(), 'target_gene']} ({df.loc[df['fused_score'].idxmax(), 'target_uniprot']})</p>
                <p><strong>Evidence Strength:</strong> {df.loc[df['fused_score'].idxmax(), 'evidence_strength']:.3f}</p>
                <p><strong>Similarity Score:</strong> {df.loc[df['fused_score'].idxmax(), 'tanimoto_max']:.3f}</p>
            </div>
        </div>
        
        <div class="section">
            <h2>📁 Data Sources & Files</h2>
            <p><strong>Results CSV:</strong> {results_csv}</p>
            <p><strong>Manifest JSON:</strong> {manifest_json}</p>
            <p><strong>Report Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
"""
    
    # Save HTML report
    html_path = output_dir / "enhanced_report.html"
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    print(f"✅ Enhanced report generated successfully!")
    print(f"📁 Report location: {html_path}")
    print(f"📊 Generated {len(plots)} visualizations")
    print(f"🧪 Analyzed {len(df)} compound-target pairs")
    print(f"🎯 Found {df['compound_id'].nunique()} unique compounds")
    
    return html_path

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("Usage: python generate_enhanced_report.py <results.csv> <manifest.json> [output_dir]")
        sys.exit(1)
    
    results_csv = sys.argv[1]
    manifest_json = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else None
    
    generate_enhanced_report(results_csv, manifest_json, output_dir)
