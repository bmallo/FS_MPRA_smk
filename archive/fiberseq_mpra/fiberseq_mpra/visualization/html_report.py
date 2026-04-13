#!/usr/bin/env python3
"""
HTML Report Generator for Fiberseq MPRA Analysis

This module generates interactive HTML reports for visualizing differential
footprint analysis results.

The report includes:
    - Summary dashboard with overview statistics
    - Interactive heatmaps for variant effects
    - Sortable/filterable results tables
    - Per-variant detail views
"""

import json
import logging
from typing import Dict, List, Optional
from pathlib import Path
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


def generate_html_report(
    results_df: pd.DataFrame,
    wt_matrix: 'FootprintMatrix',
    variant_matrices: Dict[str, 'FootprintMatrix'],
    config: 'Config',
    output_path: str,
) -> None:
    """
    Generate an interactive HTML report.
    
    Parameters:
    -----------
    results_df : pd.DataFrame
        DataFrame with all test results
    wt_matrix : FootprintMatrix
        Wild-type footprint matrix
    variant_matrices : Dict[str, FootprintMatrix]
        Dictionary of variant footprint matrices
    config : Config
        Analysis configuration
    output_path : str
        Path to write the HTML file
    """
    logger.info("Generating HTML report...")
    
    # Prepare data for JavaScript
    summary_data = _prepare_summary_data(results_df, wt_matrix, variant_matrices)
    variants_data = _prepare_variants_data(results_df, variant_matrices)
    wt_baseline_data = _prepare_wt_baseline(wt_matrix)
    
    # Generate HTML
    html_content = _generate_html_template(
        summary_data=summary_data,
        variants_data=variants_data,
        wt_baseline_data=wt_baseline_data,
        config=config,
    )
    
    # Write to file
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    logger.info(f"HTML report saved to {output_path}")


def _prepare_summary_data(
    results_df: pd.DataFrame,
    wt_matrix: 'FootprintMatrix',
    variant_matrices: Dict[str, 'FootprintMatrix']
) -> dict:
    """Prepare summary statistics for the report."""
    
    if len(results_df) == 0:
        return {
            'total_variants': 0,
            'total_tests': 0,
            'significant_tests': 0,
            'wt_read_count': wt_matrix.read_count,
        }
    
    # Count significant results by direction
    sig_df = results_df[results_df['significant']]
    gains = len(sig_df[sig_df['direction'] == 'gain'])
    losses = len(sig_df[sig_df['direction'] == 'loss'])
    
    # Get top hits
    top_hits = results_df.nsmallest(20, 'pvalue_adj')[
        ['variant_id', 'position', 'size_bin', 'log2_fc', 'pvalue_adj', 'direction']
    ].to_dict('records')
    
    return {
        'total_variants': len(variant_matrices),
        'total_tests': len(results_df),
        'significant_tests': len(sig_df),
        'significant_gains': gains,
        'significant_losses': losses,
        'wt_read_count': wt_matrix.read_count,
        'position_range': list(wt_matrix.position_range),
        'size_bins': [b.name for b in wt_matrix.size_bins],
        'top_hits': top_hits,
    }


def _prepare_variants_data(
    results_df: pd.DataFrame,
    variant_matrices: Dict[str, 'FootprintMatrix']
) -> List[dict]:
    """Prepare per-variant data for the report."""
    
    variants_data = []
    
    for var_id, var_matrix in variant_matrices.items():
        var_results = results_df[results_df['variant_id'] == var_id]
        
        sig_count = var_results['significant'].sum()
        
        # Create heatmap data (position x size_bin -> log2_fc)
        heatmap_data = []
        for _, row in var_results.iterrows():
            heatmap_data.append({
                'position': int(row['position']),
                'size_bin': row['size_bin'],
                'log2_fc': float(row['log2_fc']),
                'pvalue_adj': float(row['pvalue_adj']),
                'significant': bool(row['significant']),
                'wt_rate': float(row['wt_rate']),
                'var_rate': float(row['var_rate']),
            })
        
        variants_data.append({
            'variant_id': var_id,
            'read_count': var_matrix.read_count,
            'significant_count': int(sig_count),
            'heatmap_data': heatmap_data,
        })
    
    # Sort by number of significant hits
    variants_data.sort(key=lambda x: x['significant_count'], reverse=True)
    
    return variants_data


def _prepare_wt_baseline(wt_matrix: 'FootprintMatrix') -> dict:
    """Prepare WT baseline data for visualization."""
    
    rate_matrix = wt_matrix.get_rate_matrix()
    
    # Convert to list of dicts for JavaScript
    baseline_data = []
    for size_bin in rate_matrix.index:
        for position in rate_matrix.columns:
            baseline_data.append({
                'position': int(position),
                'size_bin': size_bin,
                'rate': float(rate_matrix.at[size_bin, position]),
            })
    
    return {
        'read_count': wt_matrix.read_count,
        'data': baseline_data,
    }


def _generate_html_template(
    summary_data: dict,
    variants_data: List[dict],
    wt_baseline_data: dict,
    config: 'Config',
) -> str:
    """Generate the complete HTML document."""
    
    # Convert data to JSON for embedding
    summary_json = json.dumps(summary_data)
    variants_json = json.dumps(variants_data)
    wt_baseline_json = json.dumps(wt_baseline_data)
    config_json = json.dumps(config.to_dict())
    
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Fiberseq MPRA Analysis Report</title>
    
    <!-- Plotly.js for interactive charts -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    
    <!-- DataTables for sortable tables -->
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
    <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
    
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background-color: #f5f5f5;
            color: #333;
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px 20px;
            margin-bottom: 30px;
            border-radius: 8px;
        }}
        
        header h1 {{
            font-size: 2em;
            margin-bottom: 10px;
        }}
        
        header p {{
            opacity: 0.9;
        }}
        
        .card {{
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 20px;
            margin-bottom: 20px;
        }}
        
        .card h2 {{
            color: #667eea;
            margin-bottom: 15px;
            padding-bottom: 10px;
            border-bottom: 2px solid #f0f0f0;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        
        .stat-box {{
            background: linear-gradient(135deg, #f5f7fa 0%, #e4e8eb 100%);
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .stat-box .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: #667eea;
        }}
        
        .stat-box .label {{
            color: #666;
            font-size: 0.9em;
        }}
        
        .stat-box.gain {{
            background: linear-gradient(135deg, #d4edda 0%, #c3e6cb 100%);
        }}
        
        .stat-box.gain .value {{
            color: #28a745;
        }}
        
        .stat-box.loss {{
            background: linear-gradient(135deg, #f8d7da 0%, #f5c6cb 100%);
        }}
        
        .stat-box.loss .value {{
            color: #dc3545;
        }}
        
        .tabs {{
            display: flex;
            border-bottom: 2px solid #e0e0e0;
            margin-bottom: 20px;
        }}
        
        .tab {{
            padding: 10px 20px;
            cursor: pointer;
            border: none;
            background: none;
            font-size: 1em;
            color: #666;
            transition: all 0.3s;
        }}
        
        .tab:hover {{
            color: #667eea;
        }}
        
        .tab.active {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            margin-bottom: -2px;
        }}
        
        .tab-content {{
            display: none;
        }}
        
        .tab-content.active {{
            display: block;
        }}
        
        #heatmap-container {{
            width: 100%;
            height: 500px;
        }}
        
        .variant-selector {{
            margin-bottom: 20px;
        }}
        
        .variant-selector select {{
            padding: 10px;
            font-size: 1em;
            border: 1px solid #ddd;
            border-radius: 4px;
            min-width: 300px;
        }}
        
        table.dataTable {{
            width: 100% !important;
        }}
        
        .significant-yes {{
            background-color: #d4edda !important;
        }}
        
        .log2fc-positive {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .log2fc-negative {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .filter-controls {{
            display: flex;
            gap: 20px;
            margin-bottom: 20px;
            flex-wrap: wrap;
        }}
        
        .filter-controls label {{
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        
        .filter-controls input[type="checkbox"] {{
            width: 18px;
            height: 18px;
        }}
        
        .legend {{
            display: flex;
            justify-content: center;
            gap: 30px;
            margin-top: 10px;
            font-size: 0.9em;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .legend-color {{
            width: 20px;
            height: 20px;
            border-radius: 3px;
        }}
        
        footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>Fiberseq MPRA Footprint Analysis</h1>
            <p>Differential chromatin architecture analysis</p>
        </header>
        
        <!-- Summary Statistics -->
        <div class="card">
            <h2>Summary Statistics</h2>
            <div class="stats-grid">
                <div class="stat-box">
                    <div class="value" id="stat-wt-reads">-</div>
                    <div class="label">WT Reads</div>
                </div>
                <div class="stat-box">
                    <div class="value" id="stat-variants">-</div>
                    <div class="label">Variants Analyzed</div>
                </div>
                <div class="stat-box">
                    <div class="value" id="stat-tests">-</div>
                    <div class="label">Statistical Tests</div>
                </div>
                <div class="stat-box">
                    <div class="value" id="stat-significant">-</div>
                    <div class="label">Significant Results</div>
                </div>
                <div class="stat-box gain">
                    <div class="value" id="stat-gains">-</div>
                    <div class="label">Footprint Gains</div>
                </div>
                <div class="stat-box loss">
                    <div class="value" id="stat-losses">-</div>
                    <div class="label">Footprint Losses</div>
                </div>
            </div>
        </div>
        
        <!-- Main Content Tabs -->
        <div class="card">
            <div class="tabs">
                <button class="tab active" data-tab="overview">Overview</button>
                <button class="tab" data-tab="variants">Variant Explorer</button>
                <button class="tab" data-tab="table">Results Table</button>
                <button class="tab" data-tab="baseline">WT Baseline</button>
            </div>
            
            <!-- Overview Tab -->
            <div id="overview" class="tab-content active">
                <h3>Top Significant Hits</h3>
                <table id="top-hits-table" class="display" style="width:100%">
                    <thead>
                        <tr>
                            <th>Variant</th>
                            <th>Position</th>
                            <th>Size Bin</th>
                            <th>Log2 FC</th>
                            <th>Adj. P-value</th>
                            <th>Direction</th>
                        </tr>
                    </thead>
                    <tbody id="top-hits-body">
                    </tbody>
                </table>
            </div>
            
            <!-- Variant Explorer Tab -->
            <div id="variants" class="tab-content">
                <div class="variant-selector">
                    <label for="variant-select">Select Variant: </label>
                    <select id="variant-select">
                        <option value="">-- Select a variant --</option>
                    </select>
                    <span id="variant-info" style="margin-left: 20px; color: #666;"></span>
                </div>
                <div id="heatmap-container"></div>
                <div class="legend">
                    <div class="legend-item">
                        <div class="legend-color" style="background: #dc3545;"></div>
                        <span>Footprint Loss (var &lt; WT)</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: #f5f5f5; border: 1px solid #ddd;"></div>
                        <span>No Change</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: #28a745;"></div>
                        <span>Footprint Gain (var &gt; WT)</span>
                    </div>
                </div>
            </div>
            
            <!-- Results Table Tab -->
            <div id="table" class="tab-content">
                <div class="filter-controls">
                    <label>
                        <input type="checkbox" id="filter-significant" checked>
                        Show only significant
                    </label>
                    <label>
                        Size bin:
                        <select id="filter-sizebin">
                            <option value="">All</option>
                        </select>
                    </label>
                </div>
                <table id="results-table" class="display" style="width:100%">
                    <thead>
                        <tr>
                            <th>Variant</th>
                            <th>Position</th>
                            <th>Size Bin</th>
                            <th>WT Rate</th>
                            <th>Var Rate</th>
                            <th>Log2 FC</th>
                            <th>P-value</th>
                            <th>Adj. P-value</th>
                            <th>Significant</th>
                        </tr>
                    </thead>
                    <tbody>
                    </tbody>
                </table>
            </div>
            
            <!-- WT Baseline Tab -->
            <div id="baseline" class="tab-content">
                <p style="margin-bottom: 15px;">
                    Wild-type footprint occupancy landscape based on <span id="baseline-reads">-</span> reads.
                </p>
                <div id="baseline-heatmap" style="width: 100%; height: 400px;"></div>
            </div>
        </div>
        
        <footer>
            Generated by Fiberseq MPRA Analysis Pipeline v0.1.0
        </footer>
    </div>
    
    <script>
        // Embedded data
        const summaryData = {summary_json};
        const variantsData = {variants_json};
        const wtBaselineData = {wt_baseline_json};
        const configData = {config_json};
        
        // Initialize on page load
        document.addEventListener('DOMContentLoaded', function() {{
            initializeSummary();
            initializeTabs();
            initializeVariantSelector();
            initializeResultsTable();
            initializeBaselineHeatmap();
        }});
        
        function initializeSummary() {{
            document.getElementById('stat-wt-reads').textContent = summaryData.wt_read_count?.toLocaleString() || '-';
            document.getElementById('stat-variants').textContent = summaryData.total_variants?.toLocaleString() || '-';
            document.getElementById('stat-tests').textContent = summaryData.total_tests?.toLocaleString() || '-';
            document.getElementById('stat-significant').textContent = summaryData.significant_tests?.toLocaleString() || '-';
            document.getElementById('stat-gains').textContent = summaryData.significant_gains?.toLocaleString() || '-';
            document.getElementById('stat-losses').textContent = summaryData.significant_losses?.toLocaleString() || '-';
            
            // Populate top hits table
            const tbody = document.getElementById('top-hits-body');
            if (summaryData.top_hits) {{
                summaryData.top_hits.forEach(hit => {{
                    const row = document.createElement('tr');
                    const fcClass = hit.log2_fc > 0 ? 'log2fc-positive' : 'log2fc-negative';
                    row.innerHTML = `
                        <td>${{hit.variant_id}}</td>
                        <td>${{hit.position}}</td>
                        <td>${{hit.size_bin}}</td>
                        <td class="${{fcClass}}">${{hit.log2_fc.toFixed(3)}}</td>
                        <td>${{hit.pvalue_adj.toExponential(2)}}</td>
                        <td>${{hit.direction}}</td>
                    `;
                    tbody.appendChild(row);
                }});
            }}
            
            $('#top-hits-table').DataTable({{
                pageLength: 10,
                order: [[4, 'asc']]
            }});
        }}
        
        function initializeTabs() {{
            document.querySelectorAll('.tab').forEach(tab => {{
                tab.addEventListener('click', function() {{
                    // Remove active class from all tabs
                    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
                    document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                    
                    // Add active class to clicked tab
                    this.classList.add('active');
                    document.getElementById(this.dataset.tab).classList.add('active');
                }});
            }});
        }}
        
        function initializeVariantSelector() {{
            const select = document.getElementById('variant-select');
            
            variantsData.forEach(variant => {{
                const option = document.createElement('option');
                option.value = variant.variant_id;
                option.textContent = `${{variant.variant_id}} (${{variant.significant_count}} significant)`;
                select.appendChild(option);
            }});
            
            select.addEventListener('change', function() {{
                if (this.value) {{
                    const variant = variantsData.find(v => v.variant_id === this.value);
                    if (variant) {{
                        document.getElementById('variant-info').textContent = 
                            `${{variant.read_count.toLocaleString()}} reads, ${{variant.significant_count}} significant hits`;
                        plotVariantHeatmap(variant);
                    }}
                }}
            }});
            
            // Auto-select first variant if available
            if (variantsData.length > 0) {{
                select.value = variantsData[0].variant_id;
                select.dispatchEvent(new Event('change'));
            }}
        }}
        
        function plotVariantHeatmap(variant) {{
            const data = variant.heatmap_data;
            
            // Get unique positions and size bins
            const positions = [...new Set(data.map(d => d.position))].sort((a, b) => a - b);
            const sizeBins = [...new Set(data.map(d => d.size_bin))];
            
            // Create z-values matrix
            const z = sizeBins.map(bin => 
                positions.map(pos => {{
                    const item = data.find(d => d.position === pos && d.size_bin === bin);
                    return item ? item.log2_fc : 0;
                }})
            );
            
            // Create hover text
            const hovertext = sizeBins.map(bin =>
                positions.map(pos => {{
                    const item = data.find(d => d.position === pos && d.size_bin === bin);
                    if (item) {{
                        return `Position: ${{pos}}<br>Size: ${{bin}}<br>Log2FC: ${{item.log2_fc.toFixed(3)}}<br>` +
                               `WT rate: ${{(item.wt_rate * 100).toFixed(2)}}%<br>Var rate: ${{(item.var_rate * 100).toFixed(2)}}%<br>` +
                               `Adj. p-value: ${{item.pvalue_adj.toExponential(2)}}<br>Significant: ${{item.significant ? 'Yes' : 'No'}}`;
                    }}
                    return '';
                }})
            );
            
            const trace = {{
                z: z,
                x: positions,
                y: sizeBins,
                type: 'heatmap',
                colorscale: [
                    [0, '#dc3545'],
                    [0.5, '#f5f5f5'],
                    [1, '#28a745']
                ],
                zmin: -3,
                zmax: 3,
                hovertext: hovertext,
                hoverinfo: 'text',
                colorbar: {{
                    title: 'Log2 Fold Change',
                    titleside: 'right'
                }}
            }};
            
            const layout = {{
                title: `Differential Footprints: ${{variant.variant_id}}`,
                xaxis: {{
                    title: 'Position',
                    tickmode: 'auto',
                    nticks: 20
                }},
                yaxis: {{
                    title: 'Footprint Size Bin'
                }},
                margin: {{ t: 50, b: 50, l: 100, r: 100 }}
            }};
            
            Plotly.newPlot('heatmap-container', [trace], layout, {{responsive: true}});
        }}
        
        function initializeResultsTable() {{
            // Populate size bin filter
            const sizeBinSelect = document.getElementById('filter-sizebin');
            const sizeBins = [...new Set(variantsData.flatMap(v => v.heatmap_data.map(d => d.size_bin)))];
            sizeBins.forEach(bin => {{
                const option = document.createElement('option');
                option.value = bin;
                option.textContent = bin;
                sizeBinSelect.appendChild(option);
            }});
            
            // Flatten all data for the table
            let allData = [];
            variantsData.forEach(variant => {{
                variant.heatmap_data.forEach(d => {{
                    allData.push({{
                        variant_id: variant.variant_id,
                        ...d
                    }});
                }});
            }});
            
            // Initialize DataTable
            const table = $('#results-table').DataTable({{
                data: allData,
                columns: [
                    {{ data: 'variant_id' }},
                    {{ data: 'position' }},
                    {{ data: 'size_bin' }},
                    {{ data: 'wt_rate', render: d => (d * 100).toFixed(2) + '%' }},
                    {{ data: 'var_rate', render: d => (d * 100).toFixed(2) + '%' }},
                    {{ 
                        data: 'log2_fc', 
                        render: d => `<span class="${{d > 0 ? 'log2fc-positive' : 'log2fc-negative'}}">${{d.toFixed(3)}}</span>`
                    }},
                    {{ data: 'pvalue_adj', render: d => d.toExponential(2) }},
                    {{ data: 'pvalue_adj', render: d => d.toExponential(2) }},
                    {{ data: 'significant', render: d => d ? 'Yes' : 'No' }}
                ],
                pageLength: 25,
                order: [[6, 'asc']],
                createdRow: function(row, data) {{
                    if (data.significant) {{
                        $(row).addClass('significant-yes');
                    }}
                }}
            }});
            
            // Filter handlers
            document.getElementById('filter-significant').addEventListener('change', function() {{
                if (this.checked) {{
                    table.column(8).search('Yes').draw();
                }} else {{
                    table.column(8).search('').draw();
                }}
            }});
            
            document.getElementById('filter-sizebin').addEventListener('change', function() {{
                table.column(2).search(this.value).draw();
            }});
            
            // Apply initial filter
            table.column(8).search('Yes').draw();
        }}
        
        function initializeBaselineHeatmap() {{
            document.getElementById('baseline-reads').textContent = wtBaselineData.read_count.toLocaleString();
            
            const data = wtBaselineData.data;
            
            // Get unique positions and size bins
            const positions = [...new Set(data.map(d => d.position))].sort((a, b) => a - b);
            const sizeBins = [...new Set(data.map(d => d.size_bin))];
            
            // Create z-values matrix
            const z = sizeBins.map(bin =>
                positions.map(pos => {{
                    const item = data.find(d => d.position === pos && d.size_bin === bin);
                    return item ? item.rate * 100 : 0;  // Convert to percentage
                }})
            );
            
            const trace = {{
                z: z,
                x: positions,
                y: sizeBins,
                type: 'heatmap',
                colorscale: 'Viridis',
                colorbar: {{
                    title: 'Occupancy (%)',
                    titleside: 'right'
                }}
            }};
            
            const layout = {{
                title: 'Wild-Type Footprint Landscape',
                xaxis: {{
                    title: 'Position',
                    tickmode: 'auto',
                    nticks: 20
                }},
                yaxis: {{
                    title: 'Footprint Size Bin'
                }},
                margin: {{ t: 50, b: 50, l: 100, r: 100 }}
            }};
            
            Plotly.newPlot('baseline-heatmap', [trace], layout, {{responsive: true}});
        }}
    </script>
</body>
</html>'''
    
    return html


if __name__ == "__main__":
    # Test the HTML generation
    print("HTML report generator module loaded successfully")
