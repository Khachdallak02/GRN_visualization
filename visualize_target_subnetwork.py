"""
Visualize Target-Centered Subnetwork from Gene Regulatory Network (GRN).

This script:
1. Loads GRN data from CSV or H5AD (GRnnData) files
2. Selects a target gene
3. Finds top N TFs (default: 10) with highest connection scores to the target
4. Finds all connections between these TFs
5. Visualizes the subnetwork interactively using Dash Cytoscape
"""

import os
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import warnings

warnings.filterwarnings('ignore')

# Configuration
DEFAULT_GRN_FILE = 'pyscenic_output/grn_heart_failure_pyscenic.csv'
OUTPUT_DIR = 'figures'
MAX_TFS = 10  # Maximum number of TFs to include

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)


def load_grn_from_h5ad(grn_file):
    """Load GRN data from H5AD file (GRnnData format)"""
    try:
        from grnndata import GRNAnnData
        import anndata as ad
    except ImportError as e:
        print("ERROR: GRnnData package not installed.")
        print("Install it with: pip install grnndata")
        print(f"Import error: {e}")
        return None
    
    # Try to import read_h5ad if available
    try:
        from grnndata import read_h5ad
        has_read_h5ad = True
    except (ImportError, AttributeError):
        has_read_h5ad = False
    
    print(f"Loading GRN from H5AD file: {grn_file}")
    try:
        # Try multiple methods to load GRNAnnData
        grn_data = None
        
        # Method 1: Try read_h5ad function from grnndata (if available)
        if has_read_h5ad:
            try:
                grn_data = read_h5ad(grn_file)
                print("✓ Loaded using grnndata.read_h5ad()")
            except Exception as e:
                print(f"  Note: read_h5ad() failed: {e}")
                grn_data = None
        
        # Method 2: Load with anndata and wrap in GRNAnnData (fallback)
        if grn_data is None:
            try:
                adata = ad.read_h5ad(grn_file)
                grn_data = GRNAnnData(adata)
                print("✓ Loaded using anndata.read_h5ad() and wrapped in GRNAnnData")
            except Exception as e2:
                print(f"ERROR: Failed to load H5AD file: {e2}")
                raise ValueError(f"Could not load H5AD file. Error: {e2}")
        
        if grn_data is None:
            print("ERROR: Failed to load GRNAnnData object")
            return None
        
        print(f"✓ Successfully loaded GRNAnnData object")
        
        # Get adjacency matrix from GRN
        adj_matrix = None
        
        # Try varp['GRN'] first (most common)
        if hasattr(grn_data, 'varp') and 'GRN' in grn_data.varp:
            adj_matrix = grn_data.varp['GRN']
            print("✓ Found GRN in varp['GRN']")
        # Try grn attribute
        elif hasattr(grn_data, 'grn') and grn_data.grn is not None:
            adj_matrix = grn_data.grn
            print("✓ Found GRN in grn attribute")
        # Try varp with different key names
        elif hasattr(grn_data, 'varp'):
            varp_keys = list(grn_data.varp.keys())
            print(f"Available varp keys: {varp_keys}")
            # Look for any matrix-like data in varp (including sparse matrices)
            for key in varp_keys:
                val = grn_data.varp[key]
                # Check if it's a matrix (numpy array or sparse matrix)
                if isinstance(val, np.ndarray) or hasattr(val, 'toarray') or hasattr(val, 'shape'):
                    adj_matrix = val
                    print(f"✓ Found GRN in varp['{key}']")
                    break
        
        if adj_matrix is None:
            print("ERROR: Could not find GRN matrix in H5AD file")
            if hasattr(grn_data, 'varp'):
                print(f"Available varp keys: {list(grn_data.varp.keys())}")
            return None
        
        # Check matrix type and shape before conversion
        print(f"✓ GRN matrix type: {type(adj_matrix)}")
        
        # Handle sparse matrices first (MUST be done before numpy array conversion)
        if hasattr(adj_matrix, 'toarray'):
            print("  Detected sparse matrix, converting to dense...")
            try:
                # Get shape before conversion for debugging
                sparse_shape = adj_matrix.shape
                sparse_format = type(adj_matrix).__name__
                print(f"  Sparse matrix shape: {sparse_shape}")
                print(f"  Sparse matrix format: {sparse_format}")
                
                # Convert sparse to dense
                adj_matrix = adj_matrix.toarray()
                print(f"  ✓ Converted sparse matrix to dense array")
            except Exception as e:
                print(f"  ERROR: Failed to convert sparse matrix: {e}")
                import traceback
                traceback.print_exc()
                return None
        
        # Now ensure it's a numpy array (after sparse conversion)
        if not isinstance(adj_matrix, np.ndarray):
            try:
                adj_matrix = np.array(adj_matrix)
                print(f"  ✓ Converted to numpy array")
            except Exception as e:
                print(f"  ERROR: Failed to convert to numpy array: {e}")
                import traceback
                traceback.print_exc()
                return None
        
        print(f"✓ GRN matrix shape: {adj_matrix.shape}")
        print(f"✓ GRN matrix dtype: {adj_matrix.dtype}")
        
        # Ensure it's 2D
        if adj_matrix.ndim != 2:
            print(f"ERROR: Expected 2D matrix, got {adj_matrix.ndim}D")
            print(f"  Matrix shape: {adj_matrix.shape}")
            print(f"  Matrix type: {type(adj_matrix)}")
            return None
        
        # Get gene names
        gene_names = grn_data.var_names.tolist()
        print(f"✓ Found {len(gene_names)} genes")
        
        # Convert adjacency matrix to edge list
        print("Converting adjacency matrix to edge list...")
        
        edges = []
        n_genes = len(gene_names)
        
        # Check if matrix dimensions match gene count
        if adj_matrix.shape[0] != n_genes or adj_matrix.shape[1] != n_genes:
            print(f"WARNING: Matrix shape {adj_matrix.shape} doesn't match gene count {n_genes}")
            print(f"  Using min dimension: {min(adj_matrix.shape[0], adj_matrix.shape[1], n_genes)}")
            n_genes = min(adj_matrix.shape[0], adj_matrix.shape[1], n_genes)
        
        for i in range(n_genes):
            for j in range(n_genes):
                weight = adj_matrix[i, j]
                # Check for non-zero, non-NaN values
                if weight != 0 and not np.isnan(weight) and not np.isinf(weight):
                    edges.append({
                        'TF': gene_names[i],
                        'Target': gene_names[j],
                        'Correlation': float(weight)
                    })
        
        grn_df = pd.DataFrame(edges)
        print(f"✓ Extracted {len(grn_df)} edges from adjacency matrix")
        print(f"\nFirst few edges:")
        print(grn_df.head())
        
        return grn_df
        
    except Exception as e:
        print(f"ERROR: Failed to load H5AD file: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_grn_data(grn_file):
    """Load GRN data from CSV or H5AD file"""
    if not os.path.exists(grn_file):
        print(f"ERROR: GRN file not found: {grn_file}")
        return None
    
    # Detect file type by extension
    file_ext = os.path.splitext(grn_file)[1].lower()
    
    if file_ext == '.h5ad':
        # Load from H5AD (GRnnData format)
        print(f"Detected H5AD file format")
        grn_df = load_grn_from_h5ad(grn_file)
    elif file_ext == '.csv':
        # Load from CSV
        print(f"Detected CSV file format")
        grn_df = pd.read_csv(grn_file)
        print(f"Loaded {len(grn_df)} regulatory relationships")
        print(f"Columns: {grn_df.columns.tolist()}")
        print(f"\nFirst few edges:")
        print(grn_df.head())
    else:
        print(f"ERROR: Unsupported file format: {file_ext}")
        print("Supported formats: .csv, .h5ad")
        return None
    
    if grn_df is None:
        return None
    
    # Validate required columns
    required_cols = ['TF', 'Target']
    if not all(col in grn_df.columns for col in required_cols):
        print(f"ERROR: File must contain columns: {required_cols}")
        print(f"Found columns: {grn_df.columns.tolist()}")
        return None
    
    # Check for Correlation column, if not present, create it from weight or set to 1.0
    if 'Correlation' not in grn_df.columns:
        if 'Weight' in grn_df.columns:
            grn_df['Correlation'] = grn_df['Weight']
            print("Note: Using 'Weight' column as 'Correlation'")
        elif 'weight' in grn_df.columns:
            grn_df['Correlation'] = grn_df['weight']
            print("Note: Using 'weight' column as 'Correlation'")
        else:
            # If no correlation/weight column, set all to 1.0
            grn_df['Correlation'] = 1.0
            print("Note: No correlation/weight column found. Setting all correlations to 1.0")
    
    # Filter out edges with correlations below 0.0001
    MIN_CORRELATION_THRESHOLD = 0.0001
    initial_count = len(grn_df)
    grn_df = grn_df[grn_df['Correlation'].abs() >= MIN_CORRELATION_THRESHOLD].copy()
    filtered_count = len(grn_df)
    removed_count = initial_count - filtered_count
    
    if removed_count > 0:
        print(f"\nFiltered out {removed_count} edges with |correlation| < {MIN_CORRELATION_THRESHOLD}")
        print(f"Remaining edges: {filtered_count} (from {initial_count} total)")
    else:
        print(f"\nAll edges have |correlation| >= {MIN_CORRELATION_THRESHOLD}")
    
    return grn_df


def create_target_centered_subnetwork(grn_df, target_gene, max_tfs=MAX_TFS):
    """
    Create a subnetwork centered on a target gene.
    
    Args:
        grn_df: DataFrame with GRN edges (TF, Target, Correlation)
        target_gene: Name of the target gene to center the subnetwork on
        max_tfs: Maximum number of TFs to include (default: 10)
    
    Returns:
        DataFrame with edges for the subnetwork, target gene name, and set of TFs
    """
    print(f"\n{'='*60}")
    print(f"Creating subnetwork centered on target gene: {target_gene}")
    print(f"{'='*60}")
    
    # Find all TFs that regulate the target gene
    # Case-insensitive matching
    target_edges = grn_df[
        grn_df['Target'].str.upper() == target_gene.upper()
    ].copy()
    
    if len(target_edges) == 0:
        print(f"ERROR: No TFs found that regulate '{target_gene}'")
        # Show sample of available target genes
        available_targets = grn_df['Target'].unique()
        print(f"\nAvailable target genes (sample of first 20):")
        for i, gene in enumerate(available_targets[:20]):
            print(f"  {i+1}. {gene}")
        if len(available_targets) > 20:
            print(f"  ... and {len(available_targets) - 20} more")
        return None, None, None
    
    print(f"Found {len(target_edges)} TFs regulating {target_gene}")
    
    # Sort by absolute correlation and select top N TFs
    target_edges['AbsCorrelation'] = target_edges['Correlation'].abs()
    target_edges = target_edges.sort_values('AbsCorrelation', ascending=False)
    
    # Select top TFs
    n_tfs_to_select = min(max_tfs, len(target_edges))
    top_tfs_edges = target_edges.head(n_tfs_to_select)
    top_tfs = set(top_tfs_edges['TF'].unique())
    
    print(f"\nSelected top {len(top_tfs)} TFs (by correlation score):")
    for idx, row in top_tfs_edges.iterrows():
        print(f"  - {row['TF']} -> {target_gene} (correlation: {row['Correlation']:.4f})")
    
    # Now find all connections between these TFs
    print(f"\nFinding connections between the {len(top_tfs)} selected TFs...")
    
    # Get all edges where both source and target are in the TF set
    tf_connections = grn_df[
        (grn_df['TF'].isin(top_tfs)) & 
        (grn_df['Target'].isin(top_tfs))
    ].copy()
    
    print(f"Found {len(tf_connections)} connections between TFs")
    
    # Show sample of TF-TF connections
    if len(tf_connections) > 0:
        print(f"\nSample TF-TF connections (top 5 by correlation):")
        tf_connections_sorted = tf_connections.sort_values('Correlation', key=abs, ascending=False)
        for idx, row in tf_connections_sorted.head(5).iterrows():
            print(f"  - {row['TF']} -> {row['Target']} (correlation: {row['Correlation']:.4f})")
    
    # Combine: TF->Target edges + TF->TF edges
    subnetwork_edges = pd.concat([
        top_tfs_edges[['TF', 'Target', 'Correlation']],
        tf_connections[['TF', 'Target', 'Correlation']]
    ], ignore_index=True)
    
    # Remove duplicates (in case there are any)
    subnetwork_edges = subnetwork_edges.drop_duplicates(subset=['TF', 'Target'])
    
    print(f"\nSubnetwork summary:")
    print(f"  - Target gene: {target_gene}")
    print(f"  - Number of TFs: {len(top_tfs)}")
    print(f"  - Total edges: {len(subnetwork_edges)}")
    print(f"    - TF -> Target edges: {len(top_tfs_edges)}")
    print(f"    - TF -> TF edges: {len(tf_connections)}")
    
    return subnetwork_edges, target_gene, top_tfs


def create_networkx_graph(subnetwork_edges, target_gene, top_tfs):
    """Create a NetworkX graph from subnetwork edges"""
    print(f"\nCreating NetworkX graph from subnetwork...")
    G = nx.DiGraph()
    
    # Ensure target gene is added as a node (even if no edges point to it)
    G.add_node(target_gene)
    
    # Add edges
    for _, row in subnetwork_edges.iterrows():
        tf = str(row['TF'])
        target = str(row['Target'])
        corr = row['Correlation']
        
        G.add_edge(tf, target,
                  correlation=corr,
                  abs_correlation=abs(corr),
                  weight=abs(corr))
    
    # Identify node types
    all_nodes = set(G.nodes())
    tfs_in_subnet = top_tfs
    
    # Add node attributes
    for node in G.nodes():
        if node.upper() == target_gene.upper():
            G.nodes[node]['type'] = 'Target'
        elif node in tfs_in_subnet:
            G.nodes[node]['type'] = 'TF'
        else:
            G.nodes[node]['type'] = 'Other'
    
    # Create target_node set for visualization
    target_node = {node for node in G.nodes() if G.nodes[node]['type'] == 'Target'}
    
    print(f"Graph created:")
    print(f"  - Nodes: {G.number_of_nodes()}")
    print(f"  - Edges: {G.number_of_edges()}")
    print(f"  - TFs: {len(tfs_in_subnet)}")
    print(f"  - Target: {len(target_node)}")
    
    return G, tfs_in_subnet, target_node


def visualize_with_dash_cytoscape(G, tfs, target_node, target_gene, max_tfs):
    """Create an interactive web-based visualization with Dash Cytoscape"""
    try:
        import dash
        from dash import html, dcc
        import dash_cytoscape as cyto
        from dash.dependencies import Input, Output
        from dash.exceptions import PreventUpdate
    except ImportError:
        print("Error: dash and dash-cytoscape are not installed.")
        print("Install them with: pip install dash dash-cytoscape")
        return False
    
    cyto.load_extra_layouts()
    
    # Get edge weights for normalization
    weights = [attrs.get('weight', 1.0) for _, _, attrs in G.edges(data=True)]
    if weights:
        min_weight = min(weights)
        max_weight = max(weights)
        weight_range = max_weight - min_weight if max_weight > min_weight else 1.0
    else:
        min_weight = 0
        weight_range = 1.0
    
    # Build elements list for Cytoscape
    elements = []
    
    # Add nodes
    for node, attrs in G.nodes(data=True):
        node_type = attrs.get('type', 'Other')
        elements.append({
            'data': {
                'id': str(node),
                'label': str(node),
                'type': node_type
            }
        })
    
    # Add edges
    for source, target, attrs in G.edges(data=True):
        weight = attrs.get('weight', 1.0)
        correlation = attrs.get('correlation', 0.0)
        # Normalize weight to 0-1 range
        norm_weight = (weight - min_weight) / weight_range if weight_range > 0 else 0.5
        
        elements.append({
            'data': {
                'source': str(source),
                'target': str(target),
                'weight': weight,
                'correlation': correlation,
                'normalized_weight': norm_weight
            }
        })
    
    app = dash.Dash(__name__)
    
    app.css.append_css({
        'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
    })
    
    css_styles = {
        'download-button': {
            'background-color': '#4CAF50',
            'border': 'none',
            'color': 'white',
            'padding': '10px 20px',
            'text-align': 'center',
            'text-decoration': 'none',
            'display': 'inline-block',
            'font-size': '16px',
            'margin': '4px 2px',
            'cursor': 'pointer',
            'border-radius': '4px'
        },
        'download-button:hover': {
            'background-color': '#45a049'
        }
    }
    
    title = f"Target-Centered Subnetwork: {target_gene} and Top {len(tfs)} TFs"
    
    app.layout = html.Div([
        html.H1(title),
        
        html.Div([
            html.Div([
                html.Label("Select Layout:", style={'fontWeight': 'bold', 'marginRight': '10px'}),
                dcc.Dropdown(
                    id='layout-dropdown',
                    options=[
                        {'label': 'Concentric', 'value': 'concentric'},
                        {'label': 'Cose (Force-directed)', 'value': 'cose'},
                        {'label': 'Circle', 'value': 'circle'},
                        {'label': 'Grid', 'value': 'grid'},
                        {'label': 'Breadthfirst', 'value': 'breadthfirst'},
                        {'label': 'Cose-Bilkent (Best for large networks)', 'value': 'cose-bilkent'},
                        {'label': 'Euler', 'value': 'euler'},
                        {'label': 'Spread', 'value': 'spread'},
                        {'label': 'Dagre (Hierarchical)', 'value': 'dagre'},
                        {'label': 'Klay (Hierarchical)', 'value': 'klay'}
                    ],
                    value='concentric',  # Good for target-centered networks
                    style={'width': '300px'}
                ),
            ], style={'width': '300px', 'display': 'inline-block', 'margin': '10px'}),
            
            html.Div([
                html.Button("Export as PNG", id="btn-download-png", 
                           style=css_styles['download-button']),
                dcc.Download(id="download-png"),
                html.Button("Export as SVG", id="btn-download-svg", 
                           style=css_styles['download-button']),
                dcc.Download(id="download-svg"),
            ], style={'display': 'inline-block', 'margin': '10px', 'vertical-align': 'bottom'})
        ], style={'display': 'flex', 'flex-wrap': 'wrap', 'align-items': 'flex-end'}),
        
        html.Div(id='node-info', style={'margin': '10px', 'padding': '10px', 'backgroundColor': '#f0f0f0'}),
        
        # The network visualization
        cyto.Cytoscape(
            id='regulon-network',
            layout={'name': 'concentric'},
            style={'width': '100%', 'height': '800px'},
            elements=elements,
            stylesheet=[
                {
                    'selector': 'node',
                    'style': {
                        'content': 'data(label)',
                        'font-size': '14px',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'color': '#333333',
                        'text-outline-width': 2,
                        'text-outline-color': '#FFFFFF',
                        'text-outline-opacity': 1,
                        'background-opacity': 0,
                        'border-width': 0,
                    }
                },
                {
                    'selector': 'node[type = "Target"]',
                    'style': {
                        'shape': 'star',
                        'width': 100,
                        'height': 100,
                        'font-size': '18px',
                        'font-weight': 'bold',
                        'color': '#FFFFFF',
                        'background-color': '#FF6B6B',
                        'background-opacity': 1,
                        'border-color': '#CC0000',
                        'border-width': 4,
                        'border-opacity': 1,
                        'text-outline-color': '#000000',
                        'text-outline-width': 3,
                        'text-outline-opacity': 1,
                        'z-index': 30
                    }
                },
                {
                    'selector': 'node[type = "TF"]',
                    'style': {
                        'shape': 'diamond',
                        'width': 70,
                        'height': 70,
                        'background-color': '#FFFF00',
                        'background-opacity': 1,
                        'border-color': '#FFD700',
                        'border-width': 3,
                        'border-opacity': 1,
                        'font-weight': 'bold',
                        'font-size': '14px',
                        'color': '#000000',
                        'text-outline-color': '#FFFFFF',
                        'text-outline-width': 2
                    }
                },
                {
                    'selector': 'edge',
                    'style': {
                        'curve-style': 'bezier',
                        'target-arrow-shape': 'triangle',
                        'line-color': '#666666',
                        'target-arrow-color': '#666666',
                        'opacity': 0.8,
                        'width': 'mapData(normalized_weight, 0, 1, 1, 10)',
                        'target-arrow-width': 'mapData(normalized_weight, 0, 1, 3, 8)',
                        'target-arrow-height': 'mapData(normalized_weight, 0, 1, 3, 8)',
                    }
                },
                {
                    'selector': 'edge[normalized_weight >= 0.7]',
                    'style': {
                        'line-color': '#333333',
                        'target-arrow-color': '#333333',
                        'opacity': 1,
                        'z-index': 20  # Show strong connections on top
                    }
                }
            ]
        ),
        
        # Hidden div for storing graph data as JSON for export
        html.Div(id='graph-data-json', style={'display': 'none'})
    ])
    
    @app.callback(
        Output('node-info', 'children'),
        Input('regulon-network', 'tapNodeData')
    )
    def display_node_info(data):
        if not data:
            return f"Click on a node to see details. Target gene: {target_gene}"
        
        node_id = data['id']
        node_type = data['type']
        
        if node_type == 'Target':
            regulators = [source for source, _ in G.in_edges(node_id)]
            regulator_count = len(regulators)
            message = f"Target Gene: {node_id}\nRegulated by {regulator_count} transcription factors"
            if regulator_count > 0:
                message += f"\nRegulators: {', '.join(regulators)}"
        elif node_type == 'TF':
            targets = [target for _, target in G.out_edges(node_id)]
            target_count = len(targets)
            message = f"Transcription Factor: {node_id}\n"
            if target_gene.upper() in [t.upper() for t in targets]:
                message += f"✓ Regulates target gene: {target_gene}\n"
            message += f"Total targets in subnetwork: {target_count}"
            if target_count > 0:
                message += f"\nTargets: {', '.join(targets)}"
        else:
            message = f"Node: {node_id}\nType: {node_type}"
        
        return html.Pre(message)
    
    @app.callback(
        Output('regulon-network', 'layout'),
        Input('layout-dropdown', 'value')
    )
    def update_layout(layout_value):
        return {'name': layout_value}
    
    @app.callback(
        Output("download-png", "data"),
        Input("btn-download-png", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_png(n_clicks):
        if n_clicks is None:
            raise PreventUpdate
        return {
            "content": "This file will be replaced by the actual PNG export from the browser",
            "filename": f"subnetwork_{target_gene}.png",
            "type": "text/plain",
            "base64": False,
        }
    
    @app.callback(
        Output("download-svg", "data"),
        Input("btn-download-svg", "n_clicks"),
        prevent_initial_call=True,
    )
    def download_svg(n_clicks):
        if n_clicks is None:
            raise PreventUpdate
        return {
            "content": "This file will be replaced by the actual SVG export from the browser",
            "filename": f"subnetwork_{target_gene}.svg",
            "type": "text/plain",
            "base64": False,
        }
    
    app.clientside_callback(
        """
        function(n_clicks) {
            if (n_clicks === undefined) return;
            const cy = document.getElementById('regulon-network')._cyreg.cy;
            const png64 = cy.png({output: 'base64uri', scale: 2, bg: '#FFFFFF'});
            const link = document.createElement('a');
            link.href = png64;
            link.download = 'subnetwork.png';
            link.click();
            
            return;
        }
        """,
        Output("btn-download-png", "n_clicks"),
        Input("btn-download-png", "n_clicks"),
        prevent_initial_call=True,
    )
    
    app.clientside_callback(
        """
        function(n_clicks) {
            if (n_clicks === undefined) return;
            const cy = document.getElementById('regulon-network')._cyreg.cy;
            const svgContent = cy.svg({scale: 2, full: true, bg: '#FFFFFF'});
            const blob = new Blob([svgContent], {type: 'image/svg+xml'});
            const url = URL.createObjectURL(blob);
            const link = document.createElement('a');
            link.href = url;
            link.download = 'subnetwork.svg';
            link.click();
            URL.revokeObjectURL(url);
            
            return;
        }
        """,
        Output("btn-download-svg", "n_clicks"),
        Input("btn-download-svg", "n_clicks"),
        prevent_initial_call=True,
    )
    
    print("Starting Dash Cytoscape visualization...")
    print("Open your web browser at http://127.0.0.1:8050/ to view the network")
    print("\nFeatures:")
    print(f"- Target gene: {target_gene} (shown as red star)")
    print(f"- Top {len(tfs)} TFs (shown as yellow diamonds)")
    print("- Click on nodes to see details")
    print("- Use the layout dropdown to change the network layout")
    print("- Use the 'Export as PNG' or 'Export as SVG' buttons to save the visualization")
    
    app.run(debug=True)
    return True


def main():
    """Main function to run the visualization"""
    parser = argparse.ArgumentParser(
        description='Visualize target-centered subnetwork from GRN',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Visualize subnetwork for PLIN1 from CSV file:
  python visualize_target_subnetwork.py --grn-file pyscenic_output/grn_heart_failure_pyscenic.csv --target-gene PLIN1
  
  # Visualize subnetwork from H5AD file with custom number of TFs:
  python visualize_target_subnetwork.py --grn-file heart_failure_grn_scprint.h5ad --target-gene PLIN1 --max-tfs 15
  
  # Visualize subnetwork for a different target:
  python visualize_target_subnetwork.py --grn-file grnndata_output/grn_heart_failure_grnndata.h5ad --target-gene VWF
        """
    )
    parser.add_argument('--grn-file', default=DEFAULT_GRN_FILE,
                        help=f'Path to the GRN file (CSV or H5AD format). Default: {DEFAULT_GRN_FILE}')
    parser.add_argument('--target-gene', required=True,
                        help='Target gene name to create a centered subnetwork (required)')
    parser.add_argument('--max-tfs', type=int, default=MAX_TFS,
                        help=f'Maximum number of TFs to include (default: {MAX_TFS})')
    
    args = parser.parse_args()
    
    # Load GRN data (supports both CSV and H5AD)
    print("="*60)
    print("Target-Centered Subnetwork Visualization")
    print("="*60)
    print(f"Target Gene: {args.target_gene}")
    print(f"Max TFs: {args.max_tfs}")
    print("="*60)
    
    grn_df = load_grn_data(args.grn_file)
    if grn_df is None:
        print("\nERROR: Failed to load GRN data. Exiting.")
        return
    
    # Create target-centered subnetwork
    result = create_target_centered_subnetwork(grn_df, args.target_gene, args.max_tfs)
    if result[0] is None:
        print("\nERROR: Failed to create target-centered subnetwork. Exiting.")
        return
    
    subnetwork_edges, target_gene, top_tfs = result
    
    # Create NetworkX graph
    G, tfs, target_node = create_networkx_graph(subnetwork_edges, target_gene, top_tfs)
    
    # Visualize with Dash Cytoscape
    visualize_with_dash_cytoscape(G, tfs, target_node, target_gene, args.max_tfs)


if __name__ == "__main__":
    main()

