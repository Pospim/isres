#!/usr/bin/env python3
"""
Create an interactive HTML file mapping motifs to their locations in sequences.
Main view shows unique motifs with occurrence counts.
Detail view shows all occurrences of a selected motif.
"""

import csv
from collections import defaultdict
import html
import re
import argparse
import pandas as pd
import plotly.express as px
from pathlib import Path


def extract_gene_name(header, id_format ):
    """
    Extract gene name from FASTA header based on format.

    Args:
        header: FASTA header string (with or without '>')
        id_format: One of 'ensembl', 'gene_ensembl', or 'gene'

    Returns:
        Extracted gene name
    """
    header = header.lstrip('>')

    if id_format == 'ensembl':
        # Format: ENSGALG00010001364 upstream 2000bp
        # Extract just the ENSEMBL ID
        match = re.match(r'(ENS\w+)', header)
        if match:
            return match.group(1)
        return header.split()[0]

    elif id_format == 'gene_ensembl':
        # Format: ZNFX1 (ENSG00000124201) upstream 2000bp
        # Extract the gene name before the parenthesis
        match = re.match(r'([^\s(]+)', header)
        if match:
            return match.group(1)
        return header.split()[0]

    elif id_format == 'gene':
        # Format: ZNFX1 or ZNFX1 upstream 2000bp
        # Extract just the gene name
        return header.split()[0]

    else:
        raise ValueError(f"Unknown id_format: {id_format}")


def read_fasta(fasta_file, id_format='gene_ensembl'):
    """
    Read FASTA file and return dictionary of {gene_name: (header, sequence)}

    Args:
        fasta_file: Path to FASTA file
        id_format: Format of FASTA headers - 'ensembl', 'gene_ensembl', or 'gene'
    """
    sequences = {}
    current_header = None
    current_gene = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene:
                    sequences[current_gene] = (current_header, ''.join(current_seq))
                current_header = line
                current_gene = extract_gene_name(line, id_format)
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_gene:
            sequences[current_gene] = (current_header, ''.join(current_seq))

    print(f"  Sample FASTA keys: {list(sequences.keys())[:3]}")
    return sequences


def read_motif_mappings(tsv_file, id_format='gene_ensembl'):
    """
    Read TSV file and return list of motif mappings

    Args:
        tsv_file: Path to TSV file with motif hits
        id_format: Format to extract gene names - 'ensembl', 'gene_ensembl', or 'gene'
    """
    mappings = []

    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Determine which column to use for FASTA key matching
            if id_format == 'ensembl' and 'ensembl_id' in row and row['ensembl_id']:
                # For ensembl format, use ensembl_id for matching FASTA keys
                fasta_key = row['ensembl_id']
            elif 'gene_name' in row and row['gene_name'] and id_format != 'ensembl':
                # For other formats, use gene_name if available
                fasta_key = row['gene_name']
            elif 'FASTA ID' in row:
                fasta_key = extract_gene_name(row['FASTA ID'], id_format)
            elif 'gene' in row:
                fasta_key = extract_gene_name(row['gene'], id_format)
            else:
                raise ValueError("TSV file must have either 'gene_name', 'ensembl_id', 'FASTA ID', or 'gene' column")

            # Store the FASTA key for sequence lookup
            row['gene'] = fasta_key
            
            # Store display name separately if available (prefer gene_name over ensembl_id)
            if 'gene_name' in row and row['gene_name']:
                row['display_name'] = row['gene_name']
            else:
                row['display_name'] = fasta_key
            
            mappings.append(row)

    print(f"  Sample mapping genes (FASTA keys): {list(set([m['gene'] for m in mappings[:100]]))[:3]}")
    print(f"  Sample display names: {list(set([m['display_name'] for m in mappings[:100]]))[:3]}")
    return mappings


def generate_heatmap(mappings, output_dir, top_n_motifs=20):
    """
    Generate an interactive plotly heatmap showing motif occurrences across genes.

    Args:
        mappings: List of motif mappings
        output_dir: Directory to save the heatmap HTML file
        top_n_motifs: Number of top motifs to display

    Returns:
        Path to the generated heatmap HTML file
    """
    # Convert mappings to DataFrame
    df = pd.DataFrame(mappings)
    
    # Use display_name for heatmap y-axis labels (human-readable names)
    df['gene_display'] = df['display_name']

    # Create a pivot table: rows = gene_display, columns = motif, values = count
    heatmap_data_full = df.groupby(['gene_display', 'motif']).size().unstack(fill_value=0)

    # Get the top N most frequent motifs (by total occurrences across all genes)
    motif_totals = heatmap_data_full.sum(axis=0).sort_values(ascending=False)
    top_motifs = motif_totals.head(top_n_motifs).index

    # Filter the heatmap to show only top motifs
    heatmap_data = heatmap_data_full[top_motifs]

    print(f"\nHeatmap Statistics:")
    print(f"  - Full dataset: {heatmap_data_full.shape[0]} genes × {heatmap_data_full.shape[1]} motifs")
    print(f"  - Filtered to top {top_n_motifs} motifs")
    print(f"  - Heatmap dimensions: {heatmap_data.shape[0]} genes × {heatmap_data.shape[1]} motifs")
    print(f"\n  Top {min(10, top_n_motifs)} most frequent motifs:")
    for i, (motif, count) in enumerate(motif_totals.head(10).items(), 1):
        print(f"    {i}. {motif}: {count} occurrences")

    # Create interactive heatmap
    fig = px.imshow(heatmap_data,
                    labels=dict(x="Motif", y="Gene Name", color="Occurrence Count"),
                    x=heatmap_data.columns,
                    y=heatmap_data.index,
                    color_continuous_scale='YlOrRd',
                    aspect='auto',
                    title=f'Interactive Gene-Motif Occurrence Heatmap (Top {top_n_motifs} Motifs)')

    fig.update_xaxes(side="bottom", tickangle=90)
    fig.update_layout(
        width=min(1600, max(1000, top_n_motifs * 20)),
        height=800,
        xaxis={'tickfont': {'size': 8}},
        yaxis={'tickfont': {'size': 6}}
    )

    # Save heatmap
    output_path = Path(output_dir)
    heatmap_file = output_path.parent / f'{output_path.stem}_heatmap_top{top_n_motifs}.html'
    fig.write_html(str(heatmap_file))

    print(f"\n✓ Interactive heatmap saved to: {heatmap_file}")
    print(f"  (Showing top {top_n_motifs} most frequent motifs out of {heatmap_data_full.shape[1]} total)")

    return str(heatmap_file)


def generate_html(sequences, mappings, output_file, heatmap_file=None):
    """Generate interactive HTML file with main view and detail view"""

    # Count motif occurrences
    motif_counts = defaultdict(int)
    for mapping in mappings:
        motif_base = mapping['motif'].split('-', 1)[1] if '-' in mapping['motif'] else mapping['motif']
        motif_counts[motif_base] += 1

    # Sort motifs by count (descending)
    sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)

    # Group mappings by motif and gene
    motif_mappings = defaultdict(list)
    gene_motifs_map = defaultdict(set)

    for mapping in mappings:
        motif_base = mapping['motif'].split('-', 1)[1] if '-' in mapping['motif'] else mapping['motif']
        motif_mappings[motif_base].append(mapping)
        gene_motifs_map[mapping['gene']].add(motif_base)

    # Start building HTML
    html_content = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Motif Mapper</title>
    <style>
        body {
            font-family: 'Courier New', monospace;
            margin: 20px;
            background-color: #f5f5f5;
        }

        .container {
            display: flex;
            gap: 20px;
        }

        #mainView .container {
            margin-left: 0;
        }

        .motif-list {
            width: 400px;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow-y: auto;
            max-height: 90vh;
        }

        #mainView .motif-list {
            position: sticky;
            top: 20px;
            height: fit-content;
        }

        #detailView .motif-list {
            position: fixed;
            left: 0;
            top: 80px;
            height: calc(100vh - 100px);
            transition: transform 0.3s ease;
            z-index: 1000;
            transform: translateX(-380px);
        }

        #detailView .motif-list:hover {
            transform: translateX(0);
        }

        .motif-list-trigger {
            position: fixed;
            left: 0;
            top: 80px;
            width: 20px;
            height: calc(100vh - 100px);
            background: rgba(76, 175, 80, 0.3);
            cursor: pointer;
            z-index: 999;
            display: none;
        }

        #detailView .motif-list-trigger {
            display: block;
        }

        .sequences {
            flex: 1;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        #mainView .sequences {
            margin-right: 20px;
            align-self: flex-start;
            position: sticky;
            top: 20px;
            max-height: 90vh;
            overflow-y: auto;
        }

        #detailView .sequences {
            margin-left: 40px;
            margin-bottom: 40px;
        }

        h1 {
            color: #333;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 10px;
        }

        h2 {
            color: #555;
            margin-top: 30px;
            font-size: 16px;
        }

        .motif-item {
            padding: 8px;
            margin: 5px 0;
            background: #f9f9f9;
            border-left: 3px solid #4CAF50;
            cursor: pointer;
            border-radius: 4px;
            transition: all 0.2s;
        }

        .motif-item:hover {
            background: #e8f5e9;
            transform: translateX(5px);
        }

        .motif-name {
            font-weight: bold;
            color: #2196F3;
        }

        .motif-count {
            background: #2196F3;
            color: white;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: bold;
            margin-left: 8px;
        }

        .motif-info {
            font-size: 12px;
            color: #666;
            margin-top: 3px;
        }

        .gene-section {
            margin-bottom: 40px;
            padding: 15px;
            background: #fafafa;
            border-radius: 8px;
        }

        .gene-header {
            font-weight: bold;
            color: #333;
            margin-bottom: 15px;
            padding: 10px;
            background: #e3f2fd;
            border-radius: 4px;
        }

        .sequence {
            word-wrap: break-word;
            line-height: 1.6;
            font-size: 11px;
            background: white;
            padding: 15px;
            border-radius: 4px;
            border: 1px solid #ddd;
            font-family: 'Courier New', monospace;
            overflow-x: auto;
        }

        .highlight {
            background-color: #ffeb3b;
            font-weight: bold;
            padding: 2px 4px;
            border-radius: 3px;
            animation: pulse 1s ease-in-out;
        }

        @keyframes pulse {
            0%, 100% { background-color: #ffeb3b; }
            50% { background-color: #ffc107; }
        }

        .strand-indicator {
            display: inline-block;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 11px;
            margin-left: 5px;
        }

        .strand-plus {
            background: #c8e6c9;
            color: #2e7d32;
        }

        .strand-minus {
            background: #ffccbc;
            color: #d84315;
        }

        .score {
            color: #9c27b0;
            font-weight: bold;
        }

        .search-box {
            margin-bottom: 15px;
            width: 100%;
            padding: 10px;
            border: 2px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }

        .search-box:focus {
            outline: none;
            border-color: #4CAF50;
        }

        .stats {
            background: #e8f5e9;
            padding: 10px;
            border-radius: 4px;
            margin-bottom: 15px;
            font-size: 12px;
        }

        .back-button {
            background: #4CAF50;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            margin-bottom: 15px;
            transition: background 0.2s;
            width: 100%;
        }

        .back-button:hover {
            background: #45a049;
        }

        .heatmap-button {
            background: #2196F3;
            color: white;
            border: none;
            padding: 10px 20px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            margin-bottom: 15px;
            transition: background 0.2s;
            width: 100%;
            text-align: center;
            text-decoration: none;
            display: inline-block;
        }

        .heatmap-button:hover {
            background: #1976D2;
        }

        #detailView {
            display: none;
            position: absolute;
            top: 0;
            left: 0;
            width: 100%;
        }

        #mainView {
            display: block;
            position: relative;
        }

        .occurrence-item {
            padding: 8px;
            margin: 5px 0;
            background: #f9f9f9;
            border-left: 3px solid #FF9800;
            cursor: pointer;
            border-radius: 4px;
            transition: all 0.2s;
        }

        .occurrence-item:hover {
            background: #fff3e0;
            transform: translateX(5px);
        }
    </style>
    <script>
        function showMotifDetails(motif) {
            console.log('Showing details for motif:', motif);

            // Hide main view and show detail view
            document.getElementById('mainView').style.display = 'none';
            document.getElementById('detailView').style.display = 'block';
            document.getElementById('currentMotif').textContent = motif;

            // Show only occurrences of this motif
            const allOccurrences = document.querySelectorAll('.occurrence-item');
            let occCount = 0;
            allOccurrences.forEach(item => {
                if (item.dataset.motif === motif) {
                    item.style.display = 'block';
                    occCount++;
                } else {
                    item.style.display = 'none';
                }
            });
            console.log('Found', occCount, 'occurrences');

            // Show only gene sections with this motif
            const allSections = document.querySelectorAll('.gene-section');
            let sectionCount = 0;
            allSections.forEach(section => {
                const motifs = section.dataset.motifs ? section.dataset.motifs.split(',') : [];
                if (motifs.includes(motif)) {
                    section.style.display = 'block';
                    sectionCount++;
                } else {
                    section.style.display = 'none';
                }
            });
            console.log('Found', sectionCount, 'gene sections with this motif');

            // Scroll to top when showing details
            window.scrollTo(0, 0);
        }

        function backToMain() {
            const detailView = document.getElementById('detailView');
            const mainView = document.getElementById('mainView');

            // Hide detail view
            detailView.style.display = 'none';

            // Reset all gene sections to hidden (they should only show in detail view)
            const allSections = document.querySelectorAll('.gene-section');
            allSections.forEach(section => {
                section.style.display = 'none';
            });

            // Show main view
            mainView.style.display = 'block';

            // Scroll to top when returning to main view
            setTimeout(() => window.scrollTo(0, 0), 0);
        }

        function scrollToMotif(geneId, start, end) {
            // Remove existing highlights
            document.querySelectorAll('.highlight').forEach(el => {
                el.classList.remove('highlight');
            });

            // Add highlight to target
            const targetId = geneId + '_' + start + '_' + end;
            const element = document.getElementById(targetId);

            if (element) {
                element.classList.add('highlight');
                element.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
        }

        function filterMotifs() {
            const searchTerm = document.getElementById('searchBox').value.toLowerCase();
            const motifItems = document.querySelectorAll('.motif-item');

            motifItems.forEach(item => {
                const text = item.textContent.toLowerCase();
                if (text.includes(searchTerm)) {
                    item.style.display = 'block';
                } else {
                    item.style.display = 'none';
                }
            });
        }

        function filterOccurrences() {
            const searchTerm = document.getElementById('searchBoxDetail').value.toLowerCase();
            const occItems = document.querySelectorAll('.occurrence-item');

            occItems.forEach(item => {
                const text = item.textContent.toLowerCase();
                if (text.includes(searchTerm)) {
                    item.style.display = 'block';
                } else {
                    item.style.display = 'none';
                }
            });
        }
    </script>
</head>
<body>
    <h1>🧬 Motif Mapper - Interactive Visualization</h1>

    <!-- Main View: Unique Motifs -->
    <div id="mainView" class="container">
        <div class="sequences">
            <h2>About</h2>
            <p>This page shows all unique motifs found in 2000bp of ISG genes.</p>
            <p><strong>Click on a motif</strong> in the list on the right to see all its occurrences across genes and their sequence locations.</p>
            <div style="margin-top: 20px; padding: 15px; background: #e3f2fd; border-left: 4px solid #2196F3; border-radius: 4px;">
                <strong>📊 Statistics:</strong><br>
                Total Motif Hits: """ + str(len(mappings)) + """<br>
                Genes Analyzed: """ + str(len(sequences)) + """<br>
                Unique Motifs Found: """ + str(len(sorted_motifs)) + """
            </div>
"""

    # Add heatmap button if available
    if heatmap_file:
        heatmap_filename = Path(heatmap_file).name
        html_content += f"""
            <a href="{heatmap_filename}" target="_blank" class="heatmap-button" style="margin-top: 20px;">
                📊 View Gene-Motif Heatmap
            </a>
"""

    html_content += """
        </div>
        <div class="motif-list">
            <h2>Unique Motifs</h2>
            <input type="text" id="searchBox" class="search-box"
                   placeholder="Search motifs..."
                   onkeyup="filterMotifs()">

            <div class="stats">
                <strong>Total Hits:</strong> """ + str(len(mappings)) + """<br>
                <strong>Genes:</strong> """ + str(len(sequences)) + """<br>
                <strong>Unique Motifs:</strong> """ + str(len(sorted_motifs)) + """
            </div>
"""

    # Add motif list (sorted by count)
    for motif_base, count in sorted_motifs:
        # Escape for JavaScript string (replace " with \")
        motif_escaped_js = motif_base.replace('\\', '\\\\').replace('"', '\\"')
        html_content += f"""
            <div class="motif-item" onclick='showMotifDetails("{motif_escaped_js}")'>
                <div class="motif-name">
                    {html.escape(motif_base)}
                    <span class="motif-count">{count}×</span>
                </div>
                <div class="motif-info">
                    Click to view all {count} occurrences
                </div>
            </div>
"""

    html_content += """
        </div>
    </div>

    <!-- Detail View: Motif Occurrences -->
    <div id="detailView" class="container">
        <div class="motif-list-trigger"></div>
        <div class="motif-list">
            <button class="back-button" onclick="backToMain()">← Back to All Motifs</button>
            <h2>Occurrences of: <span id="currentMotif"></span></h2>
            <input type="text" id="searchBoxDetail" class="search-box"
                   placeholder="Search occurrences..."
                   onkeyup="filterOccurrences()">
            <div id="occurrencesList">
"""

    # Add all occurrence items (initially hidden)
    for mapping in mappings:
        motif_base = mapping['motif'].split('-', 1)[1] if '-' in mapping['motif'] else mapping['motif']
        gene = mapping['gene']  # FASTA key for sequence lookup
        display_name = mapping.get('display_name', gene)  # Human-readable name
        start = mapping['start']
        end = mapping['end']
        strand = mapping['strand']
        score = mapping['score']
        matched_seq = mapping['matched_sequence']

        # Use gene name directly as gene_id (matches FASTA key and gene section ID)
        gene_id = gene

        # Get full header for display
        if gene in sequences:
            full_header = sequences[gene][0]
        else:
            full_header = gene

        strand_class = 'strand-plus' if strand == '+' else 'strand-minus'

        # Escape motif for HTML attribute
        motif_escaped_attr = html.escape(motif_base, quote=False)

        html_content += f"""
            <div class="occurrence-item" data-motif="{motif_escaped_attr}"
                 onclick="scrollToMotif('{gene_id}', '{start}', '{end}')" style="display:none;">
                <div class="motif-info">
                    <strong>{html.escape(display_name)}</strong> ({html.escape(full_header.lstrip('>'))})<br>
                    Position: {start}-{end}
                    <span class="strand-indicator {strand_class}">{strand}</span><br>
                    Sequence: <span style="color: #f44336;">{html.escape(matched_seq)}</span><br>
                    Score: <span class="score">{score}</span>
                </div>
            </div>
"""

    html_content += """
            </div>
        </div>

        <div class="sequences">
            <h2>Sequences with Motif Locations</h2>
"""

    # Add gene sections with full sequences (initially hidden)
    for gene in sorted(sequences.keys()):
        header, sequence = sequences[gene]
        # Use gene name directly as ID (matches TSV gene_name and FASTA key)
        gene_id = gene

        # Get all motifs for this gene
        gene_motifs_list = sorted(gene_motifs_map.get(gene, []))
        motifs_str = ','.join(gene_motifs_list)
        motifs_str_escaped = html.escape(motifs_str, quote=False)

        # Get all mappings for this gene
        gene_mappings = [m for m in mappings if m['gene'] == gene]

        if not gene_mappings:
            continue

        # Get display name from the first mapping (they should all have the same display_name for this gene)
        display_name = gene_mappings[0].get('display_name', gene)

        # Sort mappings by position
        gene_mappings.sort(key=lambda x: int(x['start']))

        html_content += f"""
            <div class="gene-section" id="{gene_id}" data-motifs="{motifs_str_escaped}" style="display:none;">
                <div class="gene-header"><strong>{html.escape(display_name)}</strong> ({html.escape(header.lstrip('>'))})</div>
                <div class="sequence">
"""

        # Build sequence with inline highlighted motifs
        seq_with_highlights = []
        last_pos = 0

        for mapping in gene_mappings:
            start = int(mapping['start']) - 1  # Convert to 0-based
            end = int(mapping['end'])

            # Add sequence before this motif
            if start > last_pos:
                seq_with_highlights.append(html.escape(sequence[last_pos:start]))

            # Add the highlighted motif
            motif_seq = sequence[start:end]
            target_id = f"{gene_id}_{mapping['start']}_{mapping['end']}"
            seq_with_highlights.append(
                f'<span id="{target_id}" style="background-color: #fff9c4; border: 2px solid #fbc02d; padding: 2px;">{html.escape(motif_seq)}</span>'
            )

            last_pos = end

        # Add remaining sequence
        if last_pos < len(sequence):
            seq_with_highlights.append(html.escape(sequence[last_pos:]))

        html_content += ''.join(seq_with_highlights)

        html_content += """
                </div>
            </div>
"""

    html_content += """
        </div>
    </div>
</body>
</html>
"""

    # Write HTML file
    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"HTML file generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Create an interactive HTML file mapping motifs to their locations in sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
FASTA ID Format Options:
  gene_ensembl : Gene names with ENSEMBL IDs (default)
                 Example: >ZNFX1 (ENSG00000124201) upstream 2000bp

  ensembl      : ENSEMBL IDs only
                 Example: >ENSGALG00010001364 upstream 2000bp

  gene         : Gene names only
                 Example: >ZNFX1 upstream 2000bp

Examples:
  python map_motifs_to_html_v2.py motif_hits.tsv sequences.fa output.html
  python map_motifs_to_html_v2.py motif_hits.tsv sequences.fa output.html --format ensembl
  python map_motifs_to_html_v2.py motif_hits.tsv sequences.fa output.html --top_motifs 30
        """
    )

    parser.add_argument('tsv_file',
                        help='TSV file with motif hits')
    parser.add_argument('fasta_file',
                        help='FASTA file with sequences')
    parser.add_argument('output_file',
                        help='Output HTML file path')
    parser.add_argument('-f', '--format',
                        choices=['ensembl', 'gene_ensembl', 'gene'],
                        default='gene_ensembl',
                        help='Format of FASTA headers (default: gene_ensembl)')
    parser.add_argument('--top_motifs',
                        type=int,
                        default=20,
                        help='Number of top motifs to display in heatmap (default: 20)')

    args = parser.parse_args()

    print(f"Reading FASTA file (format: {args.format})...")
    sequences = read_fasta(args.fasta_file, args.format)
    print(f"  Found {len(sequences)} sequences")

    print("Reading motif mappings...")
    mappings = read_motif_mappings(args.tsv_file, args.format)
    print(f"  Found {len(mappings)} motif occurrences")

    # Generate heatmap
    print("\nGenerating heatmap...")
    heatmap_file = generate_heatmap(mappings, args.output_file, args.top_motifs)

    print("\nGenerating HTML...")
    generate_html(sequences, mappings, args.output_file, heatmap_file)

    print("\n✓ Done! Open the HTML file in your browser to view the interactive map.")
    print(f"✓ Heatmap available at: {heatmap_file}")


if __name__ == '__main__':
    main()
