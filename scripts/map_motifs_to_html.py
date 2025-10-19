#!/usr/bin/env python3
"""
Create an interactive HTML file mapping motifs to their locations in sequences.
Clicking on a motif will scroll to and highlight the corresponding sequence region.
"""

import csv
from collections import defaultdict
import html


def read_fasta(fasta_file):
    """Read FASTA file and return dictionary of {header: sequence}"""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


def read_motif_mappings(tsv_file):
    """Read TSV file and return list of motif mappings"""
    mappings = []

    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            mappings.append(row)

    return mappings


def generate_html(sequences, mappings, output_file):
    """Generate interactive HTML file"""

    # Group mappings by gene
    gene_mappings = defaultdict(list)
    for mapping in mappings:
        gene = mapping['gene']
        gene_mappings[gene].append(mapping)

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

        .motif-list {
            width: 400px;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            position: sticky;
            top: 20px;
            height: fit-content;
            max-height: 90vh;
            overflow-y: auto;
        }

        .sequences {
            flex: 1;
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
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
            line-height: 1.8;
            font-size: 14px;
            background: white;
            padding: 15px;
            border-radius: 4px;
            border: 1px solid #ddd;
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

        .position-marker {
            color: #999;
            font-size: 10px;
            display: inline-block;
            margin-right: 10px;
            min-width: 40px;
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
    </style>
    <script>
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
    </script>
</head>
<body>
    <h1>🧬 Motif Mapper - Interactive Visualization</h1>

    <div class="container">
        <div class="motif-list">
            <h2>Motif Occurrences</h2>
            <input type="text" id="searchBox" class="search-box"
                   placeholder="Search motifs, genes, or sequences..."
                   onkeyup="filterMotifs()">

            <div class="stats">
                <strong>Total Hits:</strong> """ + str(len(mappings)) + """<br>
                <strong>Genes:</strong> """ + str(len(gene_mappings)) + """<br>
                <strong>Unique Motifs:</strong> """ + str(len(set(m['motif'] for m in mappings))) + """
            </div>
"""

    # Add motif list
    for mapping in mappings:
        motif_base = mapping['motif'].split('-', 1)[1] if '-' in mapping['motif'] else mapping['motif']
        gene = mapping['gene']
        gene_id = gene.replace('>', '').replace(' ', '_')
        start = mapping['start']
        end = mapping['end']
        strand = mapping['strand']
        score = mapping['score']
        sequence = mapping['reverse'] if mapping['strand'] == '-' else motif_base

        strand_class = 'strand-plus' if strand == '+' else 'strand-minus'

        html_content += f"""
            <div class="motif-item" onclick="scrollToMotif('{gene_id}', '{start}', '{end}')">
                <div class="motif-name">{html.escape(motif_base)}</div>
                <div class="motif-info">
                    {html.escape(gene)}<br>
                    Position: {start}-{end}
                    <span class="strand-indicator {strand_class}">{strand}</span><br>
                    Sequence: <span style="color: #f44336;">{html.escape(sequence)}</span><br>
                    Score: <span class="score">{score}</span>
                </div>
            </div>
"""

    html_content += """
        </div>

        <div class="sequences">
            <h2>Sequences with Motif Locations</h2>
"""

    # Add sequences with highlighted motifs
    for gene, gene_motifs in sorted(gene_mappings.items()):
        if gene not in sequences:
            continue

        sequence = sequences[gene]
        gene_id = gene.replace('>', '').replace(' ', '_')

        html_content += f"""
            <div class="gene-section" id="{gene_id}">
                <div class="gene-header">{html.escape(gene)}</div>
                <div class="sequence">
"""

        # Create a list of (position, type, data) for sorting
        events = []
        for motif in gene_motifs:
            start = int(motif['start']) - 1  # Convert to 0-based
            end = int(motif['end'])
            events.append({
                'pos': start,
                'type': 'start',
                'end': end,
                'motif': motif
            })

        events.sort(key=lambda x: x['pos'])

        # Build sequence with highlights
        last_pos = 0

        for i, event in enumerate(events):
            start = event['pos']
            end = event['end']
            motif = event['motif']

            # Add sequence before this motif
            if start > last_pos:
                # Add position markers every 60 characters
                chunk = sequence[last_pos:start]
                for j in range(0, len(chunk), 60):
                    pos_marker = last_pos + j + 1  # 1-based position
                    html_content += f'<span class="position-marker">{pos_marker:4d}:</span>'
                    html_content += html.escape(chunk[j:j+60]) + '<br>'

            # Add the highlighted motif
            motif_seq = sequence[start:end]
            target_id = f"{gene_id}_{motif['start']}_{motif['end']}"

            # Position marker for motif
            html_content += f'<span class="position-marker">{start + 1:4d}:</span>'
            html_content += f'<span id="{target_id}" style="background-color: #fff9c4; border: 2px solid #fbc02d; padding: 2px;">'
            html_content += html.escape(motif_seq)
            html_content += '</span><br>'

            last_pos = end

        # Add remaining sequence
        if last_pos < len(sequence):
            chunk = sequence[last_pos:]
            for j in range(0, len(chunk), 60):
                pos_marker = last_pos + j + 1
                html_content += f'<span class="position-marker">{pos_marker:4d}:</span>'
                html_content += html.escape(chunk[j:j+60]) + '<br>'

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
    import sys

    if len(sys.argv) < 4:
        print("Usage: python map_motifs_to_html.py <tsv_file> <fasta_file> <output_html>")
        print("\nExample:")
        print("  python map_motifs_to_html.py motif_hits.tsv sequences.fa output.html")
        sys.exit(1)

    tsv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    print("Reading FASTA file...")
    sequences = read_fasta(fasta_file)
    print(f"  Found {len(sequences)} sequences")

    print("Reading motif mappings...")
    mappings = read_motif_mappings(tsv_file)
    print(f"  Found {len(mappings)} motif occurrences")

    print("Generating HTML...")
    generate_html(sequences, mappings, output_file)

    print("\n✓ Done! Open the HTML file in your browser to view the interactive map.")


if __name__ == '__main__':
    main()
