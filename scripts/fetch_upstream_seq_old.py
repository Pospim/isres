import requests
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

UPSTREAM_REGION = 2000
SPECIES = "homo_sapiens"

def gene_id_to_ensmbl(gene_id, species=SPECIES):
    """
    Convert gene ID to Ensembl ID using Ensembl REST API.
    """
    if gene_id.startswith("ENS"):
        print(f"Already Ensembl ID: {gene_id}")
        return gene_id

    url = f"https://rest.ensembl.org/xrefs/symbol/{species}/{gene_id}?content-type=application/json"
    response = requests.get(url)

    if not response.ok:
        print(f"Error mapping gene_id '{gene_id}': {response.text}")
        return None

    data = response.json()
    for entry in data:
        if entry["type"] == "gene":
            print(f"{gene_id} mapped to {entry['id']}")
            return entry["id"]  # Return the Ensembl Gene ID
    return None

def batch_gene_to_ensembl(df, species=SPECIES) -> pd.DataFrame:
    """
    Fetch Ensembl IDs for gene IDs in a DataFrame.
    Args:
        df (pd.DataFrame): Input DataFrame with 'gene_id' column.
        species (str): Species to query Ensembl IDs from.
    Returns:
        pd.DataFrame: Updated DataFrame with a new 'ensmbl_id' column.
    """
    df = df.drop_duplicates(subset=['gene_id']).reset_index(drop=True)
    df['ensmbl_id'] = df['gene_id'].apply(lambda x: gene_id_to_ensmbl(x, species))
    print("-------------------\nEnsembl IDs fetched for gene IDs.\n-------------------")
    return df

def get_gene_coordinates(ensmbl_id) -> dict:
    """
    Fetch gene coordinates and strand information for a given Ensembl Gene ID.
    """
    url = f"https://rest.ensembl.org/lookup/id/{ensmbl_id}?content-type=application/json"
    response = requests.get(url)

    if not response.ok:
        response.raise_for_status()

    data = response.json()
    print(f"Gene {ensmbl_id} coordinates fetched.")
    return {
        "chromosome": data["seq_region_name"],
        "gene_start": data["start"],
        "gene_end": data["end"],
        "strand": data["strand"]
    }

def fetch_upstream_seq(chromosome, gene_start, gene_end, strand, species=SPECIES, upstream=UPSTREAM_REGION) -> str:
    """
    Fetch sequence with upstream region from Ensembl REST API.
    """
    # Adjust coordinates for upstream region
    if strand == 1:
        region_start = max(1, gene_start - upstream)
        region_end = gene_start - 1
    else:
        region_start = gene_end + 1
        region_end = gene_end + upstream


    # Construct API URL
    seq_url = f"https://rest.ensembl.org/sequence/region/{species}/{chromosome}:{region_start}..{region_end}:{strand}?content-type=text/x-fasta"

    #headers = {"Content-Type": "text/x-fasta"}
    #seq_response = requests.get(seq_url, headers=headers)
    headers = {"Content-Type": "text/x-fasta", "User-Agent": "pospim.inq@gmail.com"}
    response = requests.get(seq_url, headers=headers)


    if not response.ok:
        raise Exception(f"[ERROR] Failed to fetch promoter sequence: {chromosome}:{region_start}-{region_end} (strand {strand})")

    fasta = response.text
    if not fasta.startswith(">"):
        raise ValueError(f"[ERROR] Invalid FASTA returned for {chromosome}:{region_start}-{region_end}")

    print(f"[INFO] Promoter region fetched: {chromosome}:{region_start}-{region_end} (strand {strand})")
    return "".join(fasta.split("\n")[1:])

def get_sequences_to_fasta(df, output_fasta_path, species=SPECIES, upstream=UPSTREAM_REGION):
    """
    Process a DataFrame with Ensembl gene IDs and save sequences to a FASTA file.
    Args:
        df (pd.DataFrame): Input DataFrame with 'ensmbl_id' column.
        output_fasta_path (str): Path to save the FASTA file.
        species (str): Species to query sequences from.
        upstream (int): Number of nucleotides upstream to include.
    """
    # Ensure no duplicates in the Ensembl ID column
    df = df.drop_duplicates(subset=['ensmbl_id']).reset_index(drop=True)

    records = []
    for _, row in df.iterrows():
        try:
            # Fetch gene coordinates
            print(f"[INFO] Fetching coordinates for {row['ensmbl_id']}")
            coords = get_gene_coordinates(row['ensmbl_id'])
            # Fetch upstream sequence
            seq = fetch_upstream_seq(**coords, species=species, upstream=upstream)
            sanitized_seq = seq.replace("\n", "").replace(" ", "")
            # Create a SeqRecord for FASTA
            record = SeqRecord(
                Seq(sanitized_seq),
                id=row['ensmbl_id'],  # Use Ensembl ID as the record ID
                description=f"{row['ensmbl_id']} upstream {upstream}bp"
            )
            records.append(record)
            print(f"[INFO] Added record for {row['ensmbl_id']}")
            print("---------------------------------------------------------------")

        except Exception as e:
            print(f"[ERROR] Skipping {row['ensmbl_id']}: {type(e).__name__} – {e}")

    SeqIO.write(records, output_fasta_path, "fasta")
    print(f"[INFO] FASTA saved to: {output_fasta_path}")



df_path = "/home/pospim/Desktop/school/isres/human/human_isgs.csv"
out_fasta = "/home/pospim/Desktop/school/isres/human/isg_2k_upstream.fa"

# Load DataFrame (CSV expected with a column containing gene symbols)
df = pd.read_csv(df_path, sep=',')

# If the DataFrame doesn't already contain Ensembl IDs, try to find a gene-symbol
# column and map symbols -> Ensembl IDs. Supported symbol column names:
# 'ensmbl_id', 'gene_id', 'gene', 'gene_name', 'symbol', 'gene_symbol', 'Gene'
if 'ensmbl_id' not in df.columns:
    # detect possible gene symbol column
    candidate_cols = ['gene_id', 'gene', 'genes', 'gene_name', 'symbol', 'gene_symbol', 'Gene']
    found = None
    for col in candidate_cols:
        if col in df.columns:
            found = col
            break

    if found is None:
        raise ValueError(f"No gene symbol column found in {df_path}. Expected one of: {candidate_cols} or 'ensmbl_id'.")

    # normalize to 'gene_id' for the mapping function
    if found != 'gene_id':
        df = df.rename(columns={found: 'gene_id'})

    # Map gene symbols to Ensembl IDs
    print(f"[INFO] Mapping gene symbols from column '{found}' to Ensembl IDs for species '{SPECIES}'...")
    df_mapped = batch_gene_to_ensembl(df, species=SPECIES)

    # Drop entries that couldn't be mapped
    unmapped = df_mapped['ensmbl_id'].isnull().sum()
    if unmapped > 0:
        print(f"[WARN] {unmapped} gene(s) could not be mapped to Ensembl IDs and will be skipped.")

    df = df_mapped.dropna(subset=['ensmbl_id']).reset_index(drop=True)

else:
    # already contains Ensembl IDs
    print("[INFO] Input already contains 'ensmbl_id' column; skipping mapping step.")


# Now fetch upstream sequences (uses Ensembl IDs in 'ensmbl_id')
get_sequences_to_fasta(df, out_fasta)
