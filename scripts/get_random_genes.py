#!/usr/bin/env python3
"""
Create a background FASTA of random, non-ISG Gallus gallus promoters
====================================================================

• Pulls *all* Ensembl gene IDs via REST (one call per chromosome)
• Removes user-supplied ISGs
• Randomly samples N genes
• Fetches 2 kb upstream sequences in parallel
• Saves background FASTA for HOMER (-fastaBg)

Author : Marek (2025-06-23)
"""
from __future__ import annotations
import random, concurrent.futures, time
from pathlib import Path
from typing import List
import requests
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from tqdm import tqdm

# ---------- configurable constants ------------------------------------------
ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES      = "gallus_gallus"
UPSTREAM_BP  = 2_000          # promoter length
CHROMS       = [*map(str, range(1, 29)), "Z", "W", "MT"]
THREADS      = 8
SAMPLE_SIZE  = 500
PARENT_DIR   = Path("/home/pospim/Desktop/bca")
ISG_FILE     = Path(PARENT_DIR, "papers/isg_ensembl.csv")  # ISG Ensembl IDs in CSV format
OUT_FASTA    = Path(f"{PARENT_DIR}/genomes/nonISG_upstream_2kb.fa")
filtered_genes_path = Path(f"{PARENT_DIR}/genomes/filtered_gene_ids.txt")
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------#
# 1. download every gene ID in gallus_gallus (≈ 19 k)                        #
# ---------------------------------------------------------------------------#
def fetch_gene_ids(chrom: str) -> List[str]:
    """
    Fetch all Ensembl gene IDs for a given chromosome.
    """
    url = (f"{ENSEMBL_REST}/overlap/region/{SPECIES}/{chrom}"
            f"?feature=gene;content-type=application/json")
    response = requests.get(url, timeout=30)
    response.raise_for_status()  # Raise an error for bad responses
    if not response.ok:
        raise Exception(f"[ERROR] Failed to fetch gene IDs for {chrom}: {response.status_code} {response.reason}")

    data = response.json()
    return [gene['id'] for gene in data]

def get_all_gene_ids() -> List[str]:
    """
    Fetch all Ensembl gene IDs for the specified species across all chromosomes.
    """
    all_gene_ids = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = {executor.submit(fetch_gene_ids, chrom): chrom for chrom in CHROMS}
        for future in tqdm(concurrent.futures.as_completed(futures),
                            total=len(CHROMS), desc="Fetching gene IDs"):
            chrom = futures[future]
            try:
                gene_ids = future.result()
                all_gene_ids.extend(gene_ids)
            except Exception as e:
                print(f"[ERROR] Failed to fetch gene IDs for {chrom}: {e}")
    return sorted(set(all_gene_ids))

all_ids = get_all_gene_ids()
print(all_ids)
print(len(all_ids))

# ---------------------------------------------------------------------------#
# 2. remove ISGs from the list of all gene IDs                               #
# ---------------------------------------------------------------------------#
def load_isg_ids(isg_file: Path) -> set[str]:
    """
    Load ISG Ensembl IDs from a CSV file.
    """
    df = pd.read_csv(isg_file, usecols=["ensembl_id"])
    return set(df["ensembl_id"].dropna().astype(str))

def filter_isg_ids(all_gene_ids: List[str], isg_ids: set[str]) -> List[str]:
    """
    Filter out ISG Ensembl IDs from the list of all gene IDs.
    """
    all_gene_set = set(all_gene_ids)
    filtered_gene_ids = all_gene_set - isg_ids
    return sorted(filtered_gene_ids)

isgs = load_isg_ids(ISG_FILE)
print(isgs)

filtered = filter_isg_ids(all_ids, isgs)
print(filtered)


with open(filtered_genes_path, "w") as f:
    for gene_id in filtered:
        f.write(f"{gene_id}\n")

print("Filtered gene IDs saved to 'filtered_gene_ids.txt'")
print(f"FILTERED: {len(filtered)}\nALL: {len(all_ids)}")

# ---------------------------------------------------------------------------#
# 3. sample N non-ISG genes                                                  #
# ---------------------------------------------------------------------------#
def sample_background_genes(filtered_gene_ids: List[str], sample_size: int) -> List[str]:
    """
    Randomly sample a specified number of non-ISG gene IDs.
    """
    if sample_size > len(filtered_gene_ids):
        raise ValueError(f"Sample size {sample_size} exceeds available non-ISG gene IDs {len(filtered_gene_ids)}.")
    random.seed(42)  # For reproducibility
    sampled_genes = random.sample(filtered_gene_ids, sample_size)
    return sampled_genes

sampled = sample_background_genes(filtered, SAMPLE_SIZE)
print(sampled)

# ---------------------------------------------------------------------------#
# 4. fetch upstream sequences in parallel                                    #
# ---------------------------------------------------------------------------#
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

def fetch_upstream_seq(ens_id: str,
                        species: str = SPECIES,
                          upstream: int = UPSTREAM_BP) -> SeqRecord | None:
    """
    Fetch the upstream sequence for a given Ensembl gene ID.
    """
    try:
        coords = get_gene_coordinates(ens_id)
        start, end, strand = coords["gene_start"], coords["gene_end"], coords["strand"]
        chrom = coords["chromosome"]

        if strand == 1:
            region_start = max(1, start - upstream)
            region_end = start - 1
        else:
            region_start = end + 1
            region_end = end + upstream
        seq_url = (f"{ENSEMBL_REST}/sequence/region/{species}/{chrom}:"
                   f"{region_start}..{region_end}:{strand}"
                   f"?content-type=text/x-fasta")
        r = requests.get(seq_url, timeout=30,
                         headers={"User-Agent": "ISG-bg/1.0 (Marek)"})
        time.sleep(1)
        r.raise_for_status()  # Raise an error for bad responses

        seq = "".join(r.text.splitlines()[1:]).upper()
        if not seq:
            raise ValueError(f"Empty sequence for {ens_id} at {chrom}:{region_start}-{region_end}")

        return SeqRecord(
            Seq(seq),
            id=ens_id,
            description=f"{ens_id} upstream {upstream} bp")
    except Exception as e:                       # noqa: BLE001
        tqdm.write(f"[warn] {ens_id}: {e}")
        return None

def fetch_sequences_to_fasta(gene_ids: List[str], output_fasta: Path, species: str = SPECIES, upstream: int = UPSTREAM_BP):
    """
    Fetch upstream sequences for a list of gene IDs and save to a FASTA file.
    """
    records = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS) as executor:
        futures = {executor.submit(fetch_upstream_seq, ens_id, species, upstream): ens_id for ens_id in gene_ids}
        for future in tqdm(concurrent.futures.as_completed(futures),
                            total=len(gene_ids), desc="Fetching sequences"):
            ens_id = futures[future]
            try:
                record = future.result()
                if record:
                    records.append(record)
            except Exception as e:
                tqdm.write(f"[ERROR] Failed to fetch sequence for {ens_id}: {e}")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"[INFO] Saved {len(records)} sequences to {output_fasta}")

fetch_sequences_to_fasta(filtered, OUT_FASTA)

# ---------------------------------------------------------------------------#
# Main execution flow
# ---------------------------------------------------------------------------#
def main():
    t0 = time.perf_counter()

    # Step 1: Fetch all gene IDs
    print("[INFO] Fetching all gene IDs...")
    all_gene_ids = get_all_gene_ids()
    print(f"[INFO] Total gene IDs fetched: {len(all_gene_ids)}")

    # Step 2: Load ISG IDs and filter them out
    print("[INFO] Loading ISG IDs...")
    isg_ids = load_isg_ids(ISG_FILE)
    print(f"[INFO] Total ISG IDs loaded: {len(isg_ids)}")

    filtered_gene_ids = filter_isg_ids(all_gene_ids, isg_ids)
    print(f"[INFO] Total non-ISG gene IDs available: {len(filtered_gene_ids)}")

    # Step 3: Sample background genes
    print(f"[INFO] Sampling {SAMPLE_SIZE} non-ISG genes...")
    sampled_genes = sample_background_genes(filtered_gene_ids, SAMPLE_SIZE)
    print(f"[INFO] Sampled {len(sampled_genes)} non-ISG genes.")

    # Step 4: Fetch upstream sequences and save to FASTA
    print("[INFO] Fetching upstream sequences...")
    fetch_sequences_to_fasta(sampled_genes, OUT_FASTA)