#!/usr/bin/env python3
"""
Create a background FASTA of random, non-ISG Homo sapiens promoters (2 kb upstream)

• Pulls all Ensembl gene IDs & symbols via /overlap/region per chromosome
• Removes user-supplied ISGs (symbols and/or Ensembl IDs)
• Samples N genes (default N = number of records in your positive FASTA)
• Fetches 2 kb upstream sequences in parallel (strand-aware)
• Saves background FASTA for HOMER (-bg in FASTA mode)

Author: Marek (2025-11-08)
"""

from __future__ import annotations
import argparse
import concurrent.futures
import random
import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ------------------------- constants / config -------------------------------

ENSEMBL = "https://rest.ensembl.org"
SPECIES = "homo_sapiens"
CHROMS  = [*map(str, range(1, 23)), "X", "Y", "MT"]  # primary assembly
UA      = {"User-Agent": "nonISG-bg/1.0 (contact: Marek)"}
TIMEOUT = 30
RETRIES = 4
BACKOFF = 1.5

UPSTREAM_BP_DEFAULT = 2000
THREADS_DEFAULT     = 8
RANDOM_SEED         = 42

ENSG_RE = re.compile(r"^ENSG\d+", re.I)
ENST_RE = re.compile(r"^ENST\d+", re.I)
NM_RE   = re.compile(r"^NM_\d+", re.I)

# ---------------------------- HTTP helpers ----------------------------------

def _get_json(url: str, params: Optional[dict] = None) -> dict:
    """GET JSON with simple retry/backoff."""
    for attempt in range(1, RETRIES + 1):
        try:
            # Ensembl expects Accept header for content negotiation, not content-type in params
            headers = UA.copy()
            headers["Accept"] = "application/json"
            params_local = dict(params) if params else {}
            # Remove content-type from params if present
            params_local.pop("content-type", None)
            params_local.pop("content_type", None)

            r = requests.get(url, params=params_local, headers=headers, timeout=TIMEOUT)
            if r.status_code == 429:
                # rate limit – use Retry-After if present
                sleep_s = int(r.headers.get("Retry-After", 2))
                time.sleep(sleep_s)
                continue
            r.raise_for_status()
            return r.json()
        except Exception:
            if attempt == RETRIES:
                raise
            time.sleep(BACKOFF ** attempt)
    raise RuntimeError("unreachable")

def _get_text(url: str, params: Optional[dict] = None) -> str:
    """GET text (FASTA) with retry/backoff."""
    for attempt in range(1, RETRIES + 1):
        try:
            # Ensembl expects Accept header for text/x-fasta
            headers = UA.copy()
            headers["Accept"] = "text/x-fasta"
            params_local = dict(params) if params else {}
            # Remove content-type from params if present
            params_local.pop("content-type", None)
            params_local.pop("content_type", None)

            r = requests.get(url, params=params_local, headers=headers, timeout=TIMEOUT)
            if r.status_code == 429:
                sleep_s = int(r.headers.get("Retry-After", 2))
                time.sleep(sleep_s)
                continue
            r.raise_for_status()
            return r.text
        except Exception:
            if attempt == RETRIES:
                raise
            time.sleep(BACKOFF ** attempt)
    raise RuntimeError("unreachable")

# ----------------------------- core logic -----------------------------------

def fetch_gene_index(chroms: Optional[List[str]] = None) -> Dict[str, Tuple[str, str, int, int, int]]:
    """
    Return dict: ensembl_gene_id -> (symbol, chrom, start, end, strand)
    Uses /lookup/genome to get all genes in one call.
    """
    if chroms is None:
        chroms = CHROMS

    out: Dict[str, Tuple[str, str, int, int, int]] = {}

    # Use the lookup endpoint which doesn't have region size limits
    # We'll query each chromosome separately but without coordinates
    for chrom in chroms:
        print(f"[INFO] Fetching genes for chromosome {chrom}...")
        # Use the xrefs/symbol endpoint is too slow, so we use overlap with feature=gene
        # but we need to handle the 5Mb limit by chunking
        # Alternative: use the lookup/genome endpoint or FTP dump

        # Better approach: get genes via REST API lookup endpoint
        # Since overlap has 5Mb limit, we'll chunk each chromosome
        chunk_size = 5_000_000  # 5 Mb chunks
        chrom_length = 250_000_000  # Max human chromosome length (chr1 ~248Mb)

        for start in range(1, chrom_length, chunk_size):
            end = start + chunk_size - 1
            url = f"{ENSEMBL}/overlap/region/{SPECIES}/{chrom}:{start}-{end}"
            try:
                data = _get_json(url, params={"feature": "gene"})
                if not data:  # No more genes in this chromosome
                    break
                for g in data:
                    gid   = g.get("id")
                    sym   = (g.get("external_name") or "").strip()
                    g_start = int(g["start"])
                    g_end   = int(g["end"])
                    strand= int(g["strand"])
                    if gid and gid not in out:  # Avoid duplicates at chunk boundaries
                        out[gid] = (sym, chrom, g_start, g_end, strand)
            except Exception as e:
                # If we get an error about the region, we've probably gone past the chromosome end
                if "not found" in str(e).lower() or "invalid" in str(e).lower():
                    break
                # For other errors, we might want to continue
                print(f"[WARN] Error fetching {chrom}:{start}-{end}: {e}")
                continue

    print(f"[INFO] Fetched {len(out)} total genes")
    return out

def load_isg_terms(path: Path) -> Tuple[Set[str], Set[str]]:
    """
    Load ISGs from a text or CSV file (first column used).
    Returns (symbols_set, ensembl_ids_set). Auto-detects ENSG/ENST/NM_*.
    """
    symbols: Set[str] = set()
    ids: Set[str] = set()
    # read first column from CSV/TSV/TXT
    with path.open() as f:
        for line in f:
            if not line.strip():
                continue
            # split on comma or tab; take first token
            tok = line.strip().split(",")[0].split("\t")[0].strip()
            if not tok:
                continue
            if ENSG_RE.match(tok) or ENST_RE.match(tok):
                ids.add(tok.split(".")[0])
            elif NM_RE.match(tok):
                # RefSeq transcript – we do not resolve NM_ here; user may supply symbols too.
                symbols.add(tok)  # keep to allow exact exclusion if overlapping annotation carries NM_
            else:
                # assume HGNC symbol
                symbols.add(tok.upper())
    return symbols, ids

def subtract_isgs(index: Dict[str, Tuple[str,str,int,int,int]],
                  isg_symbols: Set[str], isg_ids: Set[str]) -> List[str]:
    """
    Return list of Ensembl gene IDs after removing any that match by symbol or ID.
    """
    keep: List[str] = []
    isg_ids_norm = {x.upper() for x in isg_ids}
    isg_syms_norm = {x.upper() for x in isg_symbols}
    for gid, (sym, *_rest) in index.items():
        if gid.upper() in isg_ids_norm:
            continue
        if sym and sym.upper() in isg_syms_norm:
            continue
        keep.append(gid)
    return keep

def sample_ids(ids: List[str], n: int, seed: int = RANDOM_SEED) -> List[str]:
    if n <= 0:
        return []
    if n > len(ids):
        raise ValueError(f"Requested sample size {n} > available pool {len(ids)}")
    rnd = random.Random(seed)
    return rnd.sample(ids, n)

def upstream_region(start: int, end: int, strand: int, upstream: int) -> Tuple[int,int,int]:
    """
    Return (region_start, region_end, strand) 1-based closed interval for upstream flank.
    For '+' strand: [start-upstream, start-1]; for '-' strand: [end+1, end+upstream], strand stays -1.
    Clamp to 1 if needed (left boundary).
    """
    if strand >= 0:
        s = max(1, start - upstream)
        e = max(1, start - 1)
        return s, e, 1
    else:
        s = end + 1
        e = end + upstream
        return s, e, -1

def fetch_upstream_fasta(gid: str, meta: Tuple[str,str,int,int,int], upstream_bp: int) -> Optional[SeqRecord]:
    sym, chrom, start, end, strand = meta
    rs, re, st = upstream_region(start, end, strand, upstream_bp)
    # Ensembl wants strand as 1 or -1 in the URL
    url = f"{ENSEMBL}/sequence/region/{SPECIES}/{chrom}:{rs}..{re}:{st}"
    txt = _get_text(url, params={"content-type":"text/x-fasta"})
    lines = txt.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        return None
    seq = "".join(lines[1:]).upper()
    if len(seq) == 0:
        return None
    # Some regions near chromosome starts can be shorter than upstream_bp; we still keep them.
    rec_id = gid
    desc = f"{gid}|SYMBOL={sym or 'NA'}|CHR={chrom}|STRAND={st}|UPSTREAM={upstream_bp}"
    return SeqRecord(Seq(seq), id=rec_id, description=desc)

def write_bg_fasta(sampled_ids: List[str],
                   index: Dict[str, Tuple[str,str,int,int,int]],
                   upstream_bp: int,
                   threads: int,
                   out_fa: Path) -> int:
    recs: List[SeqRecord] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(fetch_upstream_fasta, gid, index[gid], upstream_bp): gid for gid in sampled_ids}
        for fut in concurrent.futures.as_completed(futs):
            try:
                rec = fut.result()
                if rec:
                    recs.append(rec)
            except Exception:
                # skip failed; Ensembl occasionally has gaps
                continue
    if not recs:
        raise RuntimeError("No background sequences were fetched.")
    SeqIO.write(recs, out_fa, "fasta")
    return len(recs)

def count_fasta_records(path: Path) -> int:
    try:
        return sum(1 for _ in SeqIO.parse(path, "fasta"))
    except Exception:
        return 0

# ------------------------------ CLI -----------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Build human non-ISG upstream (2 kb) background FASTA for HOMER.")
    ap.add_argument("-i", "--isg-list", required=True, help="ISG list file (symbols and/or Ensembl IDs). First column used.")
    ap.add_argument("-p", "--positives-fasta", required=False, help="Positive FASTA; if -n is omitted, N = number of records here.")
    ap.add_argument("-n", "--num", type=int, default=None, help="How many background genes to sample. Default: #records in positives FASTA.")
    ap.add_argument("-u", "--upstream", type=int, default=UPSTREAM_BP_DEFAULT, help="Upstream length (bp). Default 2000.")
    ap.add_argument("-t", "--threads", type=int, default=THREADS_DEFAULT, help="Threads for fetching sequences. Default 8.")
    ap.add_argument("-c", "--chromosomes", default=None, help="Comma-separated list of chromosomes to query (e.g., '1,2,3,X'). Default: all chromosomes.")
    ap.add_argument("-o", "--out-fasta", required=True, help="Output FASTA path for background (e.g., nonISG_up2kb.bg.fa)")
    args = ap.parse_args()

    isg_file = Path(args.isg_list)
    out_fa   = Path(args.out_fasta)
    upstream = int(args.upstream)
    threads  = int(args.threads)

    # Parse chromosomes
    chroms = CHROMS
    if args.chromosomes:
        chroms = [c.strip() for c in args.chromosomes.split(',')]
        print(f"[INFO] Querying chromosomes: {', '.join(chroms)}")

    # Determine N
    N: Optional[int] = args.num
    if N is None:
        if not args.positives_fasta:
            ap.error("Either provide -n/--num or --positives-fasta to infer N.")
        posN = count_fasta_records(Path(args.positives_fasta))
        if posN <= 0:
            ap.error("Could not determine N from positives FASTA (no records). Provide -n explicitly.")
        N = posN

    # Load ISGs
    isg_symbols, isg_ids = load_isg_terms(isg_file)

    # Index all human genes
    index = fetch_gene_index(chroms)

    # Subtract ISGs by symbol or Ensembl ID
    pool_ids = subtract_isgs(index, isg_symbols, isg_ids)
    if len(pool_ids) < N:
        raise SystemExit(f"Background pool too small after ISG removal: pool={len(pool_ids)} < N={N}")

    # Sample and fetch sequences
    sampled = sample_ids(pool_ids, N)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    got = write_bg_fasta(sampled, index, upstream, threads, out_fa)

    print(f"[OK] Wrote {got} background sequences to: {out_fa}")

if __name__ == "__main__":
    main()
