"""
Fetch missing hit_sequence entries from UniProt using hit_accession,
then build a FASTA file (accession as sequence name) combining
existing + newly fetched sequences.

Usage:
    python fetch_missing_sequences.py input.csv output_prefix

Outputs:
    output_prefix_filled.csv   -> original CSV with hit_sequence filled in where possible
    output_prefix.fasta        -> FASTA file, one record per row, header = hit_accession
"""

import sys
import csv
import time
import requests
from collections import defaultdict

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
BATCH_SIZE = 100          # accessions per UniProt query
SLEEP_BETWEEN_BATCHES = 1  # seconds, be polite to the API

import csv
with open("PSblast_results.csv", newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    print("Columns:", reader.fieldnames)
    rows = list(reader)

for r in rows[:20]:
    val = r.get("hit_sequence")
    print(repr(r.get("hit_accession")), "->", repr(val))


    
def load_rows(csv_path):
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames
    return rows, fieldnames


def find_empty_seq_accessions(rows):
    """Return the accessions of rows where hit_sequence is empty/missing."""
    accs = []
    for r in rows:
        seq = (r.get("hit_sequence") or "").strip()
        if not seq:
            acc = (r.get("hit_accession") or "").strip()
            if acc:
                accs.append(acc)
    return accs


def chunked(lst, size):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]


def fetch_uniprot_sequences(accessions):
    """
    Query UniProt in batches, return dict {accession: sequence}.
    Accessions not found in UniProt will simply be absent from the result.
    """
    seq_by_acc = {}
    unique_accs = sorted(set(accessions))

    for batch in chunked(unique_accs, BATCH_SIZE):
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        params = {"query": query, "format": "fasta"}

        try:
            resp = requests.get(UNIPROT_STREAM_URL, params=params, timeout=60)
            resp.raise_for_status()
        except requests.RequestException as e:
            print(f"  ! Batch failed ({len(batch)} accessions): {e}")
            time.sleep(SLEEP_BETWEEN_BATCHES)
            continue

        # Parse simple FASTA: header format is like
        # >sp|Q6GMR7|SOME_NAME description
        current_acc = None
        current_seq_lines = []

        def flush():
            if current_acc and current_seq_lines:
                seq_by_acc[current_acc] = "".join(current_seq_lines)

        for line in resp.text.splitlines():
            if line.startswith(">"):
                flush()
                current_seq_lines = []
                parts = line[1:].split("|")
                # parts[1] is usually the accession for sp|ACC|NAME headers
                current_acc = parts[1] if len(parts) > 1 else parts[0].split()[0]
            else:
                current_seq_lines.append(line.strip())
        flush()

        found = len(seq_by_acc)
        print(f"  Batch of {len(batch)} accessions -> {found} sequences retrieved so far")
        time.sleep(SLEEP_BETWEEN_BATCHES)

    missing = [a for a in unique_accs if a not in seq_by_acc]
    if missing:
        print(f"\nWarning: {len(missing)} accession(s) not found on UniProt "
              f"(likely non-UniProt IDs, e.g. RefSeq/PDB):")
        for m in missing[:20]:
            print(f"    {m}")
        if len(missing) > 20:
            print(f"    ... and {len(missing) - 20} more")

    return seq_by_acc


def fill_in_sequences(rows, seq_by_acc):
    filled_count = 0
    for r in rows:
        seq = (r.get("hit_sequence") or "").strip()
        if not seq:
            acc = (r.get("hit_accession") or "").strip()
            if acc in seq_by_acc:
                r["hit_sequence"] = seq_by_acc[acc]
                filled_count += 1
    return filled_count


def write_csv(rows, fieldnames, out_path):
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(rows, out_path):
    """One FASTA record per row that has a sequence, accession as the name."""
    seen_names = defaultdict(int)
    with open(out_path, "w") as f:
        for r in rows:
            seq = (r.get("hit_sequence") or "").strip()
            acc = (r.get("hit_accession") or "").strip()
            if not seq or not acc:
                continue
            # disambiguate duplicate accessions (e.g. same hit appearing twice)
            seen_names[acc] += 1
            name = acc if seen_names[acc] == 1 else f"{acc}_{seen_names[acc]}"
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python fetch_missing_sequences.py input.csv output_prefix")
        sys.exit(1)

    csv_path, out_prefix = sys.argv[1], sys.argv[2]

    print(f"Loading {csv_path} ...")
    rows, fieldnames = load_rows(csv_path)
    print(f"  {len(rows)} rows loaded")

    empty_accs = find_empty_seq_accessions(rows)
    print(f"  {len(empty_accs)} rows have an empty hit_sequence")

    if empty_accs:
        print("Fetching missing sequences from UniProt ...")
        seq_by_acc = fetch_uniprot_sequences(empty_accs)
        filled = fill_in_sequences(rows, seq_by_acc)
        print(f"  Filled in {filled} / {len(empty_accs)} missing sequences")
    else:
        print("No empty hit_sequence rows found — nothing to fetch.")

    csv_out = f"{out_prefix}_filled.csv"
    fasta_out = f"{out_prefix}.fasta"

    write_csv(rows, fieldnames, csv_out)
    write_fasta(rows, fasta_out)

    print(f"\nDone.\n  CSV  -> {csv_out}\n  FASTA -> {fasta_out}")


if __name__ == "__main__":
    main()
