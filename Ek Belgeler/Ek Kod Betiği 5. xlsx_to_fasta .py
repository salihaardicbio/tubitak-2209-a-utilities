#!/usr/bin/env python3
"""
Collects 'target' (accession) and 'tSeq' (sequence) columns from every .xlsx
file in a folder and writes them all into a single FASTA file.

Only sequences with E-value < threshold are included.

Usage:
    python xlsx_to_fasta.py [folder] [-o output.fasta] [-e 0.01]
"""

import argparse
import glob
import os
import sys

import openpyxl


def find_columns(header_row):
    """Find target, tSeq and E-value columns (case-insensitive)."""

    target_idx = None
    tseq_idx = None
    evalue_idx = None

    for i, cell in enumerate(header_row):
        if cell is None:
            continue

        name = str(cell).strip().lower()

        if name == "target":
            target_idx = i

        elif name == "tseq":
            tseq_idx = i

        elif (
            "eval" in name
            or "e-value" in name
            or "evalue" in name
            or "e_value" in name
            or "expect" in name
        ):
            evalue_idx = i

    return target_idx, tseq_idx, evalue_idx


def wrap_seq(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "folder",
        nargs="?",
        default=".",
        help="Folder containing .xlsx files (default=current directory)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="combined.fasta",
        help="Output FASTA file",
    )

    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=0.01,
        help="Maximum allowed E-value (default=0.01)",
    )

    args = parser.parse_args()

    xlsx_files = sorted(glob.glob(os.path.join(args.folder, "*.xlsx")))
    xlsx_files = [
        f for f in xlsx_files
        if not os.path.basename(f).startswith("~$")
    ]

    if not xlsx_files:
        print(f"No .xlsx files found in {os.path.abspath(args.folder)}")
        sys.exit(1)

    print(f"Found {len(xlsx_files)} xlsx file(s):")
    for f in xlsx_files:
        print("  ", os.path.basename(f))

    seen_ids = {}
    seen_pairs = set()
    records = []

    for path in xlsx_files:

        try:
            wb = openpyxl.load_workbook(
                path,
                read_only=True,
                data_only=True,
            )
        except Exception as e:
            print(f"Could not open {path}: {e}")
            continue

        for ws in wb.worksheets:

            rows = ws.iter_rows(values_only=True)

            try:
                header = next(rows)
            except StopIteration:
                continue

            target_idx, tseq_idx, evalue_idx = find_columns(header)

            if (
                target_idx is None
                or tseq_idx is None
                or evalue_idx is None
            ):
                continue

            kept = 0

            for row in rows:

                if row is None:
                    continue

                if max(target_idx, tseq_idx, evalue_idx) >= len(row):
                    continue

                accession = row[target_idx]
                sequence = row[tseq_idx]
                evalue = row[evalue_idx]

                if accession is None or sequence is None or evalue is None:
                    continue

                try:
                    evalue = float(evalue)
                except Exception:
                    continue

                if evalue >= args.evalue:
                    continue

                accession = str(accession).strip()
                sequence = str(sequence).strip().upper()

                if not accession or not sequence:
                    continue

                pair = (accession, sequence)

                if pair in seen_pairs:
                    continue

                seen_pairs.add(pair)

                if accession in seen_ids:
                    seen_ids[accession] += 1
                    header_id = f"{accession}_{seen_ids[accession]}"
                else:
                    seen_ids[accession] = 0
                    header_id = accession

                records.append((header_id, sequence))
                kept += 1

            if kept:
                print(f"{os.path.basename(path)} [{ws.title}] -> {kept} sequences kept")

    if not records:
        print("No sequences passed the E-value filter.")
        sys.exit(1)

    with open(args.output, "w") as out:
        for header, sequence in records:
            out.write(f">{header}\n")
            out.write(wrap_seq(sequence))
            out.write("\n")

    print(f"\nWrote {len(records)} sequences to:")
    print(os.path.abspath(args.output))


if __name__ == "__main__":
    main()