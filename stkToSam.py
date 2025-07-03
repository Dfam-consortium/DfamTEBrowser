#!/usr/bin/env python3
"""
stkToSam.py
==================
Convert a Stockholm multiple-sequence alignment (MSA) to a SAM file and ungapped reference FASTA, *without any external libraries*.
"""

import sys
from pathlib import Path
import textwrap
import re

coord_re = re.compile(r"^(\d+)-(\d+)(?:_([+-]))?$")

# ───────────────────────────────────────────────────────── helpers ───
def collapse_ops(ops):
    """Turn ['M','M','I',…] into [('M',2), ('I',1)…]."""
    if not ops:
        return []
    out = [[ops[0], 1]]
    for op in ops[1:]:
        if op == out[-1][0]:
            out[-1][1] += 1
        else:
            out.append([op, 1])
    return out

def cigar_and_pos(ref_cols, seq_cols):
    """
    Build CIGAR and starting reference position for one sequence row.

    Returns (cigar_tuples, pos) or (None, None) if the sequence is all gaps.
    `pos` is 1‑based for SAM.
    """
    ops = []
    ref_coord = 0          # running count of reference bases seen (0-based)
    started = False
    pos = None

    for r, q in zip(ref_cols, seq_cols):
        ref_has = r not in ".-"
        seq_has = q not in ".-"

        if ref_has:
            ref_coord += 1

        if not started:
            if not seq_has:
                continue  # Still in leading gap region
            started = True
            pos = ref_coord if ref_has else ref_coord + 1

        # Emit CIGAR op starting from first aligned base
        if ref_has and seq_has:
            ops.append("M")  # Match/mismatch
        elif ref_has and not seq_has:
            ops.append("D")  # Deletion
        elif not ref_has and seq_has:
            ops.append("I")  # Insertion

    if not ops:
        return None, None  # Entire sequence is gaps

    # Collapse and remove trailing D ops
    collapsed = collapse_ops(ops)
    while collapsed and collapsed[-1][0] == "D":
        collapsed.pop()

    return collapsed, pos

def parse_orientation(sid):
    """
    Extract base ID and orientation from a Stockholm sequence ID.

    Returns:
        base_id  – ID without coordinates
        orient   – '+' or '-'
    """
    last_colon = sid.rfind(':')
    if last_colon == -1:
        return sid, '+'

    base = sid[:last_colon]
    coord_part = sid[last_colon + 1:]

    match = coord_re.match(coord_part)
    if not match:
        # Not a recognized coordinate/orient pattern; treat as pure ID
        return sid, '+'

    start, end, orient = match.groups()
    start = int(start)
    end = int(end)

    if orient in ('+', '-'):
        return base, orient
    else:
        # Infer from start and end
        return base, '+' if end >= start else '-'

# ────────────────────────────────────────── tiny Stockholm parser ───
def read_stockholm(path):
    """Return (rf_align_string, {id: aligned_string}, identifier)."""
    rf_parts = []
    seqs = {}
    identifier = "Unknown"

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("#=GC RF"):
                    rf_parts.append(line.split(None, 2)[2])
                if identifier == "Unknown" and line.startswith("#=GF ID"):
                    identifier = line.split(None, 2)[2]
                if line.startswith("#=GF AC"):
                    identifier = line.split(None, 2)[2]
                continue
            if line == "//":
                break
            sid, aln = line.split(None, 1)
            seqs.setdefault(sid, []).append(aln.strip())

    rf_align = "".join(rf_parts)
    seqs = {sid: "".join(parts) for sid, parts in seqs.items()}
    return rf_align, seqs, identifier

# ───────────────────────────────────────────── conversion driver ───
def stk_to_sam(stk_path, sam_path="out.sam", ref_fa_path="ref.fa",
               ref_name="RF"):
    """
    Convert Stockholm MSA to SAM and reference FASTA.
    """
    rf_align, seqs, identifier = read_stockholm(stk_path)
    aln_len = len(rf_align)

    # reference FASTA
    ref_seq = "".join(c for c in rf_align if c not in ".-")
    with open(ref_fa_path, "w") as fh:
        fh.write(f">{identifier}\n")
        fh.write("\n".join(textwrap.wrap(ref_seq, 60)))
        fh.write("\n")

    with open(sam_path, "w") as sam:
        sam.write("@HD\tVN:1.6\tSO:unknown\n")
        sam.write(f"@SQ\tSN:{identifier}\tLN:{len(ref_seq)}\n")
        sam.write("@PG\tID:stk_to_sam\tPN:stk_to_sam\n")

        for sid, aln in seqs.items():
            if len(aln) != aln_len:
                raise ValueError(f"ERROR: {sid} length {len(aln)} ≠ RF length {aln_len}")

            cigar, pos = cigar_and_pos(rf_align, aln)
            if cigar is None:
                continue                      # row is all gaps

            cigar_str = "".join(f"{l}{op}" for op, l in cigar)
            seq_no_gaps = aln.translate({ord('.'): None, ord('-'): None})

            base_id, orient = parse_orientation(sid)
            flag = 16 if orient == '-' else 0  # Reverse strand flag if needed
            sam.write(
                f"{base_id}\t{flag}\t{identifier}\t{pos}\t0\t{cigar_str}\t*\t0\t0\t"
                f"{seq_no_gaps}\t*\n"
            )

# ─────────────────────────────────────────────── entry point ───
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert Stockholm MSA to SAM and reference FASTA"
    )
    parser.add_argument("input", type=Path, help="Input Stockholm file")
    parser.add_argument("--sam", type=Path, default="out.sam", help="Output SAM file")
    parser.add_argument("--ref", type=Path, default="ref.fa", help="Output reference FASTA")
    args = parser.parse_args()

    stk_to_sam(args.input, args.sam, args.ref)
    print("Wrote", args.ref, "and", args.sam)
    print("Compress and index with:")
    print(f"  samtools view -C --write-index -T {args.ref} -o out.cram {args.sam}")
    print("  samtools index out.cram")

if __name__ == "__main__":
    main()
