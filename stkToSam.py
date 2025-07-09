#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    stkToSam.py

    Usage:
        ./stkToSam.py [--help]
                      --input=aln.stk
                      [--sam=out.sam]
                      [--ref=ref.fa]

    Convert a Stockholm multiple-sequence alignment (MSA) to a SAM file and
    ungapped reference FASTA file. This tool requires that the Stockholm
    file provides the reference annotation (#=GC RF line).

    The resulting SAM file may be compressed and indexed with samtools.

    Args:
        -h, --help       : Show this help message and exit
        --input          : Input Stockholm file
        --sam            : Output SAM file (default: out.sam)
        --ref            : Output reference FASTA file (default: ref.fa)

    SEE ALSO:
        samtools: https://www.htslib.org
        Dfam:     http://www.dfam.org

    AUTHOR(S):
        Robert Hubley <rhubley@systemsbiology.org>

    LICENSE:
        This code may be used in accordance with the Creative Commons
        Zero ("CC0") public domain dedication:
        https://creativecommons.org/publicdomain/zero/1.0/

    DISCLAIMER:
        This software is provided “AS IS” and any express or implied
        warranties, including, but not limited to, the implied warranties of
        merchantability and fitness for a particular purpose, are disclaimed.
        In no event shall the authors or the Dfam consortium members be
        liable for any direct, indirect, incidental, special, exemplary, or
        consequential damages (including, but not limited to, procurement of
        substitute goods or services; loss of use, data, or profits; or
        business interruption) however caused and on any theory of liability,
        whether in contract, strict liability, or tort (including negligence
        or otherwise) arising in any way out of the use of this software, even
        if advised of the possibility of such damage.
"""

import sys
import os
import argparse
from pathlib import Path
import textwrap
import re

coord_re = re.compile(r"^(\d+)-(\d+)(?:_([+-]))?$")

def _usage():
    """Print out docstring as program usage"""
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)

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
    ref_coord = 0
    started = False
    pos = None

    for r, q in zip(ref_cols, seq_cols):
        ref_has = r not in ".-"
        seq_has = q not in ".-"

        if ref_has:
            ref_coord += 1

        if not started:
            if not seq_has:
                continue
            started = True
            pos = ref_coord if ref_has else ref_coord + 1

        if ref_has and seq_has:
            ops.append("M")
        elif ref_has and not seq_has:
            ops.append("D")
        elif not ref_has and seq_has:
            ops.append("I")

    if not ops:
        return None, None

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
        return sid, '+'

    start, end, orient = match.groups()
    start = int(start)
    end = int(end)

    if orient in ('+', '-'):
        return base, orient
    else:
        return base, '+' if end >= start else '-'

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

def stk_to_sam(stk_path, sam_path="out.sam", ref_fa_path="ref.fa",
               ref_name="RF"):
    """
    Convert Stockholm MSA to SAM and reference FASTA.
    """
    rf_align, seqs, identifier = read_stockholm(stk_path)
    aln_len = len(rf_align)

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
                continue

            cigar_str = "".join(f"{l}{op}" for op, l in cigar)
            seq_no_gaps = aln.translate({ord('.'): None, ord('-'): None})
            base_id, orient = parse_orientation(sid)
            flag = 16 if orient == '-' else 0

            sam.write(
                f"{base_id}\t{flag}\t{identifier}\t{pos}\t0\t{cigar_str}\t*\t0\t0\t"
                f"{seq_no_gaps}\t*\n"
            )

def main(*args):
    class _CustomUsageAction(argparse.Action):
        def __init__(self, option_strings, dest, default=False, required=False, help=None):
            super(_CustomUsageAction, self).__init__(
                option_strings=option_strings, dest=dest,
                nargs=0, const=True, default=default,
                required=required, help=help)
        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--help", action=_CustomUsageAction)
    parser.add_argument("--input", type=Path, required=True, help="Input Stockholm file")
    parser.add_argument("--sam", type=Path, default="out.sam", help="Output SAM file")
    parser.add_argument("--ref", type=Path, default="ref.fa", help="Output reference FASTA")
    args = parser.parse_args()

    stk_to_sam(args.input, args.sam, args.ref)
    print("Wrote", args.ref, "and", args.sam)
    print("Compress and index with:")
    print(f"  samtools view -C --write-index -T {args.ref} -o out.cram {args.sam}")
    print("  samtools index out.cram")

if __name__ == '__main__':
    main(*sys.argv)

