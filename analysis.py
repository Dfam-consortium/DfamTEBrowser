import os
import subprocess
import shlex
from postprocess import (
    cluster_self_alignments,
    apply_depth_limit,
    chain_alignments,
    assign_colors_by_subject,
    assign_colors_by_pair,
)


def run_rmblastn(query, subject=None, db=None, rmblastn_path="rmblastn", params=None, matrix_dir=None, threads=1):
    """
    Run rmblastn for self/repeat alignments.
    Sets BLASTMAT environment variable.
    Returns: list of annotation dicts.
    """
    if not (subject or db):
        raise ValueError("Either subject or db must be provided.")
    if params is None:
        # self alignment uses subject=None
        #mask_level = 101 if subject else 80
        mask_level = 101
        #word_size = 7 if subject else 14
        word_size = 7
        params = (
          f"-num_alignments 9999999 -gapopen 20 -gapextend 5 "
          f"-mask_level {mask_level} -complexity_adjust -word_size {word_size} "
          "-xdrop_ungap 400 -xdrop_gap_final 200 -xdrop_gap 100 "
          "-min_raw_gapped_score 200 -dust no "
          "-outfmt \"6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen sstrand "
          "sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites qseq sseq\" "
          "-matrix comparison.matrix"
        )
        if not subject:
            params = params + f" -num_threads {threads}"

    cmd = [rmblastn_path] + shlex.split(params)
    if subject:
        cmd += ["-subject", subject]
    if db:
        cmd += ["-db", db]
    cmd += ["-query", query]

    # Set BLASTMAT environment variable
    env = os.environ.copy()
    if matrix_dir:
        env["BLASTMAT"] = matrix_dir
    else:
        raise ValueError("matrix_dir must be provided to set BLASTMAT")

    output = subprocess.check_output(cmd, text=True, env=env)
    return output


def run_blastx(query, db, blastx_path="blastx", threads=1):
    """
    Run blastx for protein alignments.
    Returns: list of annotation dicts.
    """
    outfmt = (
        "6 evalue perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen "
        "sstrand sseqid sstart send slen qframe sseq qseq"
    )
    cmd = [
        blastx_path,
        "-word_size", "2",
        "-evalue", "0.001",
        "-db", db,
        "-num_threads", str(threads),
        "-query", query,
        "-outfmt", outfmt,
    ]
    output = subprocess.check_output(cmd, text=True)
    return output

def run_ultra(query, ultra_path="ultra", threads=1):
    # TODO: Add thread support
    cmd = [ultra_path, query]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return proc.stdout.read()


def analyze_self(query, rmblastn_path, matrix_dir, threads=1):
    output = run_rmblastn(query=query, subject=query, rmblastn_path=rmblastn_path, matrix_dir=matrix_dir, threads=threads)
    annots = parse_blastn_output(output, "self")
    annots = cluster_self_alignments(annots)
    annots = apply_depth_limit(annots, "self")
    annots = assign_colors_by_pair(annots)
    return {
        "name": "Self Alignments",
        "type": "selfpair",
        "features": annots,
    }

def analyze_repeat(query, db, rmblastn_path, matrix_dir, threads=1):
    output = run_rmblastn(query=query, db=db, rmblastn_path=rmblastn_path, matrix_dir=matrix_dir, threads=threads)
    annots = parse_blastn_output(output, "repeat")
    annots = apply_depth_limit(annots, "repeat", max_depth=10)
    chains = chain_alignments(annots)
    # Assign color to chains directly (chain-level coloring)
    from colors import DISTINCT_COLORS
    color_by_subject = {}
    color_idx = 0
    for chain in chains:
        subject_id = chain['name']
        if subject_id not in color_by_subject:
            color_by_subject[subject_id] = DISTINCT_COLORS[color_idx % len(DISTINCT_COLORS)]
            color_idx += 1
        chain['color'] = color_by_subject[subject_id]
    return {
        "name": "TE DNA Homology",
        "type": "chain",
        "features": chains,
    }

def analyze_protein(query, db, blastx_path, threads=1):
    output = run_blastx(query=query, db=db, blastx_path=blastx_path, threads=threads)
    annots = parse_blastx_output(output)
    annots = apply_depth_limit(annots, "protein")
    annots = assign_colors_by_subject(annots)
    return {
        "name": "TE Protein Homology",
        "type": "annotation",
        "features": annots,
    }

def analyze_ultra(query, ultra_path, threads=1):
    output = run_ultra(query=query, ultra_path=ultra_path, threads=threads)
    annots = parse_ultra_output(output)
    for ann in annots:
        ann['color'] = "blue"
    return {
        "name": "ULTRA",
        "type": "annotation",
        "features": annots,
    }


def parse_ultra_output(output):
    """
    Parses the output of ultra into a list of annotation dicts.

    Each line is expected to have at least 4 tab-separated columns:
    [0]=name, [1]=start, [2]=end, [3]=name (or annotation label)
    Returns a list of dicts suitable for visualization tracks.
    """
    annots = []
    for line in output.strip().split('\n'):
        line = line.strip()
        if not line or not (line[0].isalnum() or line[0] == '_'):
            continue
        data = line.split('\t')
        if len(data) < 4:
            continue
        # Usually: [0]=name, [1]=start, [2]=end, [3]=name
        try:
            start = int(data[1])
            end = int(data[2])
        except Exception:
            continue
        name = data[3]
        # Shorten name if too long
        if len(name) > 14:
            name = f"{name[:4]}..{name[-4:]}[{len(name)}]"
        annots.append({
            'ref_start': start,
            'ref_end': end,
            'name': name,
            'color': 'blue',  # default, can be overridden by caller
            'type': 'ultra',
            'orient': '+'
        })
    return annots

def parse_blastn_output(output, type_):
    annotations = []
    for line in output.strip().split("\n"):
        d = line.split("\t")
        if len(d) < 20:
            continue
        score = float(d[0])
        ref_start = int(d[5])
        ref_end = int(d[6])
        cons_start = int(d[10])
        cons_end = int(d[11])
        cons_len = int(d[12])
        # NOTE: NCBI blast switches to descending order when orientation is "minus".
        #       We should correct that here.
        if d[8] and d[8] == 'minus':
            orient = '-'
            cons_start, cons_end = cons_end, cons_start
        else:
            orient = '+'
        name = d[9] or 'self'
        ref_seq = d[18]
        cons_seq = d[19]
        oseq = cons_seq
        # For self, skip perfect diagonal
        if type_ == 'self' and ref_start == cons_start and ref_end == cons_end:
            continue
        annotations.append({
            'ref_start': ref_start,
            'ref_end': ref_end,
            'cons_start': cons_start,
            'cons_end': cons_end,
            'cons_len': cons_len,
            'orient': orient,
            'name': name,
            'score': score,
            'ref_seq': ref_seq,
            'cons_seq': cons_seq,
            'oseq': oseq,
            'type': type_,
        })
    return annotations


def parse_blastx_output(blastx_lines, type_="protein"):
    annots = []
    for line in blastx_lines.strip().split("\n"):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 16:
            continue

        # Assign fields according to your outfmt
        evalue = fields[0]
        qseqid = fields[4]
        qstart = int(fields[5])
        qend = int(fields[6])
        qlen = int(fields[7])
        sseqid = fields[9]
        sstart = int(fields[10])
        send = int(fields[11])
        slen = int(fields[12])
        oseq = fields[14] # sseq
        seq = fields[15]  # qseq

        # Standardize coordinates: always ref_start <= ref_end
        if qstart <= qend:
            ref_start, ref_end, orient = qstart, qend, "+"
        else:
            ref_start, ref_end, orient = qend, qstart, "-"

        if sstart <= send:
            cons_start, cons_end = sstart, send
        else:
            cons_start, cons_end = send, sstart

        cons_len = slen

        try:
            score = float(evalue)
        except Exception:
            score = 0

        # The visualization doesn't handle insertions (relative to the reference) 
        # as of yet, so we need to filter them out.
        filtered_oseq = ''.join(o for s, o in zip(seq, oseq) if s != '-')

        annots.append({
            'ref_start': ref_start,
            'ref_end': ref_end,
            'cons_start': cons_start,
            'cons_end': cons_end,
            'cons_len': cons_len,
            'orient': orient,
            'name': sseqid,
            'score': score,
            'ref_seq': qseqid,
            'cons_seq': sseqid,
            'oseq': filtered_oseq,
            'type': type_,
        })
    return annots
