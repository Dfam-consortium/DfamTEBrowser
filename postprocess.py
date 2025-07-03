from colors import DISTINCT_COLORS
import collections

MAX_TRACK_DEPTH = 10

def cluster_self_alignments(alignments):
    aligns = sorted(alignments, key=lambda a: (a['ref_start'], -a['ref_end']))
    used_aligns = set()
    results = []
    for i, align in enumerate(aligns):
        if align['ref_start'] == align['cons_start'] and align['ref_end'] == align['cons_end']:
            continue
        key = (align['ref_start'], align['ref_end'], align['cons_start'], align['cons_end'])
        key_sym = (align['cons_start'], align['cons_end'], align['ref_start'], align['ref_end'])
        if key in used_aligns or key_sym in used_aligns:
            continue
        range_max = align['ref_end']
        cluster = [align]
        ends = [align['ref_end']]
        for j, other in enumerate(aligns):
            if j == i:
                continue
            other_key = (other['ref_start'], other['ref_end'], other['cons_start'], other['cons_end'])
            other_key_sym = (other['cons_start'], other['cons_end'], other['ref_start'], other['ref_end'])
            if other_key in used_aligns or other_key_sym in used_aligns:
                continue
            if other['ref_start'] == other['cons_start'] and other['ref_end'] == other['cons_end']:
                continue
            # Clustering condition (allowing small slop of 5 positions)
            if (other['ref_start'] > align['ref_start'] - 5 and
                other['ref_start'] < align['ref_start'] + 5 and
                other['cons_start'] < align['ref_end']):
                if other['ref_end'] > range_max:
                    range_max = other['ref_end']
                cluster.append(other)
                ends.append(other['ref_end'])
        if len(cluster) > 2:
            max_query = align['ref_end']
            for other in cluster:
                okey = (other['ref_start'], other['ref_end'], other['cons_start'], other['cons_end'])
                okey_sym = (other['cons_start'], other['cons_end'], other['ref_start'], other['ref_end'])
                used_aligns.add(okey)
                used_aligns.add(okey_sym)
                if other['ref_end'] > max_query:
                    max_query = other['ref_end']
            results.append({
                'ref_start': align['ref_start'],
                'ref_end': max_query,
                'cons_start': align['ref_start'],
                'cons_end': max_query,
                'orient': "+",
                'name': "tandem_cluster",
                'score': align['score'],
                'type': "tandem",
                'ends': ends
            })
        else:
            for other in cluster:
                okey = (other['ref_start'], other['ref_end'], other['cons_start'], other['cons_end'])
                okey_sym = (other['cons_start'], other['cons_end'], other['ref_start'], other['ref_end'])
                used_aligns.add(okey)
                used_aligns.add(okey_sym)
                results.append(other)
    return results

def apply_depth_limit(annotations, type_, max_depth=MAX_TRACK_DEPTH):
    if type_ == 'protein':
        sorted_annots = sorted(annotations, key=lambda a: a['score'])
    else:
        sorted_annots = sorted(annotations, key=lambda a: -a['score'])
    accepted = []
    depth_bins = collections.defaultdict(int)
    for ann in sorted_annots:
        start, end = ann['ref_start'], ann['ref_end']
        max_current_depth = max([depth_bins[pos] for pos in range(start, end+1)], default=0)
        if max_current_depth < max_depth:
            accepted.append(ann)
            for pos in range(start, end+1):
                depth_bins[pos] += 1
    return accepted

def compute_cigar(ref, cons):
    cigar = ''
    count = 0
    op = ''
    length = len(ref)
    for i in range(length):
        r, c = ref[i], cons[i]
        if r == '-' and c != '-':
            this_op = 'D'
        elif r != '-' and c == '-':
            this_op = 'I'
        else:
            this_op = 'M'
        if this_op == op:
            count += 1
        else:
            if op:
                cigar += f"{count}{op}"
            op = this_op
            count = 1
    if op:
        cigar += f"{count}{op}"
    return cigar

def chain_alignments(annotations, max_gap=100, max_overlap=20):
    sorted_ann = sorted(annotations, key=lambda a: (a['name'], a['orient'], a['ref_start']))
    chains = []
    current_chain = None

    for ann in sorted_ann:
        if (not current_chain or
            current_chain['name'] != ann['name'] or
            current_chain['orient'] != ann['orient']):
            if current_chain:
                chains.append(finalize_chain(current_chain))
            current_chain = start_new_chain(ann)
            continue

        ref_gap = ann['ref_start'] - current_chain['last_ref_end']
        cons_gap = ann['cons_start'] - current_chain['last_cons_end']

        colinear = (
            -max_overlap <= ref_gap <= max_gap and
            -max_overlap <= cons_gap <= max_gap and
            ((ann['ref_start'] >= current_chain['last_ref_end'] and ann['cons_start'] >= current_chain['last_cons_end']) or
             (ann['ref_start'] <= current_chain['last_ref_end'] and ann['cons_start'] <= current_chain['last_cons_end']))
        )
        if colinear:
            current_chain['ref_end'] = ann['ref_end']
            current_chain['last_ref_end'] = ann['ref_end']
            current_chain['last_cons_end'] = ann['cons_end']
            oseq = ann.get('ref_seq', '')
            seq = ann.get('cons_seq', '')
            cigar = compute_cigar(seq, oseq) if seq and oseq else ''
            current_chain['components'].append({
                'start': ann['ref_start'],
                'end': ann['ref_end'],
                'ostart': ann['cons_start'],
                'oend': ann['cons_end'],
                'osize': ann['cons_len'],
                'seq': seq.replace('-', ''),
                'oseq': oseq.replace('-', ''),
                'cigar': cigar
            })
        else:
            chains.append(finalize_chain(current_chain))
            current_chain = start_new_chain(ann)
    if current_chain:
        chains.append(finalize_chain(current_chain))
    return chains

def start_new_chain(ann):
    oseq = ann.get('ref_seq', '')
    seq = ann.get('cons_seq', '')
    cigar = compute_cigar(seq, oseq) if seq and oseq else ''
    return {
        'name': ann['name'],
        'orient': ann['orient'],
        'color': ann.get('color', 'gray'),
        'ref_start': ann['ref_start'],
        'ref_end': ann['ref_end'],
        'osize': ann['cons_len'],
        'last_ref_end': ann['ref_end'],
        'last_cons_end': ann['cons_end'],
        'components': [{
            'start': ann['ref_start'],
            'end': ann['ref_end'],
            'ostart': ann['cons_start'],
            'oend': ann['cons_end'],
            'seq': seq.replace('-', ''),
            'oseq': oseq.replace('-', ''),
            'cigar': cigar
        }]
    }

def finalize_chain(chain):
    if not chain:
        return None
    starts = [c['ostart'] for c in chain['components']]
    ends = [c['oend'] for c in chain['components']]
    chain['ostart'] = min(starts)
    chain['oend'] = max(ends)
    return chain


def assign_colors_by_pair(annotations, color_list=DISTINCT_COLORS):
    for i, annot in enumerate(annotations):
        annot["color"] = color_list[i % len(color_list)]
    return annotations


def assign_colors_by_subject(annotations, color_list=DISTINCT_COLORS):
    color_by_subject = {}
    color_idx = 0
    for annot in annotations:
        subject_id = annot['name']
        if subject_id not in color_by_subject:
            color_by_subject[subject_id] = color_list[color_idx % len(color_list)]
            color_idx += 1
        annot['color'] = color_by_subject[subject_id]
    return annotations
