from .sequence_utils import extract_splice_site

def count_flank_events(read, ref_range):
    aligned = read.get_aligned_pairs(matches_only=False, with_seq=True, with_cigar=True)
    mismatches = 0
    insertions = 0
    deletions = 0

    mismatch_tuples = []
    match_tuples = []
    insertion_tuples = []
    deletion_tuples = []

    for i, (q_pos, r_pos, base, cigar) in enumerate(aligned):
        # Mismatches and matches (inside range)
        if r_pos is not None and r_pos in ref_range:
            if base:
                if base.islower():  # mismatch
                    mismatches += 1
                    mismatch_tuples.append((q_pos, r_pos, base, cigar))
                elif base.isupper() and cigar == 0:  # match
                    match_tuples.append((q_pos, r_pos, base, cigar))

            if q_pos is None and cigar == 2:  # deletion
                deletions += 1
                deletion_tuples.append((q_pos, r_pos, base, cigar))

        # Insertions (r_pos is None, but adjacent to ref_range)
        elif r_pos is None and cigar == 1:
            left = aligned[i - 1] if i > 0 else None
            right = aligned[i + 1] if i + 1 < len(aligned) else None
            left_in_range = left and left[1] in ref_range
            right_in_range = right and right[1] in ref_range

            if left_in_range or right_in_range:
                insertions += 1
                insertion_tuples.append((q_pos, r_pos, base, cigar))

    return {
        'mismatches': mismatches,
        'insertions': insertions,
        'deletions': deletions
    }


def calculate_identity(read):
    identity_tag = read.get_tag("XI:f:")
    identity = float(identity_tag)
    return identity

def calculate_intron_position(n_pos, read):
    authorized_ops_ref = {0, 2, 3}
    authorized_ops_query = {0, 1}
    # only match (0), deletion(2) and skipped region (3) consume reference
    # match (0), insertion (1) and softclipped (4) consume query, but S are already taken into account for query start

    # Calculate intron start position on the reference
    intron_start_ref = read.reference_start
    for op, length in read.cigartuples[:n_pos] :
        if op in authorized_ops_ref:
            intron_start_ref += length

    # Calculate intron end position on the reference
    intron_end_ref = read.reference_start
    for op, length in read.cigartuples[:n_pos + 1]:
        if op in authorized_ops_ref:
            intron_end_ref += length

    # Calculate intron start position on the query
    soft_clipped_bases = sum(length for op, length in read.cigartuples if op == 4)
    intron_start_query = read.query_alignment_start
    for op, length in read.cigartuples[:n_pos]:
        if op in authorized_ops_query:
            intron_start_query += length

    return intron_start_ref+1, intron_end_ref, intron_start_query


def extract_flanks_and_intron_sequence(read, n_pos, window_size, reference, stats_collector, filter_indel_free):
    
    intron_start_ref, intron_end_ref, intron_start_query = calculate_intron_position(n_pos, read)

    identity = calculate_identity(read)
    reference_name = read.reference_name
    intron_sequence = reference.fetch(reference_name, intron_start_ref - 1, intron_end_ref)
    intron_length = len(intron_sequence)

    ups_flank_start = max(0, intron_start_ref - 1)
    ups_flank_end = max(0, intron_start_ref - 1 - window_size)
    down_flank_start = intron_end_ref + 1
    down_flank_end = intron_end_ref + window_size +1

    ups_flank = reference.fetch(reference_name, ups_flank_end, ups_flank_start)
    down_flank = reference.fetch(reference_name, down_flank_start, down_flank_end)

    ups_flank_200 = reference.fetch(reference_name, max(0, intron_start_ref - 200), intron_start_ref - 1)
    down_flank_200 = reference.fetch(reference_name, intron_end_ref, intron_end_ref + 200)

    ups_range = range(min(ups_flank_end, ups_flank_start), max(ups_flank_end, ups_flank_start))
    downs_range = range(down_flank_start, down_flank_end)

    ups_events = count_flank_events(read, ups_range)
    downs_events = count_flank_events(read, downs_range)

    mismatches = ups_events["mismatches"] + downs_events["mismatches"]
    insertions = ups_events["insertions"] + downs_events["insertions"]
    deletions = ups_events["deletions"] + downs_events["deletions"]

    orientation_alignment = "Forward" if not read.is_reverse else "Reverse"

    if filter_indel_free and (insertions > 0 or deletions > 0 or mismatches > 0):
        # Report this as an indel-containing intron to stats
        stats_collector.add_spliced_read(read, has_indel=True)
        return None 

    # Report this as a clean (no-indel) intron to stats
    stats_collector.add_spliced_read(read, has_indel=False)

    return {
        "intron_info": {
            "reference_name": reference_name,
            "query_name": read.query_name,
            "orientation_alignment": orientation_alignment,
            "intron_length": intron_length,
            "insertions": insertions,
            "deletions": deletions,
            "identity": identity,
            "mismatches": mismatches,
            "intron_sequence": intron_sequence,
            "ups_flank": ups_flank,
            "down_flank": down_flank,
            "ups_flank_200": ups_flank_200,
            "down_flank_200": down_flank_200,
            "splice_site": extract_splice_site(intron_sequence),
            "down_flank_start": down_flank_start + 1,
            "down_flank_end": down_flank_end,
            "ups_flank_start": ups_flank_start +1,
            "ups_flank_end": ups_flank_end,
            "intron_start_ref": intron_start_ref,
            "intron_end_ref": intron_end_ref,
            "intron_start_query": intron_start_query
        },
        "transcriptomic_support": 1
    }
