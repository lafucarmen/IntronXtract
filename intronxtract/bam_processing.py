from .sequence_utils import extract_splice_site, reverse_complement
from .sequence_utils import reverse_complement

CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_P, CIGAR_EQ, CIGAR_X = range(9)

# Operations consuming the reference: M, D, N, =, X
REF_CONSUMING = {CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X}
# Operations consuming the query, excluding soft clips (already accounted
# for by read.query_alignment_start): M, I, =, X
QUERY_CONSUMING = {CIGAR_M, CIGAR_I, CIGAR_EQ, CIGAR_X}
# Alignment columns used to compute identity from NM: M, I, D, =, X
ALIGNMENT_COLUMNS = {CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_EQ, CIGAR_X}

def keep_alignment(read):
    """
    Return True if this alignment record should be processed.
    Only unmapped records are skipped: primary, secondary and supplementary
    alignments are all used (filter the BAM beforehand to exclude them).
    """
    return not read.is_unmapped and bool(read.reference_name) and bool(read.cigartuples)

def alignment_type(read):
    """Return 'primary', 'secondary' or 'supplementary'."""
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"

def get_identity(read, identity_tag="XI"):
    """
    Return the alignment identity of a read.
    Uses ``identity_tag`` when present in the record. Otherwise it is computed
    from the NM tag as 100 * (1 - NM / alignment_columns), where alignment
    columns are M/I/D/=/X operations. Returns None if neither is available.
    """
    if identity_tag and read.has_tag(identity_tag):
        return float(read.get_tag(identity_tag))
    if read.has_tag("NM"):
        columns = sum(length for op, length in read.cigartuples if op in ALIGNMENT_COLUMNS)
        if columns:
            return 100.0 * (1 - read.get_tag("NM") / columns)
    return None


def get_alignment_score(read):
    """Return the AS tag of a read, or None if absent."""
    return read.get_tag("AS") if read.has_tag("AS") else None

def _md_mismatch_positions(read):
    """
    Return the set of 0-based reference positions that are mismatches according
    to the MD tag. Works without the read sequence (e.g. secondary alignments
    written with SEQ='*').
    """
    # Reference positions of aligned (M/=/X) columns, in order. MD describes
    # these columns plus deletions; introns (N) are not part of MD.
    match_columns = []
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (CIGAR_M, CIGAR_EQ, CIGAR_X):
            match_columns.extend(range(ref_pos, ref_pos + length))
        if op in REF_CONSUMING:
            ref_pos += length

    md = read.get_tag("MD")
    mismatches = set()
    column = 0
    i = 0
    while i < len(md):
        if md[i].isdigit():
            j = i
            while j < len(md) and md[j].isdigit():
                j += 1
            column += int(md[i:j])          # matching bases
            i = j
        elif md[i] == "^":
            i += 1
            while i < len(md) and md[i].isalpha():
                i += 1                      # deleted reference bases: not match columns
        else:
            if column < len(match_columns):
                mismatches.add(match_columns[column])
            column += 1                     # one mismatched base
            i += 1
    return mismatches


def get_aligned_events(read, reference):
    """
    Return a list of (query_pos, ref_pos, cigar_op, is_mismatch) for a read,
    or None if mismatches cannot be determined.

    Mismatches come from the MD tag when present (so the read sequence is not
    needed, e.g. for secondary alignments with SEQ='*'). Without MD, the read
    sequence is compared to the reference FASTA. Without MD and without
    sequence, None is returned.
    """
    pairs = read.get_aligned_pairs(matches_only=False, with_cigar=True)

    if read.has_tag("MD"):
        mismatch_positions = _md_mismatch_positions(read)
        return [
            (q_pos, r_pos, op, r_pos in mismatch_positions and op != CIGAR_D)
            for q_pos, r_pos, op in pairs
        ]

    query_seq = read.query_sequence
    if query_seq is None:  # e.g. secondary alignments without SEQ
        return None

    offset = read.reference_start
    ref_seq = reference.fetch(read.reference_name, offset, read.reference_end).upper()
    events = []
    for q_pos, r_pos, op in pairs:
        is_mismatch = False
        if q_pos is not None and r_pos is not None:
            if op == CIGAR_X:
                is_mismatch = True
            elif op == CIGAR_M:
                is_mismatch = query_seq[q_pos].upper() != ref_seq[r_pos - offset]
        events.append((q_pos, r_pos, op, is_mismatch))
    return events

def count_flank_events(events, ref_range):
    """Count mismatches, insertions and deletions inside a reference range (0-based)."""
    mismatches = insertions = deletions = 0

    for i, (q_pos, r_pos, op, is_mismatch) in enumerate(events):
        if r_pos is not None and r_pos in ref_range:
            if is_mismatch:
                mismatches += 1
            if q_pos is None and op == CIGAR_D:
                deletions += 1
        # Insertions have no reference position: count them if adjacent to the range
        elif r_pos is None and op == CIGAR_I:
            left = events[i - 1] if i > 0 else None
            right = events[i + 1] if i + 1 < len(events) else None
            if (left and left[1] in ref_range) or (right and right[1] in ref_range):
                insertions += 1

    return {"mismatches": mismatches, "insertions": insertions, "deletions": deletions}

def calculate_intron_position(n_pos, read):
    """
    Return (intron_start_ref, intron_end_ref, intron_start_query) for the N
    operation at index ``n_pos`` of the CIGAR.
    Reference coordinates are 1-based. The query position is the
    0-based offset in the read sequence as stored in the BAM.
    """
    cigar = read.cigartuples

    intron_start_ref = read.reference_start
    for op, length in cigar[:n_pos]:
        if op in REF_CONSUMING:
            intron_start_ref += length
    intron_end_ref = intron_start_ref + cigar[n_pos][1]

    intron_start_query = read.query_alignment_start
    for op, length in cigar[:n_pos]:
        if op in QUERY_CONSUMING:
            intron_start_query += length

    return intron_start_ref + 1, intron_end_ref, intron_start_query


def extract_flanks_and_intron_sequence(read, n_pos, window_size, reference, stats_collector,
                                       filter_indel_free, events, identity, alignment_score=None):
    """Build the record for one intron of one read, or None if it is filtered out."""
    intron_start_ref, intron_end_ref, intron_start_query = calculate_intron_position(n_pos, read)

    reference_name = read.reference_name
    contig_length = reference.get_reference_length(reference_name)
    intron_sequence = reference.fetch(reference_name, intron_start_ref - 1, intron_end_ref)

    # Flanks as 0-based half-open intervals
    ups_start0 = max(0, intron_start_ref - 1 - window_size)
    ups_end0 = intron_start_ref - 1
    down_start0 = intron_end_ref
    down_end0 = min(contig_length, intron_end_ref + window_size)

    ups_flank = reference.fetch(reference_name, ups_start0, ups_end0)
    down_flank = reference.fetch(reference_name, down_start0, down_end0)
    ups_flank_200 = reference.fetch(reference_name, max(0, intron_start_ref - 201), intron_start_ref - 1)
    down_flank_200 = reference.fetch(reference_name, intron_end_ref, min(contig_length, intron_end_ref + 200))

    ups_events = count_flank_events(events, range(ups_start0, ups_end0))
    down_events = count_flank_events(events, range(down_start0, down_end0))
    mismatches = ups_events["mismatches"] + down_events["mismatches"]
    insertions = ups_events["insertions"] + down_events["insertions"]
    deletions = ups_events["deletions"] + down_events["deletions"]

    has_events = insertions > 0 or deletions > 0 or mismatches > 0
    stats_collector.add_spliced_read(read, has_indel=has_events)
    if filter_indel_free and has_events:
        return None

    return {
        "intron_info": {
            "reference_name": reference_name,
            "query_name": read.query_name,
            "orientation_alignment": "Reverse" if read.is_reverse else "Forward",
            "intron_length": len(intron_sequence),
            "insertions": insertions,
            "deletions": deletions,
            "identity": identity,
            "alignment_score": alignment_score,
            "mismatches": mismatches,
            "intron_sequence": intron_sequence,
            "ups_flank": ups_flank,
            "down_flank": down_flank,
            "ups_flank_200": ups_flank_200,
            "down_flank_200": down_flank_200,
            "splice_site": extract_splice_site(intron_sequence),
            # 1-based, inclusive coordinates
            "ups_flank_start": ups_start0 + 1,
            "ups_flank_end": ups_end0,
            "down_flank_start": down_start0 + 1,
            "down_flank_end": down_end0,
            "intron_start_ref": intron_start_ref,
            "intron_end_ref": intron_end_ref,
            "intron_start_query": intron_start_query,
        },
        "transcriptomic_support": 1,
    }

def extract_unmapped_region(read, min_softclip=10, mapped_bases=5):
    """
    Return (seq, qstart, qend) for the 5' soft-clipped region plus the first
    ``mapped_bases`` mapped nucleotides, in the read's original orientation.
    Only reported if the soft clip is at least ``min_softclip`` nt.
    """
    seq = read.query_sequence
    cigar = read.cigartuples

    if read.is_reverse:
        seq = reverse_complement(seq)
        cigar = cigar[::-1]

    if cigar and cigar[0][0] == CIGAR_S:
        softclip_len = cigar[0][1]
        if softclip_len >= min_softclip:
            end = softclip_len + mapped_bases
            return seq[:end], 1, end  # 1-based coordinates of the trimmed sequence
    return None, None, None