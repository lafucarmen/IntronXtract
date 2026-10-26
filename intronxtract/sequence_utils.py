_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence"""
    return sequence.translate(_COMPLEMENT)[::-1]


def extract_splice_site(intron_sequence):
    """Return the dinucleotide splice site of an intron"""
    donor_splice = intron_sequence[:2].upper()
    acceptor_splice = intron_sequence[-2:].upper()
    return f"{donor_splice}-{acceptor_splice}"


def duplicate_intron_data(intron_data):
    """
    Return a copy of an intron record on the opposite strand
    (reverse complement of the intron and swapped/reverse-complemented flanks).
    """
    info = intron_data["intron_info"]
    duplicated_data = intron_data.copy()
    dup = duplicated_data["intron_info"] = info.copy()

    dup["intron_sequence"] = reverse_complement(info["intron_sequence"])
    dup["ups_flank"] = reverse_complement(info["down_flank"])
    dup["down_flank"] = reverse_complement(info["ups_flank"])
    dup["ups_flank_200"] = reverse_complement(info["down_flank_200"])
    dup["down_flank_200"] = reverse_complement(info["ups_flank_200"])
    dup["splice_site"] = extract_splice_site(dup["intron_sequence"])
    dup["orientation_alignment"] = (
        "Reverse" if info["orientation_alignment"] == "Forward" else "Forward"
    )
    duplicated_data["is_duplicate"] = True
    return duplicated_data