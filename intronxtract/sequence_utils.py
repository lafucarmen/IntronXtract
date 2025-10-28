from Bio.Seq import Seq

def extract_splice_site(intron_sequence):
    donor_splice=intron_sequence[:2].upper()
    acceptor_splice=intron_sequence[-2:].upper()
    return f"{donor_splice}-{acceptor_splice}"

def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.
    """
    return str(Seq(sequence).reverse_complement())

def duplicate_intron_data(intron_data):
    """
    Duplicates an intron's data using the reverse complement of sequences.
    """
    # Create a copy of the dictorionary containing intron info
    duplicated_data = intron_data.copy()
    
    # Add a new tag to distinguish the copied intron
    duplicated_data["intron_info"] = duplicated_data["intron_info"].copy()

    # Do the reverse complement of intron sequence and flanks
    duplicated_data["intron_info"]["intron_sequence"] = reverse_complement(intron_data["intron_info"]["intron_sequence"])
    duplicated_data["intron_info"]["ups_flank"] = reverse_complement(intron_data["intron_info"]["down_flank"])
    duplicated_data["intron_info"]["down_flank"] = reverse_complement(intron_data["intron_info"]["ups_flank"])
    duplicated_data["intron_info"]["ups_flank_200"] = reverse_complement(intron_data["intron_info"]["down_flank_200"])
    duplicated_data["intron_info"]["down_flank_200"] = reverse_complement(intron_data["intron_info"]["ups_flank_200"])
    
    # Extract correct splice site
    duplicated_data["intron_info"]["splice_site"] = extract_splice_site(duplicated_data["intron_info"]["intron_sequence"])

    # Invert the orientation of the alignment
    duplicated_data["intron_info"]["orientation_alignment"] = (
        "Reverse" if intron_data["intron_info"]["orientation_alignment"] == "Forward" else "Forward"
    )
    duplicated_data['is_duplicate'] = True
    
    return duplicated_data
