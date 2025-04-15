#!/usr/bin/python3

"""
Extract intron and flank sequences from BAM/SAM files using reference genome FASTA.
Computes insertions, deletions, and mismatches around intron junctions using the cs tag.
"""

import os
import argparse
import pysam
from Bio.Seq import Seq
import re

def validate_files(args):
    """Checks if input files exist and are accessible."""
    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError(f"Error: The reference file '{args.fasta_file}' does not exist.")
    
    if not os.path.exists(args.input_file):
        raise FileNotFoundError(f"Error: The BAM/SAM file '{args.input_file}' does not exist.")
    
    if not args.output_file:
        raise ValueError("Error: Output file (-o) must be specified.")
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        raise FileNotFoundError(f"Error: The output directory '{output_dir}' does not exist.")

def cs_split(cs_tag):
    """
    Split the CS tag into its individual operations (match, insertion, deletion, mismatch, intron)
    with the corresponding length for each operation.

    Example: CTC+a=G*ag=A-g=GC~ct87ac=CCTTGCCCTT 
    will be split into:
                        [(0, 2),(2, 1), (0, 1),(4, 1),(0, 1),(1, 1),(0, 3)] 
    and
                        [(0, 10)]
    """
    token_pattern = re.compile(r"=([ACGTN]+)|\+([acgtn]+)|-([acgtn]+)|\*([acgtn])([acgtn])|~([acgtn]{2})(\d+)([acgtn]{2})")
    
    tokens = list(token_pattern.finditer(cs_tag))
    operations = []

    for token in tokens:
        if token.group(1):  # Match
            operations.append((0, len(token.group(1))))  # 0 for match
        elif token.group(2):  # Insertion
            operations.append((1, len(token.group(2))))  # 1 for insertion
        elif token.group(3):  # Deletion
            operations.append((2, len(token.group(3))))  # 2 for deletion
        elif token.group(4):  # Mismatch
            operations.append((4, 1))  # 4 for mismatch (always 1 base length)
        elif token.group(6):  # Intron
            intron_length = int(token.group(7))
            splice_sites = (token.group(6) + token.group(8)).upper()
            operations.append((3, intron_length))  # 3 for intron
    
    return operations, splice_sites

def parse_cs_tag(cs_tag, window_size):
    """
    Parses the CS tag and counts insertions, deletions, and mismatches in the flanks
    """
    insertions = 0
    deletions = 0
    mismatches = 0
    introns = []
    splice_sites = []

    operations, splice_sites = cs_split(cs_tag)  

    for i, (op, length) in enumerate(operations):
        if op == 3:  
            intron_length = length
            left_bases = []
            right_bases = []
            base_count = 0

            for j in range(i - 1, -1, -1):
                op_type, op_length = operations[j]
                if base_count + op_length > window_size:
                    left_bases.append((op_type, window_size - base_count))
                    break
                left_bases.append((op_type, op_length))
                base_count += op_length

            base_count = 0
            for j in range(i + 1, len(operations)):
                op_type, op_length = operations[j]
                if base_count + op_length > window_size:
                    right_bases.append((op_type, window_size - base_count))
                    break
                right_bases.append((op_type, op_length))
                base_count += op_length

            insertions = sum(length for op, length in left_bases + right_bases if op == 1)
            deletions = sum(length for op, length in left_bases + right_bases if op == 2)
            mismatches = sum(length for op, length in left_bases + right_bases if op == 4)

            introns.append(intron_length)
            return insertions, deletions, mismatches, introns, splice_sites

    return 0, 0, 0, introns, splice_sites

def get_intron_position(intron_pos, read):
    """
    Computes intron start and end positions in reference and query sequences.
    """
    valid_ops_ref = {0, 2, 3}
    valid_ops_query = {0, 1}

    intron_start_ref = read.reference_start
    for op, length in read.cigartuples[:intron_pos]:
        if op in valid_ops_ref:
            intron_start_ref += length

    intron_end_ref = read.reference_start
    for op, length in read.cigartuples[:intron_pos + 1]:
        if op in valid_ops_ref:
            intron_end_ref += length

    intron_start_query = read.query_alignment_start
    for op, length in read.cigartuples[:intron_pos]:
        if op in valid_ops_query:
            intron_start_query += length

    return intron_start_ref + 1, intron_end_ref, intron_start_query

def extract_sequences(read, intron_pos, window_size, reference):
    """
    Extracts the intron, upstream, and downstream sequences.
    """
    intron_start_ref, intron_end_ref, intron_start_query = get_intron_position(intron_pos, read)
    reference_name = read.reference_name
    intron_sequence = reference.fetch(reference_name, intron_start_ref - 1, intron_end_ref)

    # Flank regions
    ups_flank_start = max(0, intron_start_ref - 1)
    ups_flank_end = max(0, intron_start_ref - 1 - window_size)
    down_flank_start = intron_end_ref
    down_flank_end = intron_end_ref + window_size

    upstream_flank = reference.fetch(reference_name, ups_flank_end, ups_flank_start)
    downstream_flank = reference.fetch(reference_name, down_flank_start, down_flank_end)

    return (
        intron_sequence, upstream_flank, downstream_flank,
        ups_flank_start, ups_flank_end, down_flank_start, down_flank_end, 
        intron_start_ref, intron_end_ref, intron_start_query
    )

def analyze_read(read, intron_pos, window_size, reference):
    cs_tag = dict(read.get_tags()).get("cs", "")
    if not cs_tag:
        raise ValueError(f"Error: Read {read.query_name} is missing the CS field.")
    
    intron_start_ref, intron_end_ref, intron_start_query  = get_intron_position(intron_pos, read)
    insertions, deletions, mismatches, introns, splice_sites = parse_cs_tag(cs_tag, window_size)

    return insertions, deletions, mismatches, introns, splice_sites

def analyze_intron_flank(read, intron_pos, window_size, reference):
    """
    Extract sequences and compute statistics for intron regions.
    """
    (
        intron_sequence, 
        upstream_flank, 
        downstream_flank,
        ups_flank_start, 
        ups_flank_end, 
        down_flank_start, 
        down_flank_end,
        intron_start_ref, 
        intron_end_ref, 
        intron_start_query
    ) = extract_sequences(read, intron_pos, window_size, reference)

    insertions, deletions, mismatches, introns, splice_sites = analyze_read(read, intron_pos, window_size, reference)

    return (
        insertions, 
        deletions, 
        mismatches, 
        intron_sequence, 
        upstream_flank, 
        downstream_flank, 
        ups_flank_start, 
        ups_flank_end, 
        down_flank_start, 
        down_flank_end,
        intron_start_ref, 
        intron_end_ref, 
        intron_start_query, 
        introns, 
        splice_sites
    )

def main():
    """
    Main function to parse arguments and process BAM file.
    """
    parser = argparse.ArgumentParser(
        description="Extract intron and flanks sequence, indels, and mismatches in BAM file"
    )
    parser.add_argument(
        "-i", "--input_file", 
        help="Input SAM/BAM file"
    )
    parser.add_argument(
        "-o", "--output_file", 
        help="Output file"
    )
    parser.add_argument(
        "-w", "--window_size", 
        type=int, 
        default=10, 
        help="Window size around the intron (default = 10nt)"
    )
    parser.add_argument(
        "-f", "--fasta_file", 
        help="FASTA file for the reference genome"
    )

    args = parser.parse_args()
    try:
        validate_files(args)
    except FileNotFoundError as e:
        print(e)
        return

    try:
        with pysam.FastaFile(args.fasta_file) as reference:
            with pysam.AlignmentFile(args.input_file, "rb") as bamfile, open(args.output_file, "w") as outfile:
                for read in bamfile:
                    if read.is_unmapped:
                        continue
                    try:
                        cs_tag = dict(read.get_tags()).get("cs", "")
                        if not cs_tag:
                            raise ValueError(f"Error: Read {read.query_name} is missing the cs field.")
                    except ValueError as e:
                        print(e)

                    n_positions = [pos for pos, op in enumerate(read.cigartuples) if op[0] == 3]
                    for intron_pos in n_positions:
                        try:
                            (
                                insertions,
                                deletions,
                                mismatches,
                                intron_sequence,
                                upstream_flank,
                                downstream_flank,
                                ups_flank_start,
                                ups_flank_end,
                                down_flank_start,
                                down_flank_end,
                                intron_start_ref,
                                intron_end_ref,
                                intron_start_query,
                                introns,
                                splice_sites
                            ) = analyze_intron_flank(read, intron_pos, args.window_size, reference)

                            orientation_alignment = "Forward" if not read.is_reverse else "Reverse"

                            outfile.write(
                                f">Intron,{read.reference_name},"
                                f"RefPos:{intron_start_ref}_{intron_end_ref},"
                                f"{read.query_name},{orientation_alignment}," 
                                f"I:{insertions},"
                                f"D:{deletions},"
                                f"X:{mismatches}," 
                                f"QueryPos:{intron_start_query},"
                                f"SpliceSite:{splice_sites},"
                                f"Intron_Length:{introns}\n"
                            )

                            outfile.write(f"{str(intron_sequence)}\n")

                            outfile.write(
                                f">Downstream_Flank,{read.reference_name},{down_flank_start+1}_{down_flank_end+1},{read.query_name},{orientation_alignment}\n"
                                f"{downstream_flank}\n"
                            )

                            outfile.write(
                                f">Upstream_Flank,{read.reference_name},{ups_flank_start}_{ups_flank_end},{read.query_name},{orientation_alignment}\n"
                                f"{upstream_flank}\n\n"
                            )
                        except Exception as e:
                            print(f"Error processing intron in read {read.query_name}: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    main()
