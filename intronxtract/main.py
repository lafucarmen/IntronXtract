import os
import sys
import argparse
import pysam
from .stats import StatsCollector
from .bam_processing import extract_flanks_and_intron_sequence
from .analysis import count_transcriptomic_support, remove_redundant_introns
from .sequence_utils import duplicate_intron_data
from .output_writer import (
    write_gff, write_main_fasta_output, write_iic_file, 
    write_correct_reads_list, write_stats_summary
)

def create_parser():
    """Creates the argument parser."""
    parser = argparse.ArgumentParser(description="Extract intron and flanks sequence, indels, and mismatches in BAM file")
    parser.add_argument("-i", "--input_file", help="Input SAM/BAM file", required=True)
    parser.add_argument("-o", "--output_file", help="Output file for results", required=True)
    parser.add_argument("-w", "--window_size", type=int, default=10, help="Window size around the intron (default = 10nt)")
    parser.add_argument("-f", "--fasta_file", help="FASTA file", required=True)
    parser.add_argument("--filter_indel_free", action="store_true", help="Extract only introns without insertions or deletions")
    parser.add_argument("--transcriptomic_support", type=int, help="Minimum number of metaT reads supporting each intron")
    parser.add_argument("--remove_redundancy", action="store_true", help="Remove redundant introns based on position and identity")
    parser.add_argument("--intron_iic", action="store_true", help="Add additional iic file")
    parser.add_argument("--duplicate", action="store_true", help="Duplicate intron sequences and flanks to have both strands (useful when we don't know the orientation of the read)")
    parser.add_argument("--stats",action="store_true",help="Add metaT and intron statistics")
    return parser

def main():
    """Main function that orchestrates the intron extraction."""
    parser = create_parser()
    args = parser.parse_args()

    # 1. Initialize statistics collector
    stats = StatsCollector()
    introns_data = []

    # 2. BAM Processing
    print(f"Starting intron extraction from {args.input_file}...")
    try:
        with pysam.FastaFile(args.fasta_file) as reference:
            with pysam.AlignmentFile(args.input_file, "rb") as bamfile:
                for read in bamfile:
                    stats.total_metaT.add(read.query_name)
                    if read.reference_name:
                         stats.total_metaG.add(read.reference_name)
                    
                    if read.is_unmapped or not read.reference_name:
                        continue
                    
                    n_positions = [pos for pos, op in enumerate(read.cigartuples) if op[0] == 3]
                    for n_pos in n_positions:
                        result = extract_flanks_and_intron_sequence(
                            read, n_pos, args.window_size, reference, 
                            stats, args.filter_indel_free
                        )
                        if result is not None:
                            introns_data.append(result)
    except FileNotFoundError as e:
        print(f"Error: File not found. {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"Error processing BAM/FASTA file: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        return 1

    print(f"Found {stats.total_introns} total introns ({stats.noindels_introns} passing filters).")

    # 3. Analysis and Filtering
    introns_data = count_transcriptomic_support(introns_data)

    if args.remove_redundancy:
        introns_data, _ = remove_redundant_introns(introns_data)
        stats.noindels_nonredundant_introns = len(introns_data)
        stats.noindels_noredundant_metaT = set(i['intron_info']['query_name'] for i in introns_data)
        stats.noindels_noredundant_metaG = set(i['intron_info']['reference_name'] for i in introns_data)
        print(f"Filtered down to {len(introns_data)} non-redundant introns.")

    # Save final data before support filter (for stats)
    introns_data_final_for_stats = introns_data.copy()

    if args.transcriptomic_support:
        introns_data = [i for i in introns_data if i['transcriptomic_support'] >= args.transcriptomic_support]
        stats.noindels_nonredundant_supported_introns = len(introns_data)
        stats.noindels_noredundant_supported_metaT = set(i['intron_info']['query_name'] for i in introns_data)
        stats.noindels_noredundant_supported_metaG = set(i['intron_info']['reference_name'] for i in introns_data)
        print(f"Filtered down to {len(introns_data)} supported introns (support >= {args.transcriptomic_support}).")
        # Use supported introns for splice site stats if filter is applied
        introns_data_final_for_stats = introns_data.copy()


    # 4. Write statistics
    if args.stats:
        print("Writing statistics files...")
        base_output = os.path.splitext(args.output_file)[0]
        stats_output_file = base_output + "_stats"
        correct_reads_output_file = base_output + "_correct_structure_reads"
        gff_output_file = base_output + "_correct_intron_exon_structures.gff"

        correct_structure_metaT = list(stats.total_spliced_metaT - stats.indels_spliced_metaT)
        bam_filename = os.path.basename(args.input_file)
        
        write_stats_summary(stats_output_file, stats, introns_data_final_for_stats, bam_filename, correct_structure_metaT)
        write_correct_reads_list(correct_reads_output_file, correct_structure_metaT)
        write_gff(args.input_file, gff_output_file, set(correct_structure_metaT))

    # 5. Duplication (optional)
    if args.duplicate:
        print("Duplicating introns for reverse complement...")
        duplicated_introns = []
        for intron in introns_data:
            duplicated_introns.append(intron)
            duplicated_introns.append(duplicate_intron_data(intron))
        introns_data = duplicated_introns

    # 6. Write main outputs
    print(f"Writing main output to {args.output_file}...")
    write_main_fasta_output(args.output_file, introns_data, args.duplicate)
    
    if args.intron_iic:
        iic_output_file = os.path.splitext(args.output_file)[0] + ".iic"
        print(f"Writing IIC file to {iic_output_file}...")
        write_iic_file(iic_output_file, introns_data, args.duplicate)

    print("Intron extraction complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
