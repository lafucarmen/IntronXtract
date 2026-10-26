"""Command-line interface of intronxtract."""

import argparse
import logging
import os
import sys
from collections import Counter

import pysam

from . import __version__
from .stats import StatsCollector
from .bam_processing import (
    extract_flanks_and_intron_sequence, get_aligned_events, get_identity, get_alignment_score,
    keep_alignment, alignment_type, CIGAR_N,
)
from .analysis import count_transcriptomic_support, remove_redundant_introns
from .sequence_utils import duplicate_intron_data
from .output_writer import (
    open_alignment_file, write_gff, write_main_fasta_output, write_iic_file,
    write_correct_reads_list, write_stats_summary, write_sl_output,
)

logger = logging.getLogger("intronxtract")


class InputError(Exception):
    """Problem with the user's input files, reported without a traceback."""


def create_parser():
    """Create the argument parser."""
    parser = argparse.ArgumentParser(
        prog="intronxtract",
        description=(
            "Extract intron sequences, their flanking regions and alignment "
            "events (mismatches/indels in the flanks) from spliced RNA-seq "
            "alignments (SAM/BAM/CRAM) and a reference FASTA."
        ),
        epilog="Documentation: https://github.com/lafucarmen/IntronXtract",
    )
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}")

    req = parser.add_argument_group("required arguments")
    req.add_argument("-i", "--input", "--input_file", dest="input_file", required=True, metavar="ALIGNMENTS",
                     help="Spliced alignments in SAM, BAM or CRAM format (e.g. minimap2 -ax splice)")
    req.add_argument("-f", "--fasta", "--fasta_file", dest="fasta_file", required=True, metavar="FASTA",
                     help="Reference FASTA used for the alignment (plain or bgzip-compressed)")
    req.add_argument("-o", "--output", "--output_file", dest="output_file", required=True, metavar="OUTPUT",
                     help="Main output file (FASTA-like). Other outputs use the same prefix")

    opt = parser.add_argument_group("extraction")
    opt.add_argument("-w", "--window", "--window_size", dest="window_size", type=int, default=10, metavar="N",
                     help="Length (nt) of the exonic flanks reported and checked for indels/mismatches "
                          "(default: %(default)s)")
    opt.add_argument("--identity-tag", default="XI", metavar="TAG",
                     help="SAM tag holding the read identity (default: %(default)s). "
                          "If absent, identity is computed from the NM tag")

    flt = parser.add_argument_group("filters")
    flt.add_argument("--filter-indel-free", "--filter_indel_free", dest="filter_indel_free",
                     action="store_true",
                     help="Keep only introns with no insertions, deletions or mismatches in the flanks")
    flt.add_argument("--min-support", "--transcriptomic_support", dest="transcriptomic_support",
                     type=int, metavar="MIN_READS",
                     help="Keep only introns supported by at least MIN_READS alignments")
    flt.add_argument("--remove-redundancy", "--remove_redundancy", dest="remove_redundancy",
                     action="store_true",
                     help="Report each intron once: the read with the highest identity "
                          "(or highest AS if there is no identity, otherwise one at random)")

    out = parser.add_argument_group("additional outputs")
    out.add_argument("--duplicate", action="store_true",
                     help="Also write the reverse complement of each intron and flanks "
                          "(useful when the read strand is unknown)")
    out.add_argument("--intron-iic", "--intron_iic", dest="intron_iic", action="store_true",
                     help="Write PREFIX.iic: 200 nt upstream, intron, 200 nt downstream (tab-separated)")
    out.add_argument("--stats", action="store_true",
                     help="Write PREFIX_stats, PREFIX_correct_structure_reads and a GFF3 of read structures")
    out.add_argument("--sl-output", "--sl_output", dest="sl_output", action="store_true",
                     help="Write PREFIX_sl_output.tsv with 5' soft-clipped sequences "
                          "(spliced-leader / trans-splicing detection)")

    misc = parser.add_argument_group("misc")
    misc.add_argument("-q", "--quiet", action="store_true", help="Only print warnings and errors")
    return parser


def validate_inputs(args):
    """Check input files and options before starting; raise InputError with a clear message."""
    for label, path in (("Alignment file", args.input_file), ("FASTA file", args.fasta_file)):
        if not os.path.isfile(path):
            raise InputError(f"{label} not found: {path}")
    if args.window_size < 0:
        raise InputError("--window must be >= 0")
    if args.transcriptomic_support is not None and args.transcriptomic_support < 1:
        raise InputError("--min-support must be >= 1")
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    if not os.path.isdir(out_dir):
        raise InputError(f"Output directory does not exist: {out_dir}")


def check_reference_names(bamfile, reference):
    """Make sure the alignment and the FASTA describe the same sequences."""
    bam_refs = set(bamfile.references)
    fasta_refs = set(reference.references)
    if bam_refs and not bam_refs & fasta_refs:
        raise InputError(
            "None of the sequence names in the alignment header are present in the FASTA. "
            "Did you use the same reference for the alignment? "
            f"(e.g. alignment: '{sorted(bam_refs)[0]}', FASTA: '{sorted(fasta_refs)[0] if fasta_refs else 'empty'}')"
        )
    missing = bam_refs - fasta_refs
    if missing:
        logger.warning("%d sequence(s) of the alignment header are not in the FASTA; "
                       "reads on them will be skipped (e.g. '%s')", len(missing), sorted(missing)[0])
    return fasta_refs


def extract_introns(args, stats):
    """Read the alignment file and return the list of intron records."""
    introns_data = []
    skipped_no_seq = 0
    used = Counter()  # mapped alignments used, by type

    with pysam.FastaFile(args.fasta_file) as reference, \
            open_alignment_file(args.input_file, args.fasta_file) as bamfile:
        fasta_refs = check_reference_names(bamfile, reference)

        for read in bamfile.fetch(until_eof=True):
            stats.total_metaT.add(read.query_name)
            if read.reference_name:
                stats.total_metaG.add(read.reference_name)

            if not keep_alignment(read):
                continue
            if read.reference_name not in fasta_refs:
                continue
            used[alignment_type(read)] += 1

            n_positions = [pos for pos, (op, _) in enumerate(read.cigartuples) if op == CIGAR_N]
            if not n_positions:
                continue

            events = get_aligned_events(read, reference)
            if events is None:
                skipped_no_seq += 1
                continue
            identity = get_identity(read, args.identity_tag)
            alignment_score = get_alignment_score(read)

            for n_pos in n_positions:
                result = extract_flanks_and_intron_sequence(
                    read, n_pos, args.window_size, reference, stats,
                    args.filter_indel_free, events, identity, alignment_score,
                )
                if result is not None:
                    introns_data.append(result)

    report_alignment_types(used)
    if skipped_no_seq:
        logger.warning("%d spliced alignment(s) skipped: no read sequence (SEQ='*') and no MD tag",
                       skipped_no_seq)
    return introns_data


def report_alignment_types(used):
    """Tell the user which kinds of mapped alignments were used."""
    message = ("Mapped alignments used: %d primary, %d secondary, %d supplementary."
               % (used["primary"], used["secondary"], used["supplementary"]))
    if used["secondary"] or used["supplementary"]:
        logger.warning("%s Secondary and supplementary alignments are included. If you do not "
                       "want them, remove them from the alignment file first, e.g. "
                       "'samtools view -b -F 0x900 in.bam > primary.bam' (0x100 = secondary, "
                       "0x800 = supplementary).", message)
    else:
        logger.info("%s", message)


def run(args):
    validate_inputs(args)
    stats = StatsCollector()
    prefix = os.path.splitext(args.output_file)[0]

    logger.info("intronxtract %s: reading %s", __version__, args.input_file)
    introns_data = extract_introns(args, stats)
    logger.info("Found %d introns (%d without indels/mismatches in the flanks).",
                stats.total_introns, stats.noindels_introns)

    # Analysis and filtering
    introns_data = count_transcriptomic_support(introns_data)

    if args.remove_redundancy:
        introns_data, _ = remove_redundant_introns(introns_data)
        stats.noindels_nonredundant_introns = len(introns_data)
        stats.noindels_noredundant_metaT = {i["intron_info"]["query_name"] for i in introns_data}
        stats.noindels_noredundant_metaG = {i["intron_info"]["reference_name"] for i in introns_data}
        logger.info("%d non-redundant introns.", len(introns_data))

    if args.transcriptomic_support:
        introns_data = [i for i in introns_data if i["transcriptomic_support"] >= args.transcriptomic_support]
        stats.noindels_nonredundant_supported_introns = len(introns_data)
        stats.noindels_noredundant_supported_metaT = {i["intron_info"]["query_name"] for i in introns_data}
        stats.noindels_noredundant_supported_metaG = {i["intron_info"]["reference_name"] for i in introns_data}
        logger.info("%d introns with support >= %d.", len(introns_data), args.transcriptomic_support)

    if args.stats:
        stats_output_file = prefix + "_stats"
        correct_reads_output_file = prefix + "_correct_structure_reads"
        gff_output_file = prefix + "_correct_intron_exon_structures.gff"
        logger.info("Writing statistics to %s", stats_output_file)

        correct_structure_metaT = list(stats.total_spliced_metaT - stats.indels_spliced_metaT)
        write_stats_summary(stats_output_file, stats, introns_data,
                            os.path.basename(args.input_file), correct_structure_metaT)
        write_correct_reads_list(correct_reads_output_file, correct_structure_metaT)
        write_gff(args.input_file, gff_output_file, set(correct_structure_metaT), args.fasta_file)

    # Duplication: splice sites are not oriented when the read strand is unknown
    if args.duplicate:
        duplicated = []
        for intron in introns_data:
            duplicated.append(intron)
            duplicated.append(duplicate_intron_data(intron))
        introns_data = duplicated

    logger.info("Writing main output to %s", args.output_file)
    write_main_fasta_output(args.output_file, introns_data, args.duplicate)

    if args.intron_iic:
        iic_output_file = prefix + ".iic"
        logger.info("Writing IIC file to %s", iic_output_file)
        write_iic_file(iic_output_file, introns_data, args.duplicate)

    if args.sl_output:
        sl_output_file = prefix + "_sl_output.tsv"
        logger.info("Writing SL output to %s", sl_output_file)
        write_sl_output(args.input_file, sl_output_file, args.fasta_file, args.identity_tag)

    logger.info("Done.")
    return 0


def main(argv=None):
    """Entry point of the `intronxtract` command."""
    args = create_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.WARNING if args.quiet else logging.INFO,
        format="[intronxtract] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )
    try:
        return run(args)
    except InputError as e:
        logger.error("%s", e)
        return 1
    except (OSError, ValueError, KeyError) as e:
        logger.error("Could not process the input files: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())