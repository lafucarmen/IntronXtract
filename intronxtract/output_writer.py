"""Writers for all output files."""

import pysam

from .stats import group_splice_site
from .sequence_utils import extract_splice_site
from .bam_processing import (
    extract_unmapped_region, keep_alignment, get_identity, get_alignment_score,
    CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X,
)

def open_alignment_file(path, fasta_file=None):
    """Open a SAM/BAM/CRAM file (format auto-detected)."""
    return pysam.AlignmentFile(path, reference_filename=fasta_file)

def write_gff(bam_file, output_file, correct_structures, fasta_file=None):
    """
    Write a GFF3 file with the exon/intron structure of the reads whose introns
    have no indels/mismatches in their flanks. The BAM does not need to be indexed.
    """
    with open_alignment_file(bam_file, fasta_file) as bam, open(output_file, "w") as gff_file:
        gff_file.write("##gff-version 3\n")
        mrna_counter = 0

        for read in bam.fetch(until_eof=True):
            if read.query_name not in correct_structures:
                continue
            if not keep_alignment(read):
                continue

            reference_name = read.reference_name
            start_ref = read.reference_start + 1
            end_ref = read.reference_end
            strand = "."
            mrna_counter += 1
            mrna_id = f"aln{mrna_counter}"

            gff_file.write(
                f"{reference_name}\tintronxtract\tmRNA\t{start_ref}\t{end_ref}\t.\t{strand}\t.\t"
                f"ID={mrna_id};Name={read.query_name}\n"
            )

            exon_positions = []
            intron_positions = []
            exon_start = None
            current_pos = start_ref

            for op, length in read.cigartuples:
                if op in (CIGAR_M, CIGAR_D, CIGAR_EQ, CIGAR_X):  # exon
                    if exon_start is None:
                        exon_start = current_pos
                    current_pos += length
                elif op == CIGAR_N:  # intron
                    if exon_start is not None:
                        exon_positions.append((exon_start, current_pos - 1))
                    intron_positions.append((current_pos, current_pos + length - 1))
                    current_pos += length
                    exon_start = None

            if exon_start is not None:
                exon_positions.append((exon_start, current_pos - 1))

            for n, (exon_start, exon_end) in enumerate(exon_positions, 1):
                gff_file.write(
                    f"{reference_name}\tintronxtract\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t"
                    f"ID={mrna_id}.exon{n};Parent={mrna_id}\n"
                )
            for n, (intron_start, intron_end) in enumerate(intron_positions, 1):
                gff_file.write(
                    f"{reference_name}\tintronxtract\tintron\t{intron_start}\t{intron_end}\t.\t{strand}\t.\t"
                    f"ID={mrna_id}.intron{n};Parent={mrna_id}\n"
                )


def _format_identity(identity):
    return "NA" if identity is None else round(identity, 3)


def write_main_fasta_output(output_file, introns_data, duplicate_flag):
    """Writes the main FASTA-like output file."""
    with open(output_file, "w") as outfile:
        for intron in introns_data:
            info = intron["intron_info"]

            if duplicate_flag:
                if 'is_duplicate' in intron and intron['is_duplicate']:
                    intron_label = "_intron2"
                else:
                    intron_label = "_intron1"
            else:
                intron_label = ""

            outfile.write(f">Intron{intron_label},"
                          f"{info['reference_name']},"
                          f"RefPos:{info['intron_start_ref']}_{info['intron_end_ref']},"
                          f"{info['query_name']},"
                          f"Alignment:{info['orientation_alignment']},"
                          f"Intron_Length:{info['intron_length']},"
                          f"I:{info['insertions']},"
                          f"D:{info['deletions']},"
                          f"X:{info['mismatches']},"
                          f"QueryPos:{info['intron_start_query']},"
                          f"Identity:{_format_identity(info['identity'])},"
                          f"Splice_Site:{info['splice_site']},"
                          f"transcriptomic_support:{intron['transcriptomic_support']}\n"
                          )

            outfile.write(f"{info['intron_sequence']}\n"
                          f">Downstream_Flank{intron_label},"
                          f"{info['reference_name']},"
                          f"{info['down_flank_start']}_{info['down_flank_end']},"
                          f"{info['query_name']},{info['orientation_alignment']}\n"
                          f"{info['down_flank']}\n"
                          f">Upstream_Flank{intron_label},"
                          f"{info['reference_name']},"
                          f"{info['ups_flank_start']}_{info['ups_flank_end']},"
                          f"{info['query_name']},"
                          f"{info['orientation_alignment']}\n"
                          f"{info['ups_flank']}\n\n")

def write_iic_file(output_file, introns_data, duplicate_flag):
    """Writes the .iic output file."""
    with open(output_file, "w") as iic_file:
        for intron in introns_data:
            info = intron["intron_info"]
            
            if duplicate_flag:
                intron_label = "_intron1" if not ('is_duplicate' in intron and intron['is_duplicate']) else "_intron2"
            else:
                intron_label = ""

            iic_file.write(f"{info['reference_name']}_{info['query_name']},"
                           f"RefPos:{info['intron_start_ref']}_{info['intron_end_ref']},"
                           f"{intron_label}\t{info['ups_flank_200']}\t"
                           f"{info['intron_sequence']}\t"
                           f"{info['down_flank_200']}\n")

def write_correct_reads_list(output_file, correct_reads_list):
    """Writes the list of correct structure read IDs."""
    with open(output_file, "w") as correct_reads_file:
        correct_reads_file.write("\n".join(sorted(correct_reads_list)))
        if correct_reads_list:
            correct_reads_file.write("\n")

def _count(value):
    """Size of a set / value of a counter, or NA if that filtering step was not run."""
    if value is None:
        return "NA"
    return len(value) if isinstance(value, set) else value


def write_stats_summary(output_file, stats_collector, introns_data_final, alignment_filename,
                        correct_structure_reads):
    """
    Write the run statistics as a two-line TSV (header + values).

    reads      = query sequences (e.g. RNA-seq reads)
    references = reference sequences (contigs/chromosomes) with alignments
    clean_flanks = no insertion, deletion or mismatch in the -w flanks
    Columns of a filtering step that was not requested are NA.
    """
    st = stats_collector
    splice_sites = [extract_splice_site(i["intron_info"]["intron_sequence"]) for i in introns_data_final]
    ss = group_splice_site(splice_sites)

    columns = [
        ("alignment_file", alignment_filename),
        # reads
        ("reads_total", _count(st.total_metaT)),
        ("reads_spliced", _count(st.total_spliced_metaT)),
        ("reads_spliced_clean_flanks", _count(st.noindels_spliced_metaT)),
        ("reads_after_redundancy_removal", _count(st.noindels_noredundant_metaT)),
        ("reads_after_min_support", _count(st.noindels_noredundant_supported_metaT)),
        ("reads_all_introns_clean_flanks", len(correct_structure_reads)),
        # introns
        ("introns_total", _count(st.total_introns)),
        ("introns_clean_flanks", _count(st.noindels_introns)),
        ("introns_after_redundancy_removal", _count(st.noindels_nonredundant_introns)),
        ("introns_after_min_support", _count(st.noindels_nonredundant_supported_introns)),
        # reference sequences
        ("references_with_alignments", _count(st.total_metaG)),
        ("references_with_spliced_reads", _count(st.total_spliced_metaG)),
        ("references_clean_flanks", _count(st.noindels_spliced_metaG)),
        ("references_after_redundancy_removal", _count(st.noindels_noredundant_metaG)),
        ("references_after_min_support", _count(st.noindels_noredundant_supported_metaG)),
        # splice sites (each pooled with its reverse complement)
        ("gt_ag_count", ss["gt_ag_count"]),
        ("gc_ag_count", ss["gc_ag_count"]),
        ("at_ac_count", ss["at_ac_count"]),
        ("nci_count", ss["nci_count"]),
        ("gt_ag_pct", ss["gt_ag_pct"]),
        ("gc_ag_pct", ss["gc_ag_pct"]),
        ("at_ac_pct", ss["at_ac_pct"]),
        ("nci_pct", ss["nci_pct"]),
    ]
    for n in (1, 2, 3):
        columns += [
            (f"top{n}_ss", ss[f"top{n}_ss"].replace(" ", "/")),
            (f"top{n}_ss_count", ss[f"top{n}_ss_count"]),
            (f"top{n}_ss_pct", ss[f"top{n}_ss_pct"]),
        ]

    with open(output_file, "w") as stats_file:
        stats_file.write("\t".join(name for name, _ in columns) + "\n")
        stats_file.write("\t".join(str(value) for _, value in columns) + "\n")

def write_sl_output(input_bam, output_tsv, fasta_file=None, identity_tag="XI"):
    """
    Extracts 5' softclipped regions (potential SL sites) 
    and writes them to a TSV file.
    """
    with open(output_tsv, "w") as softclip_file:
        softclip_file.write("\t".join([
            "query_name", "query_length", "query_mapped_start", "query_mapped_end",
            "reference_name", "ref_start", "ref_end",
            "percent_identity", "sequence", "alignment_score", "strand"
        ]) + "\n")

        with open_alignment_file(input_bam, fasta_file) as bamfile:
            for read in bamfile.fetch(until_eof=True):
                if not keep_alignment(read):
                    continue
                if read.query_sequence is None:
                    continue
                if read.has_tag("SA"):  # exclude any read with SA:Z field
                    continue

                strand = "-" if read.is_reverse else "+"
                identity = get_identity(read, identity_tag)
                pid = "NA" if identity is None else identity
                score = get_alignment_score(read)
                score = "NA" if score is None else score

                trimmed_seq, _, _ = extract_unmapped_region(read)

                if trimmed_seq is None or len(trimmed_seq) < 15:  # 10 softclip + 5 mapped
                    continue

                query_mapped_start = read.query_alignment_start + 1
                query_mapped_end = read.query_alignment_end

                softclip_file.write("\t".join(map(str, [
                    read.query_name,
                    len(read.query_sequence),
                    query_mapped_start,  
                    query_mapped_end,    
                    read.reference_name,
                    read.reference_start + 1,
                    read.reference_end,
                    pid,
                    trimmed_seq,
                    score,
                    strand
                ])) + "\n")