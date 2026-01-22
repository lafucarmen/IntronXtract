import pysam
import os
from .stats import group_splice_site
from .sequence_utils import extract_splice_site
from .bam_processing import extract_unmapped_region

def get_intron_header_attributes(intron_info, transcriptomic_support, intron_label=""):
    """
    Constructs the GFF header attributes for an intron.
    """
    identity_value = intron_info.get('identity', "NA")
    if isinstance(identity_value, (int, float)):
        identity_str = str(round(identity_value, 3))
    else:
        identity_str = identity_value
        
    as_score = intron_info.get('alignment_score', "NA")

    header_attrs = [
        f"ID=Intron{intron_label}.{intron_info['query_name']}.{intron_info['intron_start_ref']}_{intron_info['intron_end_ref']}",
        f"Parent={intron_info['query_name']}",
        f"RefPos={intron_info['intron_start_ref']}_{intron_info['intron_end_ref']}",
        f"Alignment={intron_info['orientation_alignment']}",
        f"Intron_Length={intron_info['intron_length']}",
        f"I={intron_info.get('insertions', 'NA')}",
        f"D={intron_info.get('deletions', 'NA')}",
        f"X={intron_info.get('mismatches', 'NA')}",
        f"QueryPos={intron_info.get('intron_start_query', 'NA')}",
        f"Identity={identity_str}", 
        f"AS={as_score}", 
        f"Splice_Site={intron_info.get('splice_site', 'NA_NA').replace('-', '_')}",
        f"Support={transcriptomic_support}"
    ]
    return ";".join(header_attrs)

def write_gff(bam_file, output_file, introns_data):
    introns_lookup = {}
    for intron in introns_data:
        info = intron['intron_info']
        key = (info['query_name'], info['intron_start_ref'])
        introns_lookup[key] = intron

    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_file, "w") as gff_file:
        gff_file.write("##gff-version 3\n")

        for read in bam.fetch():
            if read.is_unmapped or read.seq is None:
                continue
            if 3 not in [op for op, _ in read.cigartuples]:
                continue
            
            if read.is_supplementary:
                aln_type = "supp"
            elif read.is_secondary:
                aln_type = "secondary"
            else:
                aln_type = "primary"

            ref_name = read.reference_name
            strand = "-" if read.is_reverse else "+"
            parent_id = read.query_name
            
            unique_parent_id = f"{parent_id}_{aln_type}"
            
            gff_file.write(
                f"{ref_name}\tIntronXtract\tmRNA\t{read.reference_start + 1}\t{read.reference_end}\t.\t{strand}\t.\tID={unique_parent_id};Name={parent_id};AlnType={aln_type}\n"
            )

            current_pos = read.reference_start + 1
            exon_start = current_pos
            exon_positions = []
            intron_positions = []

            for op, length in read.cigartuples:
                if op in [0, 2, 7, 8]: 
                    current_pos += length
                elif op == 3: # Intrón
                    exon_positions.append((exon_start, current_pos - 1))
                    intron_positions.append((current_pos, current_pos + length - 1))
                    current_pos += length
                    exon_start = current_pos
            exon_positions.append((exon_start, read.reference_end))

            for idx, (i_start, i_end) in enumerate(intron_positions, 1):
                key = (parent_id, i_start)
                target = introns_lookup.get(key)

                if target:
                    attributes = get_intron_header_attributes(
                        target['intron_info'],
                        target['transcriptomic_support'],
                        intron_label=f"_intron{idx}"
                    )
                else:
                    as_tag = read.get_tag("AS") if read.has_tag("AS") else "NA"
                    attributes = f"ID={unique_parent_id}.int{idx};Parent={unique_parent_id};Intron_Length={i_end - i_start + 1};AS={as_tag};AlnType={aln_type}"

                gff_file.write(f"{ref_name}\tIntronXtract\tintron\t{i_start}\t{i_end}\t.\t{strand}\t.\t{attributes}\n")

            for e_start, e_end in exon_positions:
                gff_file.write(f"{ref_name}\tIntronXtract\texon\t{e_start}\t{e_end}\t.\t{strand}\t.\tParent={unique_parent_id}\n")

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
                          f"Identity:{round(info['identity'], 3)},"
                          f"AS:{info['alignment_score']},"
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
        correct_reads_file.write("\n".join(correct_reads_list))

def write_stats_summary(output_file, stats_collector, introns_data_final, bam_filename, correct_structure_metaT):
    """Calculates and writes the final statistics summary file."""
    with open(output_file, "w") as stats_file:
        stats_file.write("bamfile\t"
                         "total_metaT\t"
                         "total_spliced_metaT\t"
                         "spliced_metaT_no_indels\t"
                         "spliced_metaT_no_indels_non_redundant\t"
                         "spliced_metaT_no_indels_non_redundant_supported\t"
                         "total_introns\t"
                         "introns_no_indels\t"
                         "introns_no_indels_non_redundant\t"
                         "introns_no_indels_non_redundant_supported\t"
                         "metaT_correct_structure\t"
                         "total_metaG\ttotal_spliced_metaG\t"
                         "spliced_metaG_no_indels\t"
                         "spliced_metaG_no_indels_non_redundant\t"
                         "spliced_metaG_no_indels_non_redundant_supported\t"
                         "gt_ag_count\t"
                         "gc_ag_count\t"
                         "at_ac_count\t"
                         "nci_count\t"
                         "gt_ag_pct\t"
                         "gc_ag_pct\t"
                         "at_ac_pct\t"
                         "nci_pct\t"
                         "top1_ss\t"
                         "top1_ss_count\t"
                         "top1_ss_pct\t"
                         "top2_ss\t"
                         "top2_ss_count\t"
                         "top2_ss_pct\t"
                         "top3_ss\t"
                         "top3_ss_count\t"
                         "top3_ss_pct\n")

        intron_sequences = [intron["intron_info"]["intron_sequence"] for intron in introns_data_final]
        splice_sites = [extract_splice_site(seq) for seq in intron_sequences]
        splice_stats = group_splice_site(splice_sites)

        stats_file.write(f"{bam_filename}\t{len(stats_collector.total_metaT)}\t{len(stats_collector.total_spliced_metaT)}\t"
                         f"{len(stats_collector.noindels_spliced_metaT)}\t{len(stats_collector.noindels_noredundant_metaT)}\t"
                         f"{len(stats_collector.noindels_noredundant_supported_metaT)}\t{stats_collector.total_introns}\t"
                         f"{stats_collector.noindels_introns}\t{stats_collector.noindels_nonredundant_introns}\t"
                         f"{stats_collector.noindels_nonredundant_supported_introns}\t{len(correct_structure_metaT)}\t"
                         f"{len(stats_collector.total_metaG)}\t{len(stats_collector.total_spliced_metaG)}\t{len(stats_collector.noindels_spliced_metaG)}\t"
                         f"{len(stats_collector.noindels_noredundant_metaG)}\t{len(stats_collector.noindels_noredundant_supported_metaG)}\t"
                         f"{splice_stats['gt_ag_count']}\t"
                         f"{splice_stats['gc_ag_count']}\t"
                         f"{splice_stats['at_ac_count']}\t"
                         f"{splice_stats['nci_count']}\t"
                         f"{splice_stats['gt_ag_pct']}\t"
                         f"{splice_stats['gc_ag_pct']}\t"
                         f"{splice_stats['at_ac_pct']}\t"
                         f"{splice_stats['nci_pct']}\t"
                         f"{splice_stats['top1_ss'].replace(' ', '/')}\t"
                         f"{splice_stats['top1_ss_count']}\t"
                         f"{splice_stats['top1_ss_pct']}\t"
                         f"{splice_stats['top2_ss'].replace(' ', '/')}\t"
                         f"{splice_stats['top2_ss_count']}\t"
                         f"{splice_stats['top2_ss_pct']}\t"
                         f"{splice_stats['top3_ss'].replace(' ', '/')}\t"
                         f"{splice_stats['top3_ss_count']}\t"
                         f"{splice_stats['top3_ss_pct']}\n")


def write_sl_output(input_bam, output_tsv):
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

        with pysam.AlignmentFile(input_bam, "rb") as bamfile:
            for read in bamfile:
                if read.is_unmapped or read.query_sequence is None:
                    continue
                if read.has_tag("SA"):  # exclude any read with SA:Z field
                    continue

                strand = "-" if read.is_reverse else "+"
                pid = "NA"
                score = "NA"

                for tag in read.get_tags():
                    if tag[0] == "XI":
                        pid = tag[1]
                    elif tag[0] == "AS":
                        score = tag[1]

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