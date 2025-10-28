from Bio.Seq import Seq

class StatsCollector:
    """A class to collect statistics during the run."""
    def __init__(self):
        self.total_introns = 0
        self.noindels_introns = 0
        self.indels_introns = 0
        self.noindels_nonredundant_introns = 0
        self.noindels_nonredundant_supported_introns = 0
        
        self.total_metaT = set()
        self.total_spliced_metaT = set()
        self.indels_spliced_metaT = set()
        self.noindels_spliced_metaT = set()
        self.noindels_noredundant_metaT = set()
        self.noindels_noredundant_supported_metaT = set()

        self.total_metaG = set()
        self.total_spliced_metaG = set()
        self.noindels_spliced_metaG = set()
        self.noindels_noredundant_metaG = set()
        self.noindels_noredundant_supported_metaG = set()
    
    def add_spliced_read(self, read, has_indel=False):
        """Updates counters for a processed intron."""
        self.total_introns += 1
        self.total_spliced_metaT.add(read.query_name)
        self.total_spliced_metaG.add(read.reference_name)

        if has_indel:
            self.indels_introns += 1
            self.indels_spliced_metaT.add(read.query_name)
        else:
            self.noindels_introns += 1
            self.noindels_spliced_metaT.add(read.query_name)
            self.noindels_spliced_metaG.add(read.reference_name)

def group_splice_site(splice_sites):
    group_counts = {}
    total_count = len(splice_sites)
    for sequence in splice_sites:
        seq = Seq(sequence)
        reverse_comp = str(seq.reverse_complement())
        sorted_group = " ".join(sorted([sequence, reverse_comp]))
        group_counts[sorted_group] = group_counts.get(sorted_group, 0) + 1

    # Initialize counts
    gt_ag_count = 0
    gc_ag_count = 0
    at_ac_count = 0
    nci_count = 0

    # Categorize splice site types
    for group, count in group_counts.items():
        if "GT-AG" in group or "CT-AC" in group: 
            gt_ag_count += count
        elif "GC-AG" in group or "CT-GC" in group: 
            gc_ag_count += count
        elif "AT-AC" in group or "GT-AT" in group: 
            at_ac_count += count
        else:
            nci_count += count
            
    # Re-logic to match original script (AT-AC was also counted as NCI)
    nci_count_corrected = 0
    for group, count in group_counts.items():
        is_canonical = ("GT-AG" in group or "CT-AC" in group) or \
                       ("GC-AG" in group or "CT-GC" in group)
        is_atac = "AT-AC" in group or "GT-AT" in group
        
        if not is_canonical and not is_atac:
             nci_count_corrected += count
    
    nci_count = nci_count_corrected + at_ac_count


    # Compute percentages
    gt_ag_percentage = round((gt_ag_count / total_count) * 100, 2) if total_count else 0
    gc_ag_percentage = round((gc_ag_count / total_count) * 100, 2) if total_count else 0
    at_ac_percentage = round((at_ac_count / total_count) * 100, 2) if total_count else 0
    nci_percentage = round((nci_count / total_count) * 100, 2) if total_count else 0

    # Sort splice site groups by frequency
    sorted_groups = sorted(group_counts.items(), key=lambda x: x[1], reverse=True)

    # Extract top 3 most common splice sites
    top_splice_sites = [(sorted_groups[i] if i < len(sorted_groups) else ("N/A", 0)) for i in range(3)]

    return {
        "gt_ag_count": gt_ag_count,
        "gc_ag_count": gc_ag_count,
        "at_ac_count": at_ac_count,
        "nci_count": nci_count,
        "gt_ag_pct": gt_ag_percentage,
        "gc_ag_pct": gc_ag_percentage,
        "at_ac_pct": at_ac_percentage,
        "nci_pct": nci_percentage,
        "top1_ss": top_splice_sites[0][0],
        "top1_ss_count": top_splice_sites[0][1],
        "top1_ss_pct": round((top_splice_sites[0][1] / total_count) * 100, 2) if total_count else 0,
        "top2_ss": top_splice_sites[1][0],
        "top2_ss_count": top_splice_sites[1][1],
        "top2_ss_pct": round((top_splice_sites[1][1] / total_count) * 100, 2) if total_count else 0,
        "top3_ss": top_splice_sites[2][0],
        "top3_ss_count": top_splice_sites[2][1],
        "top3_ss_pct": round((top_splice_sites[2][1] / total_count) * 100, 2) if total_count else 0
    }
