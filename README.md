# IntronXtract - Intron Extraction and Analysis

**IntronXtract** is a bioinformatics tool for extracting intron sequences and their exonic flanks from spliced RNA-seq alignments.
For every intron it reports the coordinates, the splice site, the identity of the read and the number of mismatches, insertions and deletions in the flanking exons, so that only introns with a reliable exon-intron boundary can be kept. 

## Contents

- [Installation](#installation)
- [How it works](#how-it-works)
- [Input requirements](#input-requirements)
- [Usage](#usage)
- [Outputs](#outputs)

## Installation

IntronXtract needs Python ≥ 3.8 and pysam.

**With conda/mamba**:

```bash
git clone https://github.com/lafucarmen/IntronXtract.git
cd IntronXtract
conda env create -f environment.yml
conda activate intronxtract
```
**With pip:**

```bash
pip install git+https://github.com/lafucarmen/IntronXtract.git
```

Check the installation:

```bash
intronxtract --version
```

## How it works

For every alignment, each `N` operation of the CIGAR is an intron. IntronXtract
takes the intron and `-w` nucleotides of exon on each side from the reference,
and counts the insertions (I), deletions (D) and mismatches (X) of the read inside those
two flanks:

```
reference      ...exon 1....|GT.............intron..............AG|....exon 2...
1-based pos.     41 ..... 50 51 ............................... 90 91 ..... 100
                [ upstream  ]                                      [ downstream ]
                  -w nt                                                -w nt
read           ...CGGCAACGACA------------- (N in CIGAR) -------------TTTCGCCGTT...
                    ^ I/D/X counted here                            ^ and here
```

An intron is *clean* when both flanks have no insertion, deletion or mismatch:
the read supports that exon–intron boundary base by base. With
`--filter-indel-free` only clean introns are reported.

## Input requirements

| Input | Notes |
|---|---|
| Alignments (`-i`) | SAM, BAM or CRAM from a **spliced aligner** (introns must appear as `N` in the CIGAR), e.g. `minimap2 -ax splice` (ONT) or `minimap2 -ax splice:hq` (PacBio). Sorting and indexing are not required but recommended. |
| Reference (`-f`) | The **same** FASTA used for the alignment, plain or bgzip-compressed. A `.fai` index is created automatically if missing. |

Example:

```bash
minimap2 -ax splice:hq -t 8 genome.fasta reads.fastq.gz | samtools view -b -o aln.bam
intronxtract -i aln.bam -f genome.fasta -o introns.fa -w 15 --filter-indel-free --stats
```

### Primary, secondary and supplementary alignments

Only unmapped records are ignored: primary, secondary and supplementary
alignments are all used, and the number of each is reported in the log.
Secondary alignments are alternative placements of the same read, so they can
add introns and increase the support of an intron (in one ONT test set they
added ~12 % introns). To keep only some of them, filter the file first:

```bash
samtools view -b -F 0x900 aln.bam > aln.primary.bam        # primary only
samtools view -b -F 0x100 aln.bam > aln.no_secondary.bam   # primary + supplementary
```

### SAM tags 

| Tag | Used for | If absent |
|---|---|---|
| `XI` or `--identity-tag`| Read identity | Computed from `NM` |
| `NM` | Read identity | Identity is `NA` |
| `MD` | Mismatches in the flanks | Read compared to the FASTA (needs the read sequence) |

**Identity.** When computed from `NM`, identity is the standard 
alignment identity: `100 × (1 − NM / alignment columns)`, where alignment
columns are the M/I/D/=/X bases. The `XI` tag can be added to the BAM file using
[bamFilters](https://github.com/institut-de-genomique/bamFilters)(option `-p`):

```bash
bamFilters -b aln.bam -o aln.xi.bam -p          # add the XI tag
bamFilters -b aln.bam -o aln.xi90.bam -p -i 90  # add XI and keep alignments with identity > 90 %
```

**Secondary alignments without sequence.** minimap2 writes secondary
alignments with `SEQ='*'` by default. IntronXtract still processes them if they
have an `MD` tag, and skips them with a warning otherwise. Both secondary sequences and 
MD tags should be requested when aligning: 

```bash
minimap2 -ax splice:hq --MD --secondary-seq -t 8 genome.fasta reads.fastq.gz \
    | samtools sort -o aln.bam
```


## Usage

```
intronxtract -i ALIGNMENTS -f FASTA -o OUTPUT [options]
```

### Examples

Extract all introns with 15 nt flanks:

```bash
intronxtract -i alignments.bam -f genome.fasta -o my_introns.fa -w 15
```

High-confidence introns: no indels/mismatches in the flanks, supported by ≥ 3
reads, reported once each, with statistics:

```bash
intronxtract -i alignments.bam -f genome.fasta -o my_introns.fa -w 15 \
    --filter-indel-free --min-support 3 --remove-redundancy --stats
```

Unknown read strand (e.g. unstranded cDNA): also write the reverse complement of
every intron, plus the `.iic` file:

```bash
intronxtract -i alignments.bam -f genome.fasta -o my_introns.fa --duplicate --intron-iic
```

Spliced-leader candidates (5′ soft clips ≥ 10 nt):

```bash
intronxtract -i alignments.bam -f genome.fasta -o my_introns.fa --sl-output
```

On a SLURM cluster, keep the log with the job:

```bash
intronxtract -i aln.bam -f genome.fasta -o out/introns.fa --stats 2> out/intronxtract.log
```

### Options

| Option | Description |
|---|---|
| `-i`, `--input` | Alignment file (SAM/BAM/CRAM). **Required** |
| `-f`, `--fasta` | Reference FASTA. **Required** |
| `-o`, `--output` | Main output file; other outputs share its prefix. **Required** |
| `-w`, `--window` | Flank length (nt) reported and checked for indels/mismatches (default 10) |
| `--identity-tag` | SAM tag with the read identity (default `XI`; falls back to `NM`) |
| `--filter-indel-free` | Keep only introns with no indels or mismatches in either flank |
| `--min-support N` | Keep only introns supported by ≥ N alignments |
| `--remove-redundancy` | Keep one record per intron: highest identity read (highest `AS` if there is no identity, otherwise a reproducible random choice) |
| `--duplicate` | Also write the reverse-complemented intron and flanks |
| `--intron-iic` | Write `PREFIX.iic` |
| `--stats` | Write `PREFIX_stats`, `PREFIX_correct_structure_reads` and `PREFIX_correct_intron_exon_structures.gff` |
| `--sl-output` | Write `PREFIX_sl_output.tsv` |
| `-q`, `--quiet` | Only print warnings and errors |
| `-V`, `--version` | Print version |



## Outputs

`PREFIX` is the output file name without its extension (`-o out/introns.fa`
→ `out/introns`). The examples below come from the [Quick start](#quick-start) run.

| File | Option | Content |
|---|---|---|
| `OUTPUT` | always | Intron and flank sequences, one block per intron and read |
| `PREFIX.iic` | `--intron-iic` | 200 nt exon + intron + 200 nt exon, one line per intron |
| `PREFIX_stats` | `--stats` | Run summary (one line) |
| `PREFIX_correct_structure_reads` | `--stats` | Reads whose introns all have clean flanks |
| `PREFIX_correct_intron_exon_structures.gff` | `--stats` | GFF3 of those reads |
| `PREFIX_sl_output.tsv` | `--sl-output` | 5′ soft-clipped sequences |

### Main output (`-o`)

Three FASTA records per intron and read, followed by a blank line:

```
>Intron,chr1,RefPos:51_90,r_perfect,Alignment:Forward,Intron_Length:40,I:0,D:0,X:0,QueryPos:30,Identity:99.5,Splice_Site:GT-AG,transcriptomic_support:8
GTAAGTTGTGCGGCGACCCTTGCGACAGTGACGCTTTCAG
>Downstream_Flank,chr1,91_100,r_perfect,Forward
TTTCGCCGTT
>Upstream_Flank,chr1,41_50,r_perfect,Forward
GGCAACGACA

```

| Field | Meaning |
|---|---|
| `chr1` | Reference sequence (contig) |
| `RefPos` | Intron start_end on the reference |
| `r_perfect` | Read name |
| `Alignment` | `Forward`/`Reverse`: strand of the alignment (not of the transcript) |
| `I`, `D`, `X` | Insertions, deletions, mismatches in both flanks (`-w` nt each side) |
| `QueryPos` | 0-based position, in the read sequence as stored in the BAM, of the first base after the splice junction |
| `Identity` | Read identity (see [SAM tags](#sam-tags-all-optional)); `NA` if unavailable |
| `Splice_Site` | First and last two intron nucleotides on the reference forward strand |
| `transcriptomic_support` | Number of alignments with the same intron (contig, start, end) |
| Flank headers | Reference, start_end, read, alignment strand |


With `--duplicate`, each intron also appears on the opposite strand, labelled
`_intron2` (the original is `_intron1`). Its intron is reverse-complemented and
the flanks swap sides, with their coordinates:

```
>Intron_intron2,chr1,RefPos:51_90,r_perfect,Alignment:Reverse,Intron_Length:40,I:0,D:0,X:0,QueryPos:30,Identity:99.5,Splice_Site:CT-AC,transcriptomic_support:8
CTGAAAGCGTCACTGTCGCAAGGGTCGCCGCACAACTTAC
>Downstream_Flank_intron2,chr1,41_50,r_perfect,Reverse
TGTCGTTGCC
>Upstream_Flank_intron2,chr1,91_100,r_perfect,Reverse
AACGGCGAAA
```

### `PREFIX.iic` (`--intron-iic`)

Input for [intronIC](https://github.com/glarue/intronIC) (classification of
U12-type / minor introns). Tab-separated, no header, one line per intron:
`name`, 200 nt upstream exon, intron, 200 nt downstream exon (shorter at contig
ends). This is the 4-column `.iic` format read by `intronIC -q`. Sequences are
shortened here:

```
chr1_r_perfect,RefPos:51_90,	AAGCCCAATAAACCACTCTG...	GTAAGTTGTGCGGCGACCCT...	TTTCGCCGTTGCCTAAACCT...
```

#### Classifying U2/U12-type introns with intronIC

```bash
intronxtract -i aln.bam -f genome.fasta -o introns.fa \
    --filter-indel-free --remove-redundancy --duplicate --intron-iic
intronIC -q introns.iic -n my_species
```

- `--filter-indel-free` keeps only introns whose junctions are supported base
  by base, so the splice-site motifs scored by intronIC are reliable.
- `--remove-redundancy` writes each intron once instead of once per read.
- `--duplicate` writes each intron on both strands (`_intron1` / `_intron2`),
  because unstranded reads do not tell which strand the transcript comes from.
  intronIC scores both; usually only one orientation shows canonical splice
  sites (e.g. `GT-AG` or `AT-AC` rather than `CT-AC` or `GT-AT`). With stranded
  data, `--duplicate` is not needed.

### `PREFIX_stats` (`--stats`)

Two-line tab-separated file (header + values). *Clean flanks* means no
insertion, deletion or mismatch in the `-w` nt on either side of the intron.
Columns of a filtering step that was not requested are `NA`. Transposed for
readability:

```
alignment_file                       reads.sam
reads_total                          10
reads_spliced                        8
reads_spliced_clean_flanks           5
reads_after_redundancy_removal       NA
reads_after_min_support              NA
reads_all_introns_clean_flanks       5
introns_total                        8
introns_clean_flanks                 5
introns_after_redundancy_removal     NA
introns_after_min_support            NA
references_with_alignments           1
references_with_spliced_reads        1
references_clean_flanks              1
references_after_redundancy_removal  NA
references_after_min_support         NA
gt_ag_count                          8
...
top1_ss                              CT-AC/GT-AG
top1_ss_count                        8
top1_ss_pct                          100.0
```

To view it this way: `paste <(head -1 PREFIX_stats | tr '\t' '\n') <(tail -1 PREFIX_stats | tr '\t' '\n')`

| Column | Meaning |
|---|---|
| `alignment_file` | Input file name |
| `reads_total` | Distinct read names in the file (mapped or not) |
| `reads_spliced` | Reads with at least one intron |
| `reads_spliced_clean_flanks` | Reads with at least one intron with clean flanks |
| `reads_after_redundancy_removal` | Reads still represented after `--remove-redundancy` |
| `reads_after_min_support` | Reads still represented after `--min-support` |
| `reads_all_introns_clean_flanks` | Reads whose introns *all* have clean flanks (listed in `PREFIX_correct_structure_reads`) |
| `introns_total` | Introns found (one per read and intron) |
| `introns_clean_flanks` | Introns with clean flanks |
| `introns_after_redundancy_removal` | Introns after `--remove-redundancy` |
| `introns_after_min_support` | Introns after `--min-support` |
| `references_with_alignments` | Reference sequences (contigs/chromosomes) with at least one alignment |
| `references_with_spliced_reads` | Reference sequences with at least one intron |
| `references_clean_flanks` | Reference sequences with at least one intron with clean flanks |
| `references_after_redundancy_removal` | Reference sequences after `--remove-redundancy` |
| `references_after_min_support` | Reference sequences after `--min-support` |
| `gt_ag_*`, `gc_ag_*`, `at_ac_*`, `nci_*` | Count and % of the final introns per splice-site class, each pooled with its reverse complement (`nci` = non-canonical; AT-AC introns are counted both in `at_ac_*` and `nci_*`) |
| `top{1,2,3}_ss*` | Three most frequent splice sites, their count and % |

### `PREFIX_correct_structure_reads` (`--stats`)

Names of the reads whose introns all have clean flanks, one per line, sorted:

```
r_1
r_2
r_3
r_4
r_5
```

### `PREFIX_correct_intron_exon_structures.gff` (`--stats`)

GFF3 with the exon/intron structure of those reads, ready for IGV, `gffread`
or AGAT. Each alignment is one `mRNA` (`Name` = read name); strand is `.`
because reads are not stranded:

```
##gff-version 3
chr1	intronxtract	mRNA	21	120	.	.	.	ID=aln1;Name=r_1
chr1	intronxtract	exon	21	50	.	.	.	ID=aln1.exon1;Parent=aln1
chr1	intronxtract	exon	91	120	.	.	.	ID=aln1.exon2;Parent=aln1
chr1	intronxtract	intron	51	90	.	.	.	ID=aln1.intron1;Parent=aln1
```

### `PREFIX_sl_output.tsv` (`--sl-output`)

Alignments without supplementary parts (`SA` tag) whose 5′ end (in read
orientation) is soft-clipped by ≥ 10 nt. `sequence` is the soft clip plus the
first 5 mapped nucleotides, e.g. to look for spliced-leader sequences:

```
query_name	query_length	query_mapped_start	query_mapped_end	reference_name	ref_start	ref_end	percent_identity	sequence	alignment_score	strand
r_sl	72	13	72	chr1	21	120	100.0	ACGTTTGACCAAACTGG	100	+
```
