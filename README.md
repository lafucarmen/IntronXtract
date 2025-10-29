# intronxtract - Intron Extraction and Analysis

`intronxtract` is a bioinformatics tool for extracting intron sequences and their flanking regions from `BAM` alignment files.

## Usage

```
intronxtract [-h] -i INPUT -f FASTA -o OUTPUT -w WINDOW 
             [--sl-output] [--filter_indel_free] 
             [--transcriptomic_support MIN_READS] [--remove_redundancy] 
             [--duplicate] [--stats] [--intron_iic]
```

## Examples

The script requires, at minimum, a BAM file, the reference genome FASTA, and an output file.

### Basic usage
Extract introns and flanking regions of 15nt length
```bash
intronxtract \
    -i alignments.bam \
    -f genome.fasta \
    -o my_introns \
    -w 15
```
### Advanced usage 
Extract high-confidence introns by applying different filters and generates a stats file.

```bash
intronxtract \
    -i alignments.bam \
    -f genome.fasta \
    -o my_introns.fasta \
    -w 15 \
    --filter_indel_free \
    --transcriptomic_support 3 \
    --remove_redundancy
    --stats 
