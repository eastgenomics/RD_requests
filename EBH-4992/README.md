# extract_split_reads.py

Extracts soft-clipped reads from a BAM file at a specified genomic region and writes them to FASTA for BLAST analysis.

Written to investigate a TSC1 exon 1 deletion detected by MLPA but not called by GATK gCNV, where a coverage artefact (double dip) and pile-up of soft-clipped reads was visible in IGV at chr9:132,944,000-132,946,000.

## Background

Soft-clipped reads occur when the end of a read cannot be aligned to the reference. A cluster of soft-clipped reads at a single genomic position is a strong indicator of a structural variant breakpoint. BLASTing the clipped sequences reveals where the other end of the breakpoint maps — identifying the SV type (deletion, inversion, mobile element insertion, etc.).

## Setup

```bash
python3 -m venv ~/venvs/splitreads
source ~/venvs/splitreads/bin/activate
pip install pysam
```

## Usage

The BAM must be indexed before running:

```bash
samtools index sample1.bam
```

Edit the variables at the top of the script. Note: chromosome name must match the contig names in the BAM header — this BAM uses `9` not `chr9`:

```python
BAM = "sample1.bam"          # path to indexed BAM
REGION = "9:132943000-132946000"
OUTPUT = "tsc1_clipped_reads.fa"
```

Then run:

```bash
python extract_split_reads.py
```

## Output

**To terminal:** a position summary table showing the number of soft-clipped reads at each position in the region, sorted by count. The position with the highest count is the likely breakpoint.

```
Position             Count   5prime   3prime
----------------------------------------------
132945412               87       12       75
132944891               14        9        5
...
```

**To file:** a FASTA file of unique clipped sequences from positions with >= `MIN_SUPPORTING_READS` reads.

## BLASTing the output

1. Go to https://blast.ncbi.nlm.nih.gov/Blast.cgi and click **Nucleotide BLAST**
2. Under "Enter Query Sequence", click **"Choose File"** and upload the output FASTA (e.g. `tsc1_clipped_reads.fa`)
3. Database: leave as **"Standard databases (nr etc.)"** with **"Nucleotide collection (nr/nt)"**
4. Under "Program Selection", select **"Somewhat similar sequences (blastn)"** — the clipped sequences may be short so this is more sensitive than megablast
5. Click **BLAST**

### Interpreting results

Each query sequence is a clipped read end from the suspected breakpoint. The top hits reveal where that sequence maps in the genome:

| Result | Interpretation |
|---|---|
| Hits to TSC1 region on chr9 | Inversion or complex rearrangement with both breakpoints in the gene |
| Hits to a repetitive element (LINE-1, Alu, SVA) | Mobile element insertion at the breakpoint |
| Hits to a different chromosome | Translocation |
| Hits to chr9 but far from TSC1 | Large deletion or inversion with a distal breakpoint |
| No hits | Sequence is in a poorly assembled region; try lowering `MIN_CLIP_LEN` to get longer clips |

The genomic coordinates of the top hit give you the second breakpoint of the structural variant.

## Tunable parameters

| Parameter | Default | Description |
|---|---|---|
| `MIN_CLIP_LEN` | 20 | Minimum soft-clip length in bp to consider |
| `MIN_MAPQ` | 20 | Minimum read mapping quality |
| `MIN_SUPPORTING_READS` | 3 | Minimum reads at a position to include in FASTA |
