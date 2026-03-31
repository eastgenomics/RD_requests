# TSC1 Exon 1 Structural Variant Analysis

**Sample:** 141351816-26013R0074-26NGCEN7-9527-M-99347387
**Gene:** TSC1 (chr9, minus strand, GRCh38)
**Date:** 2026-03-31

---

## Background

Routine MLPA and panel NGS (targeted sequencing, ~3000x depth) were performed for TSC1. A discrepancy was identified between the two methods.

---

## Initial Findings

### MLPA
Heterozygous deletion of **exon 1** detected (probe: chr9:132,944,553–132,944,592).

### Panel NGS — GATK gCNV
**No deletion called.**

### IGV
A **double-dip coverage artefact** was visible in Sample 1 at approximately chr9:132,944,000–132,946,000, not observed in a comparison sample. A dense pile-up of soft-clipped reads was visible at the right edge of the coverage drop.

---

## Hypothesis

The MLPA–gCNV discrepancy and the IGV artefact are consistent with a structural variant (SV) at TSC1 exon 1. A soft-clipped read analysis was performed to identify the SV breakpoint(s).

---

## Methods

### Soft-clipped read extraction

A custom Python script (`extract_split_reads.py`) was written to extract soft-clipped reads from the BAM file using `pysam`. Soft-clipped reads (minimum clip length 20bp, minimum MAPQ 20) were extracted across the region of interest. Clipped sequences were written to FASTA and submitted to NCBI nucleotide BLAST (blastn, somewhat similar sequences, core_nt database).

Two search regions were used:
1. `9:132,943,000–132,946,000` (initial, narrow)
2. `9:132,935,000–132,946,000` (extended, after considering MLPA probe positions)

### MLPA probe positions

The 23 MLPA probes span chr9:132,896,694–132,944,592. Key probes:

| Probe | Position | MLPA result |
|---|---|---|
| Probe 22 | chr9:132,935,017–132,935,056 | Normal (2 copies) |
| Probe 23 (exon 1) | chr9:132,944,553–132,944,592 | **Deletion (1 copy)** |

---

## Results

### Soft-clip position summary

| Position | Side | Supporting reads | Interpretation |
|---|---|---|---|
| chr9:132,945,047 | 3' clip | **181** | Confirmed SV breakpoint |
| chr9:132,944,965 | 5' clip | 120 | Adjacent to upstream breakpoint |
| chr9:132,944,816 | 5' clip | 60 | Near panel coverage edge |
| chr9:132,944,295 | 5' clip | 7 | Panel edge artefact (see below) |
| chr9:132,944,138 | 3' clip | 11 | Panel edge artefact (see below) |

### Confirmed breakpoint: chr9:132,945,047

The dominant signal of 181 soft-clipped reads at chr9:132,945,047 is located 1,129bp above the panel coverage gap edge and is considered a genuine SV breakpoint.

BLAST analysis of the clipped sequences identified two distinct populations of hits against the TSC1 RefSeqGene (NG_012386.1, positions ~4,468–4,585):

| Clip length | BLAST strand | Identity | Interpretation |
|---|---|---|---|
| Short (25–66bp) | Plus/Plus | 100% | Normal TSC1 sequence clipping at the breakpoint edge |
| Long (56–121bp) | Plus/Minus | 73–87% | Clipped sequence crosses the breakpoint junction; maps to TSC1 in reverse complement with imperfect identity |

The longer clips with imperfect identity mapping to TSC1 on the reverse strand are **not consistent with a simple clean deletion**, in which junction-spanning clips would be expected to map to the other side of the deletion with ~100% identity. This pattern is more consistent with a complex rearrangement at the breakpoint junction (e.g., inversion, deletion-inversion, or microhomology-mediated repair).

### Panel coverage gap

A critical constraint was identified: **there are no reads between chr9:132,929,294 and chr9:132,943,918** (a 14,624bp gap in panel coverage).

```
chr9 position:

132,929,294                132,943,918  132,944,553  132,945,047
      |                          |            |            |
      |<====== NO READS =========|            |            |
                                 |<-- 635bp ->|<-- 494bp ->|
                          gap edge        exon 1 probe   breakpoint
```

Positions of the weak downstream signals relative to the gap:

| Signal | Position | Distance from gap edge |
|---|---|---|
| 11-read signal | 132,944,138 | 220bp above gap |
| 7-read signal | 132,944,295 | 377bp above gap |

Both signals fall within 377bp of the coverage gap edge. At 3,000x depth, soft-clipping at panel capture boundaries is a recognised artefact. **These signals are not considered genuine breakpoints.**

### Downstream breakpoint location

The downstream breakpoint is constrained by the MLPA result:
- Must be **below** 132,944,553 (to include the exon 1 probe in the SV)
- Must be **above** 132,935,056 (probe 22 is not affected)
- This 9,497bp window is **93% within the panel coverage gap**

The downstream breakpoint almost certainly falls within the coverage gap and **cannot be detected from this panel NGS dataset.**

### Why gCNV failed to call the deletion

GATK gCNV requires read depth signal across the affected region. As the SV extends into a panel coverage gap, there are no reads in the deleted/rearranged region to register a depth reduction. The exon 1 probe region (132,944,553–132,944,592), while covered, is situated within a CpG island (GC-rich TSC1 promoter), which causes high baseline coverage variability. This combination — partial coverage gap overlap and GC-driven noise at the only covered probe — is sufficient to explain the gCNV failure.

---

## Conclusions

### What is known with confidence

1. There is a structural variant affecting TSC1 with **one confirmed breakpoint at chr9:132,945,047** (181 supporting soft-clipped reads).
2. The SV includes the MLPA exon 1 probe region (chr9:132,944,553–132,944,592), consistent with the MLPA deletion signal.
3. The downstream breakpoint is within the panel coverage gap (chr9:132,929,294–132,943,918) and is undetectable from this data.
4. The BLAST signature at the upstream breakpoint (imperfect identity, reverse-strand mapping of longer clips) suggests the junction is **not a simple clean deletion** and may involve a complex rearrangement such as an inversion or deletion-inversion.

### What remains unresolved

| Question | Status |
|---|---|
| Exact downstream breakpoint coordinate | Unknown — in coverage gap |
| SV size | Unknown — between ~1.5kb and ~10kb |
| SV type (deletion, inversion, complex) | Unknown — cannot distinguish from panel NGS |

---

## Recommended Next Steps

### Definitive
**Oxford Nanopore long-read sequencing** is the most appropriate follow-up. Long reads (typically 10–50kb) will span both the upstream breakpoint and the coverage gap in a single read, directly revealing:
- The exact downstream breakpoint coordinate
- The SV type (deletion, inversion, complex rearrangement)
- The full junction sequence

### Alternative
**Sanger sequencing with junction-spanning primers** — if the downstream breakpoint can be narrowed further (e.g., by array CGH), primers can be designed flanking both breakpoints. A deletion would produce a short band; an inversion would only amplify from one orientation.

### Not recommended
Further analysis of this panel NGS BAM is unlikely to yield additional information. The limiting factor is the coverage gap, which cannot be overcome by lowering analysis thresholds or widening the search region.

---

## Files

| File | Description |
|---|---|
| `extract_split_reads.py` | Script used for soft-clip extraction |
| `tsc1_clipped_reads.fa` | FASTA of clipped sequences (initial region) |
| `WR14D0TP014-Alignment.txt` | NCBI BLAST results (initial region, 835 queries) |
| `WR2KK08V016-Alignment-wide.txt` | NCBI BLAST results (extended region, 835 queries) |
