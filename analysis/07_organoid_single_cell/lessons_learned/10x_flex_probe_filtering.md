# 10X Flex Probe Set: Gene Exclusion via `included=FALSE`

## Problem

HSD17B7 — the primary gene of interest for this project — was completely absent from Cell Ranger filtered output despite being a well-annotated protein-coding gene in GRCh38.

## Root Cause

In the **Chromium Human Transcriptome Probe Set v1.1.0** (GRCh38-2024-A), HSD17B7 has only **1 probe** and it is marked `included=FALSE`:

```
ENSG00000132196,GAAGACC...,ENSG00000132196|HSD17B7|d3308cd,FALSE,unspliced,HSD17B7
```

The `included` column in the probe set CSV controls whether Cell Ranger uses a probe's UMI counts in the **filtered** gene expression matrix. Probes marked `FALSE` are excluded because 10X predicts they have off-target activity to homologous genes/sequences.

For comparison, all other HSD17B family members (HSD17B1-6, HSD17B8-14) have 2-3 probes each, all marked `TRUE`.

## Key Insight

- The **physical probe IS in the kit** and was hybridized — this is not a wet-lab issue
- HSD17B7 **does appear in the raw matrix** (`raw_feature_bc_matrix`) for all samples
- It is only excluded from the **filtered matrix** (`filtered_feature_bc_matrix`) by Cell Ranger's computational filter
- The GEM-X Flex v1 kit with Human Transcriptome Probe Kit (PN-1000785) correctly maps to probe set v1.1.0

## Fix

Add `filter-probes,false` to the `[gene-expression]` section of the Cell Ranger multi config CSV:

```csv
[gene-expression]
reference,/path/to/reference
probe-set,/path/to/probe_set.csv
filter-probes,false
create-bam,false
```

This tells Cell Ranger to include UMI counts from all non-deprecated probes, including those with predicted off-target activity.

**Trade-off**: Genes with `included=FALSE` probes may have some off-target signal. Interpret with appropriate caution, but for a well-expressed target gene like HSD17B7 in a relevant tissue, the signal is likely genuine.

## Lesson

**Always check if your key genes of interest are in the probe set AND marked `included=TRUE`** before running Cell Ranger. A gene can be "in" the probe set but still excluded from output.

To check:
```bash
grep "YOUR_GENE" Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv
```

Look at the 4th column — if it says `FALSE`, you need `filter-probes=false` in your config.

## Date

2026-03-09
