# Rat Bulk RNA-seq Data

## Source

Raw FASTQ files are available from the corresponding author. Rat bulk RNA-seq data will be deposited in GEO as part of the study submission.

## Samples

6 rat mammary tumor samples, comparing Young vs Elderly:

| Sample | Group | Description |
|--------|-------|-------------|
| 102-FF-Tumor | CONTROL | Elderly rat tumor |
| 107-L4t2-Tumor | CONTROL | Elderly rat tumor |
| 116-R4-Tumor | CONTROL | Elderly rat tumor |
| 157-R2-Tumor | TEST | Young rat tumor |
| 158-L2-Tumor | TEST | Young rat tumor |
| 167-L2-Tumor | TEST | Young rat tumor |

## Reference Genome

- **Assembly:** mRatBN7.2
- **GTF:** `Rattus_norvegicus.mRatBN7.2.112.gtf`

A STAR index must be built locally from the mRatBN7.2 assembly before running the alignment pipeline.

## Files

- `sample_metadata.csv` — Sample annotations for DESeq2
