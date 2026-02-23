# Rat Bulk RNA-seq Data

## Source
- **Raw FASTQs:** `/ix1/alee/LO_LAB/General/Lab_Data/20240622_Neil_Rahul_RNASeqRat/`
- **Original analysis:** `/ix1/alee/LO_LAB/Personal/Rahul/Neil_RNAseq/`

## Samples
6 rat mammary tumor samples, comparing Young (TEST) vs Elderly (CONTROL):

| Sample | Group | Description |
|--------|-------|-------------|
| 102-FF-Tumor | TEST | Young rat tumor |
| 107-L4t2-Tumor | TEST | Young rat tumor |
| 116-R4-Tumor | TEST | Young rat tumor |
| 157-R2-Tumor | CONTROL | Elderly rat tumor |
| 158-L2-Tumor | CONTROL | Elderly rat tumor |
| 167-L2-Tumor | CONTROL | Elderly rat tumor |

## Reference Genome
- **Assembly:** mRatBN7.2
- **STAR index:** `/ix1/alee/LO_LAB/Personal/Rahul/Reference_genome/Rat/`
- **GTF:** `Rattus_norvegicus.mRatBN7.2.112.gtf`

## Files
- `sample_metadata.csv` - Sample annotations for DESeq2
