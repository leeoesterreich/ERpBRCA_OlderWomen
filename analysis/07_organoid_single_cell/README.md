# Section 07: Organoid Single-Cell RNA-seq Analysis

## Analysis Overview

**Experiment:** PDO-296 breast cancer organoid treated with 7 conditions across 2 Flex pools.

| Sample | Treatment         | Pool |
|--------|-------------------|------|
| OS01   | Vehicle           | 1    |
| OS02   | E1                | 1    |
| OS03   | E1 + ICI          | 1    |
| OS04   | E1 + HSD17B7i     | 1    |
| OS05   | E2                | 2    |
| OS06   | E2 + ICI          | 2    |
| OS07   | E2 + HSD17B7i     | 2    |

### Pipeline Order

| Step | Script | Description |
|------|--------|-------------|
| Pre  | `build_reference.sh` | Build custom GRCh38 reference with HSD17B7 probe retention |
| 01   | `01_cellranger_multi_pool{1,2}.sh` | Cell Ranger multi (Flex) demultiplexing and counting |
| 02   | `02_qc.py` | Load per-sample Cell Ranger outputs, QC metrics, doublet detection (Scrublet), filtering |
| 03   | `03_preprocess.py` | Normalize, HVG selection, PCA, batch integration (Harmony), UMAP, Leiden clustering |
| 04   | `04_pseudobulk_de.py` | Pseudobulk differential expression (PyDESeq2) across treatment contrasts |
| 05   | `05_pathways.py` | GSEA (gseapy) and pathway activity inference (decoupleR/PROGENy) |
| 06   | `06_cell_cycle.py` | Cell cycle phase scoring and distribution across treatments |
| 07   | `07_single_cell_pathways.py` | Single-cell pathway scoring (AUCell-style) for hallmark and custom gene sets |
| 08   | `08_inhibitor_mechanism_exploration.py` | HSD17B7 inhibitor mechanism: target genes, dose-response signatures |
| 09   | `09_pathway_enrichment_analysis.py` | Pathway enrichment analysis on DE gene lists (Enrichr, MSigDB) |
| 10   | `10_heterogeneity_analysis.py` | Intra-treatment heterogeneity: entropy, dispersion, subcluster analysis |
| 11   | `11_proliferation_analysis.py` | Proliferation scoring (MKI67, TOP2A, cell cycle genes) across treatments |

### Running the Pipeline

```bash
# Submit as SLURM job (recommended)
sbatch run_analysis.sbatch

# Or run interactively (requires active conda env)
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging
bash run_analysis.sh
```

Cell Ranger (step 01) and reference building must be run separately before the main pipeline.

## Data File Locations & Update Checklist

| File | Current Location | Needs GEO Upload? | Notes |
|------|-----------------|-------------------|-------|
| Pool 1 FASTQs | `/ix1/alee/LO_LAB/General/Lab_Data/20260304_NeilOrganoidSingleCell_Alex/Lee_020326_AL1_FlexPool_1_OES9112A1` | **Yes** | Raw sequencing data |
| Pool 2 FASTQs | `/ix1/alee/LO_LAB/General/Lab_Data/20260304_NeilOrganoidSingleCell_Alex/Lee_020326_AL1_FlexPool_2_OES9112A2` | **Yes** | Raw sequencing data |
| Pool 1 Cell Ranger output | `/ix1/.../NeilOrganoidSingleCell/pool1_flex/outs/` | No | Derived from FASTQs |
| Pool 2 Cell Ranger output | `/ix1/.../NeilOrganoidSingleCell/pool2_flex/outs/` | No | Derived from FASTQs |
| Processed h5ad files | `/ix1/.../NeilOrganoidSingleCell/data/processed/` | No | Derived; reproducible from pipeline |
| Custom GRCh38 reference | `/ix1/.../NeilOrganoidSingleCell/references/refdata-gex-...` | No | Buildable via `build_reference.sh` |
| Probe set CSV | `/ix1/.../NeilOrganoidSingleCell/references/Chromium_...csv` | No | 10x Genomics public download |
| Sample metadata | `configs/sample_metadata.csv` (in repo) | N/A | Included in repo |

When paths change (e.g., after GEO upload), update `DATA_PATHS` in `_config.py`.

## GEO Submission Deliverables

| Deliverable | Format | Source |
|-------------|--------|--------|
| Raw count matrix (per sample) | MEX (matrix.mtx.gz + barcodes/features) | Cell Ranger `outs/per_sample_outs/` |
| Filtered count matrix (merged) | h5ad or CSV | `data/processed/qc_filtered.h5ad` |
| Cell metadata | TSV | `.obs` from final h5ad (sample, treatment, pool, leiden, cell_cycle_phase) |
| Sample metadata | CSV | `configs/sample_metadata.csv` |

**GEO series metadata:**

- Platform: 10x Chromium Flex
- Organism: Homo sapiens
- Cell line: PDO-296
- Reference genome: GRCh38-2024-A (custom, with HSD17B7 probe retention)
- Probe set: v1.1.0 with `filter-probes,false`
- Cell Ranger version: 9.0.1

## Probe Filtering Note

HSD17B7 is marked `included=FALSE` in the 10x probe set v1.1.0. We use `filter-probes,false` in Cell Ranger multi configs to retain all probes, ensuring HSD17B7 expression is captured. See `lessons_learned/10x_flex_probe_filtering.md` for details.

## Environment & Reproducibility

- **Environment:** `erp_brca_aging` at `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/erp_brca_aging`
- **Seeds:** Scripts use `seed=42` / `random_state=42` (deviation from project-wide convention of 12345, preserved for consistency with original analysis results)
- **Configuration:** All external paths centralized in `_config.py`

## Directory Layout

```
07_organoid_single_cell/
├── _config.py                          # Shared paths and constants
├── configs/
│   └── sample_metadata.csv             # Sample-to-treatment mapping
├── lessons_learned/
│   └── 10x_flex_probe_filtering.md     # Probe filtering documentation
├── 01_cellranger_multi_pool{1,2}.sh    # Cell Ranger multi configs
├── build_reference.sh                  # Custom reference builder
├── 02_qc.py ... 11_proliferation_analysis.py  # Analysis scripts
├── run_analysis.sh                     # Pipeline runner (interactive)
├── run_analysis.sbatch                 # Pipeline runner (SLURM)
├── submit_*.sh                         # Individual job submission scripts
├── outputs/                            # Analysis outputs (CSVs, h5ad)
├── figures/                            # Generated figures (PDF + PNG)
└── logs/                               # SLURM and runtime logs
```
