# Section 07: Organoid Single-Cell RNA-seq Analysis

## Analysis Overview

**Experiment:** PDO-296 breast cancer organoid treated with 7 conditions across 2 Flex pools.

| Sample | Treatment         | Pool |
|--------|-------------------|------|
| OS01   | Vehicle           | 1    |
| OS02   | E1                | 1    |
| OS03   | E1 + fulvestrant  | 1    |
| OS04   | E1 + HSD17B7i     | 1    |
| OS05   | E2                | 2    |
| OS06   | E2 + fulvestrant  | 2    |
| OS07   | E2 + HSD17B7i     | 2    |

## Pipeline Order

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

## Running the Pipeline

```bash
# Submit as SLURM job (recommended)
sbatch run_analysis.sbatch

# Or run interactively
bash run_analysis.sh
```

Cell Ranger (step 01) and reference building must be run separately before the main pipeline. See `environment.yml` in the repository root for conda environment setup.

## Probe Filtering Note

HSD17B7 is marked `included=FALSE` in the 10x probe set v1.1.0. The Cell Ranger multi configurations use `filter-probes,false` to retain all probes, ensuring HSD17B7 expression is captured. See `lessons_learned/10x_flex_probe_filtering.md` for details.

## Reproducibility

- **Seeds:** Scripts use `seed=42` / `random_state=42`
- **Configuration:** Shared paths and constants are centralized in `_config.py`

## Directory Layout

```
07_organoid_single_cell/
├── _config.py                        # Shared paths and constants
├── _figure_config.py                 # Figure styling configuration
├── configs/
│   ├── cellranger_multi_config_pool1.csv
│   ├── cellranger_multi_config_pool2.csv
│   └── sample_metadata.csv           # Sample-to-treatment mapping
├── lessons_learned/
│   └── 10x_flex_probe_filtering.md
├── build_reference.sh                # Custom reference builder
├── 01_cellranger_multi_pool{1,2}.sh  # Cell Ranger multi configs
├── 02_qc.py ... 11_proliferation_analysis.py  # Analysis scripts
├── run_analysis.sh                   # Pipeline runner (interactive)
├── run_analysis.sbatch               # Pipeline runner (SLURM)
└── submit_*.sh                       # Individual step submission scripts
```
