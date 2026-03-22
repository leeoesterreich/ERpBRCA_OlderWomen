# ERpBRCA_OlderWomen CodeOcean Capsule Design

**Date:** 2026-03-22
**Target journal:** Nature Aging
**Reference capsule:** CITEgeistNeilAnalysis/CITEgeist_CodeOcean (1.7 GB total, 1.4 GB data, Python-only, 6 steps)

## Goal

Create a CodeOcean capsule for the ERpBRCA_OlderWomen publication that provides:
1. **Peer reviewer verification** — reviewers regenerate figures at a glance
2. **Full computational reproducibility** — pipeline recomputes from count matrices
3. **Code archive** — citable, inspectable, DOI-linked

Nature Aging has a CodeOcean partnership; capsules are provisioned during peer review and published with a DOI.

## Constraints

| Resource | Limit | Notes |
|----------|-------|-------|
| Git repo (code/) | < 1 GB | Scripts only, no data |
| Individual git file | < 100 MB | |
| Data assets | Tier-dependent | See Data Strategy |
| Free academic storage | 20 GB total | Shared across all capsules; ~4 GB available |
| Free academic compute | 10 hr/month | May be expanded by Nature partnership |

## Architecture: Section-Isolated Pipeline (Approach B)

Each of the 7 analysis sections is a self-contained module with its own `run.sh`. A master orchestrator runs them sequentially, but sections are independently runnable. Precomputed fallbacks enable fast figure regeneration.

### Capsule Directory Structure

```
ERpBRCA_CodeOcean/
├── README.md
├── environment/
│   ├── Dockerfile                     # bioconductor/bioconductor_docker:3.19 + miniconda
│   └── postInstall                    # R/Python package installation
├── code/
│   ├── run_all.sh                     # Full pipeline from count matrices
│   ├── run_codeocean.sh               # Figure regen from precomputed intermediates
│   ├── run_section.sh                 # Single section: ./run_section.sh 02_rat_snrnaseq
│   ├── config/
│   │   ├── paths.R                    # Centralized path config with fallback logic
│   │   ├── paths.py                   # Python equivalent
│   │   ├── figure_config.R            # Shared plot styling
│   │   └── figure_config.py           # Shared plot styling
│   ├── 01_human_bulk_rnaseq/
│   │   ├── run.sh
│   │   ├── mica/                      # MICA subpipeline (7 scripts)
│   │   └── *.R                        # 5 main scripts
│   ├── 02_rat_snrnaseq/
│   │   ├── run.sh
│   │   └── *.R                        # 9 scripts
│   ├── 03_comparison/
│   │   ├── run.sh
│   │   └── *.R                        # 1 script
│   ├── 03_rat_wes/
│   │   ├── run.sh
│   │   └── *.py                       # 3 scripts (post-VEP)
│   ├── 04_human_scrnaseq/
│   │   ├── run.sh
│   │   └── *.R, *.py                  # 21 scripts (Xu et al. 2024 atlas)
│   ├── 05_rat_bulk_rnaseq/
│   │   ├── run.sh
│   │   └── *.R                        # 3 scripts (DESeq2 onward)
│   ├── 06_spatial_biopsies/
│   │   ├── run.sh
│   │   └── *.py                       # 2 analysis + 2 utility modules
│   └── 07_organoid_single_cell/
│       ├── run.sh
│       └── *.py                       # 10 analysis + 2 config modules
├── data/                              # Mounted read-only data asset (gitignored)
│   ├── README.md
│   ├── human_bulk_rnaseq/
│   ├── rat_snrnaseq/
│   ├── rat_wes/
│   ├── human_scrnaseq/
│   ├── rat_bulk_rnaseq/
│   ├── spatial_biopsies/
│   ├── organoid_single_cell/
│   ├── gmt/
│   ├── precomputed/                   # Intermediate RDS/h5ad (Tier 2 only)
│   │   └── {01..07}_outputs/
│   └── figure_regen/                  # Lightweight CSVs for figure regeneration
│       └── {01..07}_outputs/
└── results/
    ├── figures/{01..07}_*/
    ├── tables/
    ├── intermediate/{01..07}_outputs/
    └── logs/
```

## Per-Section Script Inclusion Lists

Defines exactly which scripts are included in the capsule vs excluded.

### Section 01: Human Bulk RNA-seq (R)
**Included:**
- `01_preprocess.R`, `02_run_gsva.R`, `03_run_progeny.R`, `04_correlations.R`, `05_visualize.R`
- `mica/01_prep_data.R`, `mica/02_run_gsva.R`, `mica/03_run_mica.R`, `mica/04_fig4_estrogen.R`, `mica/05_fig5_immune.R`, `mica/07_statistical_audit.R`

**Excluded:**
- `mica/00_setup_env.sh` (HPC env setup)
- `mica/06_extract_claims.py` (validation tooling, not analysis)
- `mica/08_compile_report.py` (validation tooling)

**MICA dependency:** MICA is an R package. Must be installable via `install.packages()` or Bioconductor in the Docker container. If not available on CRAN/Bioconductor, provide precomputed MICA output and run scripts 04-05 (figure generation) only. **Decision point during implementation.**

### Section 02: Rat snRNA-seq (R)
**Included:**
- `01_qc_filter.R`, `02_normalize_integrate.R`, `03_cluster_annotate.R`, `03b_sctype_annotate.R`, `03c_sctype_cluster.R`, `08_differential_expression.R`, `09_differential_abundance.R`, `method_comparison.R`, `spot_check_markers.R`

**Excluded:** None (all 9 scripts included)

### Section 03_comparison: Cross-Species Comparison (R)
**Included:**
- `01_generate_comparison_report.R`

**Note:** This script compares `results/original/` vs `results/corrected/` directories (pre- vs post-biostatistical corrections). It does NOT depend on sections 01/02 outputs. Both `original/` and `corrected/` result sets must be included in the data asset under `data/03_comparison/original/` and `data/03_comparison/corrected/`. The `paths.R` config maps these explicitly rather than using the generic `resolve_input()` pattern.

### Section 03_rat_wes: Rat WES (Python)
**Included:**
- `01_parse_vep.py`, `02_cosmic_signatures.py`, `03_generate_oncoplot.py`

**Excluded:**
- `00_run_vep.sbatch` (requires Ensembl VEP cache)
- `run_analysis.sbatch` (SLURM wrapper)

### Section 04: Human scRNA-seq (R + Python)
**Included (pipeline scripts):**
- `00b_preprocess_seurat.R`, `01_load_subset_data.R`, `03_cell_type_annotation.R`, `04_cell_fractions.R`, `05_gene_expression_violin.R`, `06_run_gsva.R`, `07_run_progeny.R`, `08_run_wcsea.R`, `09_cellphonedb_prep.R`, `10_load_macrophage_seurat.R`, `11_macrophage_deg_analysis.R`, `12_macrophage_pathway_enrichment.R`, `12_enrichr_pathway_analysis.py`, `13_macrophage_cellphonedb.R`, `14_multicelltype_deg.R`, `15_multicelltype_pathway.R`, `16_cellchat_analysis.R`, `16_run_cellphonedb.py`, `17_cellphonedb_dotplot.R`

**Excluded:**
- `00a_download_geo.sh` (replaced with capsule-specific `00a_download_xu_geo.sh` for Xu et al. 2024 data)
- `15b_pathway_polarity_diagnostic.R` (diagnostic, not pipeline)
- `16b_run_cellphonedb_v4.py` (legacy v4 variant)
- `regenerate_fig7bc_fonts.R` (one-off font fix)
- All `sbatch_*.sh` and `run_analysis.sbatch` SLURM wrappers

**CellPhoneDB note:** Scripts 09, 13, 16, 17 depend on CellPhoneDB. Decision point: install v5 in Docker, or provide precomputed results + run visualization (script 17) only.

### Section 05: Rat Bulk RNA-seq (R)
**Included:**
- `03_deseq2.R`, `04_pam50_subtyping.R`, `05_validate_vs_rahul.R`

**Excluded:**
- `00_fastp_trim.sh`, `01_alignment.sh`, `02_htseq_count.sh` (alignment pipeline)
- `install_r_packages.sh` (HPC setup)
- `run_analysis.sh` (SLURM orchestrator)

### Section 06: Spatial Biopsies (Python)
**Included:**
- `01_immune_secretion.py`, `02_immune_pathways.py`, `figure_config.py`, `utils.py`

**Excluded:** None

### Section 07: Organoid Single Cell (Python)
**Included:**
- `02_qc.py`, `03_preprocess.py`, `04_pseudobulk_de.py`, `05_pathways.py`, `06_cell_cycle.py`, `07_single_cell_pathways.py`, `08_inhibitor_mechanism_exploration.py`, `09_pathway_enrichment_analysis.py`, `10_heterogeneity_analysis.py`, `11_proliferation_analysis.py`, `_config.py`, `_figure_config.py`

**Excluded:**
- `build_reference.sh` (Cell Ranger reference build)
- `01_cellranger_multi_pool1.sh`, `01_cellranger_multi_pool2.sh` (Cell Ranger)
- `run_analysis.sh` and all `submit_*.sh` SLURM wrappers

## Section Dependency Graph

```
03_comparison    (independent — reads results/original/ + results/corrected/)
01               (independent)
02               (independent)
03_rat_wes       (independent)
04               (independent — Xu et al. 2024 data)
05               (independent)
06               (independent — HCC22-088 spatial)
07               (independent)
```

**All sections are independent.** Section 03_comparison reads pre/post correction comparison data, NOT outputs from sections 01 or 02. This means `run_all.sh` can run them in any order, though sequential is simplest.

## Execution Modes

### Mode 1: `run_all.sh` — Full Computation

Runs all sections sequentially from count matrices. Each section writes intermediates to `results/intermediate/XX_outputs/`.

**Estimated runtimes** (approximate, depends on CodeOcean machine specs):

| Section | Estimated Runtime | Bottleneck |
|---------|------------------|------------|
| 01 | ~20 min | GSVA, PROGENy scoring |
| 02 | ~45 min | Seurat integration, DoubletFinder |
| 03_comparison | ~2 min | Simple CSV comparison |
| 03_rat_wes | ~10 min | Signature extraction |
| 04 | ~2-4 hours | Seurat SCT on 115K cells, CellChat, CellPhoneDB |
| 05 | ~15 min | DESeq2, PAM50 |
| 06 | ~10 min | Spatial pathway scoring |
| 07 | ~30 min | PyDESeq2, GSEA |
| **Total** | **~4-6 hours** | |

**Error handling:** Each section's `run.sh` uses `set -e` (fails fast on individual script error within a section). `run_all.sh` does NOT use `set -e` — instead it wraps each section call in a trap, logs success/failure, and continues to the next section since all are independent. Final exit code reflects whether any section failed.

### Mode 2: `run_codeocean.sh` — Figure Regeneration (~15-30 min)

Reads from precomputed data and regenerates all figures. In Tier 1, reads lightweight CSVs from `data/figure_regen/`. In Tier 2, can also read full RDS/h5ad from `data/precomputed/`.

### Mode 3: `run_section.sh <section_name>` — Single Section

Runs one section in isolation. Uses full directory name to avoid ambiguity between `03_comparison` and `03_rat_wes`:

```bash
./run_section.sh 03_comparison    # Not "./run_section.sh 03"
./run_section.sh 03_rat_wes
./run_section.sh 04_human_scrnaseq
```

### Input Resolution Logic

Every script resolves inputs through centralized config (`paths.R` / `paths.py`):

```
Priority:
1. results/intermediate/XX_outputs/   (from current run_all.sh run)
2. data/precomputed/XX_outputs/       (Tier 2 large objects, if present)
3. data/figure_regen/XX_outputs/      (Tier 1 lightweight CSVs)
```

Both `data/precomputed/` and `data/figure_regen/` exist in the directory structure. Tier 1 populates only `figure_regen/`. Tier 2 populates both.

Each script has a guard:
```r
# R pattern — centralized in paths.R as resolve_input()
resolve_input <- function(section, filename) {
  candidates <- c(
    file.path(RESULTS_DIR, "intermediate", section, filename),
    file.path(DATA_DIR, "precomputed", section, filename),
    file.path(DATA_DIR, "figure_regen", section, filename)
  )
  for (path in candidates) {
    if (file.exists(path)) return(path)
  }
  stop(paste("No input found for", section, filename))
}
```

**Pipeline markers:** The main repo uses `analysis/.pipeline_markers/` for SLURM job tracking. These are NOT used in the capsule. The capsule `run.sh` scripts handle sequencing directly.

## Two-Tier Data Strategy

Decision pending CodeOcean support response on expanded quotas.

### Tier 1: Lean (~2.3 GB compressed) — fits free academic quota

| Component | Compressed Size | Notes |
|-----------|----------------|-------|
| Sec 01: bulk counts + TPM + metadata | ~30 MB | GSE276755 |
| Sec 02: 6 x 10X matrices | ~600 MB | GSE276758 |
| Sec 03_comparison: original + corrected results | ~10 MB | Pre/post correction CSVs |
| Sec 03_wes: VEP outputs | ~1 MB | GSE276759 |
| Sec 04: Xu et al. 2024 mtx + metadata | ~1 GB | Compressed mtx (5:1 ratio verified) |
| Sec 05: htseq count matrices | ~1 MB | GSE276757 |
| Sec 06: spatial biopsy data | ~500 MB | HCC22-088, Carleton et al. |
| Sec 07: post-Cell Ranger matrices | ~80 MB | GEO pending |
| GMT files | ~1 MB | Estrogen signatures |
| Figure-regen CSVs (all sections) | ~50 MB | In data/figure_regen/ |
| **Total** | **~2.3 GB** | |

**Section 04 alternative:** If compressed Xu data is too large, download from GEO at runtime. The Xu et al. 2024 paper deposited data at GEO (accession TBD — need to check key resources table of doi:10.1016/j.xcrm.2024.101500). Code at github.com/ChanLab-UTSW/BreastCancer_Integrated.

**Figure-regen mode:** During `run_all.sh`, each section's plotting scripts export plot-ready CSVs/TSVs to `results/intermediate/`. The `prepare_data_asset.sh --tier 1` script copies these to `data/figure_regen/`. The `run_codeocean.sh` mode reads them to regenerate figures without heavy objects.

### Tier 2: Comfortable (~20-35 GB) — if expanded quota granted

Everything in Tier 1, plus:

| Additional Component | Size | Benefit |
|---------------------|------|---------|
| Sec 04 preprocessed Seurat RDS | ~7 GB | Skip Xu download + SCT normalization |
| Sec 02 integrated.rds | ~3 GB | Figure-regen with full Seurat object |
| Sec 04 CellChat/pathway objects | ~2 GB | Full intermediate reproducibility |
| Sec 07 precomputed h5ad | ~2 GB | |
| Sec 06 spatial uncompressed | ~2.5 GB | |
| **Total** | **~20-35 GB** | |

### Data Asset Packaging

```bash
# From ERpBRCA_OlderWomen/
./scripts/prepare_data_asset.sh --tier 1    # Lean
./scripts/prepare_data_asset.sh --tier 2    # Full
```

This script:
- Copies correct inputs from main repo
- Compresses where beneficial (mtx, pkl)
- Exports figure-regen CSVs from existing outputs (Tier 1)
- Copies precomputed RDS/h5ad (Tier 2)
- Generates SHA256 checksums for all files
- Reports total size
- Validates completeness (all sections have required inputs)

## Environment

Single Docker container with both R and Python.

**Base image:** `bioconductor/bioconductor_docker:3.19` (R 4.3.3, Ubuntu 22.04)

**Added via Dockerfile:** miniconda with Python 3.10

**R packages (postInstall):**
- Bioconductor: DESeq2, GSVA, PROGENy, SingleR, SingleCellExperiment, msigdbr, qusage, speckle
- CRAN: Seurat (>=5.0), harmony, DoubletFinder, CellChat, genefu, patchwork, pheatmap, ggplot2, data.table, lme4, dittoseq

**Python packages (postInstall):**
- scanpy, anndata, squidpy, gseapy, pydeseq2, decoupler, scrublet, statsmodels, adjustText, psutil, scikit-misc
- commot is NOT installed by default (torch dependency chain). Sec 06 scripts that require commot will use precomputed results. Check during implementation which Sec 06 scripts actually import commot vs just read its outputs.

**Potentially problematic packages (decision points during implementation):**

| Package | Issue | Fallback |
|---------|-------|----------|
| CellPhoneDB v5 | Complex dependencies, naming issues | Precomputed results + visualization only |
| commot | Requires torch, pyg, torch-scatter | Precomputed results (only used in Sec 06) |
| MICA | May not be on CRAN/Bioconductor | Precomputed MICA output + figure scripts |

**Excluded tools (not installable / not needed):**

| Tool | Reason | Data provided instead |
|------|--------|----------------------|
| Cell Ranger 9.0.1 | Proprietary, large install | Post-CR filtered matrices |
| STAR | Requires reference genome indices | htseq count matrices |
| VEP (Ensembl) | Requires annotation cache (~20 GB) | VEP-annotated outputs |
| SigProfilerMatrixGenerator | Complex install | Precomputed signatures |

## Excluded Analysis Steps

These steps are excluded from capsule execution but documented with exact parameters in the README:

| Step | Section | Tool | Why excluded |
|------|---------|------|-------------|
| FASTQ trimming | 05 | fastp | Requires raw FASTQs |
| STAR alignment | 05 | STAR 2.7.x | Requires mRatBN7.2 genome index |
| Feature counting | 05 | htseq-count | Requires aligned BAMs |
| Cell Ranger multi | 07 | Cell Ranger 9.0.1 | Requires raw FASTQs + reference |
| Reference build | 07 | Cell Ranger | Requires GRCh38 + probe set |
| VEP annotation | 03 | Ensembl VEP | Requires annotation cache |

## Sync Script: `scripts/sync_to_capsule.sh`

One-directional sync from main repo to capsule. Main repo is source of truth.

```bash
./scripts/sync_to_capsule.sh                          # Sync all sections
./scripts/sync_to_capsule.sh --section 04_human_scrnaseq  # Sync one section
./scripts/sync_to_capsule.sh --dry-run                # Preview changes
./scripts/sync_to_capsule.sh --diff                   # Show full diff
```

**Behavior:**

1. **Copies analysis scripts** from `analysis/XX_*/` to `ERpBRCA_CodeOcean/code/XX_*/`, using the per-section inclusion lists defined above
2. **Excludes** per the exclusion lists above: alignment scripts, SLURM wrappers (`sbatch_*.sh`, `submit_*.sh`, `run_analysis.sbatch`, `run_analysis.sh`), `.pipeline_markers/`, `__pycache__/`
3. **Replaces config files only** — swaps HPC-specific `_config.py`/paths with capsule-aware versions. Analysis scripts stay identical to main repo.
4. **Generates `code/MANIFEST.md`** — lists every script with source commit hash and modification date
5. **Validates completeness** — checks every script referenced in `run.sh` exists in capsule

**Does NOT:** copy data, modify analysis logic, auto-commit.

**Capsule location:** `../ERpBRCA_CodeOcean/` (sibling directory), configurable via `--target` or `CAPSULE_DIR` env var.

## README Structure

```
# Title + authors + journal
## Quick Start (3 commands: run_codeocean.sh, run_all.sh, run_section.sh)
## Analysis Sections (table: section, description, language, runtime, manuscript figures)
## Data Provenance (table: dataset, GEO accession, citation, preprocessing)
## Excluded Steps (table with exact commands/parameters for full reproducibility)
## Environment (R/Python versions, key packages with versions)
## Output Map (results/ path -> manuscript figure mapping)
## Seeding & Reproducibility (set.seed(12345), Section 07 exception: seed=42)
```

Each section's `run.sh` includes a header comment with: language, runtime estimate, inputs, outputs, manuscript figure references.

`data/README.md` documents provenance for every file with SHA256 checksums (auto-generated by `prepare_data_asset.sh`).

## Key Differences from CITEgeist Capsule

| Aspect | CITEgeist | ERpBRCA |
|--------|-----------|---------|
| Languages | Python only | R + Python |
| Sections | 6 | 7 (+ mica subpipeline) |
| Data size | 1.4 GB (pkl) | 2.3-35 GB (tier-dependent) |
| Docker base | conda only | bioconductor + conda |
| External deps | Gurobi (optional) | CellPhoneDB, commot, MICA (all TBD) |
| Section isolation | No (sequential only) | Yes (all sections independent) |
| Figure-regen data | Full pkl objects | Lightweight CSVs (Tier 1) or full objects (Tier 2) |
| Sync mechanism | Manual | `sync_to_capsule.sh` |
| Section runner | N/A | `run_section.sh <name>` |

## Open Questions for Implementation

1. **CellPhoneDB v5:** Clean install in Docker, or precomputed results only?
2. **commot:** Clean install (torch dependency chain), or precomputed results for Sec 06?
3. **MICA:** Available on CRAN/Bioconductor, or precomputed output + figure scripts?
4. **CodeOcean quota:** Tier 1 or Tier 2? Pending support response.
5. **Section 04 Xu data GEO accession:** Need to find from paper's key resources table (doi:10.1016/j.xcrm.2024.101500). NOT GSE176078 (old Wu data).
