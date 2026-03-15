# Methods: Rat Whole Exome Sequencing Analysis (`03_rat_wes`)

## 1. Data Acquisition

Whole exome sequencing (WES) was performed on mammary tumors from six rats, comprising three old (sample IDs: 102, 107, 116) and three young (sample IDs: 157, 158, 167) animals. Matched spleen tissue served as the germline reference for each animal. Somatic variants were called using Mutect2 (upstream of this pipeline); input VCF files were located in a shared directory (`NeilRatWES/RatWES_Mutect2_VCF_Input/`). Sample naming in the input VCFs followed the format `Rat_O_<ID>` (old) or `Rat_Y_<ID>` (young), which was parsed via regex to extract the three-digit numeric sample identifier.

> **[WARNING]** Hardcoded external input path for VCF files (`/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input`) reduces portability. This path is duplicated in both `00_run_vep.sbatch` and `02_cosmic_signatures.py`.

## 2. Variant Annotation with VEP

Somatic variants were annotated using Ensembl Variant Effect Predictor (VEP) v114.2 (module `ensembl-vep/114.2`), run in offline mode against a local cache (`00_run_vep.sbatch`). VEP was executed with the following parameters:

| Parameter | Value |
|-----------|-------|
| `--species` | `rattus_norvegicus` |
| `--assembly` | `Rnor_6.0` |
| `--cache_version` | `95` |
| `--offline` | enabled |
| `--force_overwrite` | enabled |
| `--tab` | tab-delimited output |
| `--fork` | `4` (parallel threads) |

Output files were named `<ID>_tumor_vs_<ID>_spleen.mutect2.txt` and organized into per-sample subdirectories. The SLURM job requested 4 CPUs and 16 GB memory with a 2-hour time limit.

> **[WARNING]** VEP cache version 95 (Ensembl 95, January 2019) paired with VEP software v114.2 represents a substantial version mismatch. Cache version 95 corresponds to the Rnor_6.0 assembly, but annotations may lack updates from later Ensembl releases. This should be documented as a limitation or the cache should be updated.

## 3. VEP Output Parsing and Filtering

Annotated variants were parsed and filtered using a custom Python script (`01_parse_vep.py`). The VEP tab-delimited output was read with custom header parsing: the script identified the header line as the first line beginning with a single `#` (but not `##`), stripped the leading `#`, and used the resulting fields as column names. All subsequent non-comment lines were loaded via `pd.read_csv()` with `comment='#'`.

Variants were retained if they met two criteria:
1. **Ensembl gene filter**: The `Gene` column value started with the prefix `'ENS'` (string match via `str.startswith()`).
2. **Impact filter**: The `IMPACT` field (either as a standalone column or extracted from the `Extra` column) was `'HIGH'` or `'MODERATE'`.

Only subfolders whose names ended in `'spleen'` (matching the `<ID>_tumor_vs_<ID>_spleen` convention) were processed. Filtered variants from each sample were saved as individual CSV files (`<sample>_filtered.csv`) and concatenated into a combined file (`all_samples_filtered.csv`) using `pd.concat()` with sample names as hierarchical index keys.

> **[WARNING]** Python scripts do not set random seeds despite project requirement (`set.seed(12345)` / `np.random.seed(12345)`). While `01_parse_vep.py` is deterministic (no stochastic operations), this violates the project-wide reproducibility standard stated in `CLAUDE.md`.

## 4. Mutational Signature Analysis

COSMIC mutational signatures were assigned using SigProfilerAssignment (`02_cosmic_signatures.py`). The `Analyzer.cosmic_fit()` function was called with the following parameters:

| Parameter | Value |
|-----------|-------|
| `input_type` | `"vcf"` |
| `context_type` | `"96"` (trinucleotide context, SBS) |
| `genome_build` | `"rn6"` (rat genome, Rnor_6.0) |
| `cosmic_version` | `3.4` |

The function performed signature decomposition by fitting observed trinucleotide mutation spectra to the COSMIC v3.4 single base substitution (SBS) reference signatures. Input VCF files were read directly from the shared input directory.

After fitting, the script parsed the SigProfiler output file (`Assignment_Solution/Activities/Assignment_Solution_Activities.txt`), removed signatures with zero activity across all samples, and annotated each sample with its age group (Old or Young) based on a hardcoded ID-to-group mapping. Results were saved to `cosmic_signatures.csv`.

A stacked bar chart was generated showing signature activities per sample. Samples were ordered with Young first (descending numeric ID: 167, 158, 157) then Old (descending: 116, 107, 102). Display labels appended an age suffix (e.g., `"157_Y"`, `"102_O"`). The figure used the `tab20` colormap, `figsize=(14, 8)`, axis label font size 18, tick font size 14/12, and was saved as both SVG and PNG (300 DPI).

> **[WARNING]** Young/old sample suffix labeling is inconsistent between scripts, risking interpretation errors. In `02_cosmic_signatures.py`, Young samples (157, 158, 167) receive the suffix `_Y` and Old samples (102, 107, 116) receive `_O`. However, in `03_generate_oncoplot.py` (line 217), the comment states the opposite convention: "Per legend: Young (102, 107, 116) = '_Y', Old (157, 158, 167) = '_O'". The actual code in `03_generate_oncoplot.py` assigns `_Y` to samples in the `YOUNG_SAMPLES` list (157, 158, 167), which is consistent with `02_cosmic_signatures.py`. The misleading comment could cause confusion during review or future edits.

> **[WARNING]** Python scripts do not set random seeds despite project requirement. SigProfilerAssignment may use non-deterministic optimization internally; without explicit seed setting, exact reproducibility of signature activities is not guaranteed.

> **[WARNING]** Hardcoded external input path for COSMIC reference data reduces portability. The VCF input path (`/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/NeilRatWES/RatWES_Mutect2_VCF_Input`) is hardcoded rather than parameterized.

## 5. Oncoplot Generation

An oncoplot was generated from the filtered variant data (`03_generate_oncoplot.py`) through the following steps:

### 5.1 Rat-to-Human Gene Mapping

Rat Ensembl gene IDs were mapped to human gene symbols via a two-step BioMart query using the `pybiomart` Python package:

1. **Homolog lookup**: Rat gene IDs (`rnorvegicus_gene_ensembl` dataset) were queried for human homologs (`hsapiens_homolog_ensembl_gene` attribute). Queries were batched in chunks of 200 genes to avoid API timeouts.
2. **Symbol resolution**: The resulting human Ensembl IDs were queried against the `hsapiens_gene_ensembl` dataset for `external_gene_name` attributes, also in chunks of 200.

Results were cached to `homolog_cache.csv` to avoid redundant API calls on re-runs. Genes without a human homolog or without a resolved symbol were excluded.

> **[WARNING]** BioMart queries depend on the live Ensembl REST API (`http://www.ensembl.org`), which introduces a non-reproducibility risk. The BioMart database is updated with each Ensembl release; gene symbol mappings and homolog assignments may change over time. The cache file mitigates this for repeated runs but the initial query results are not version-locked.

### 5.2 Cancer Gene Filtering

Mapped variants were filtered against a curated breast cancer gene list (`data/brca_genelist.csv`, 229 genes plus header). The gene list contained two columns: `Gene` (human gene symbol) and `Cancer` (cancer type annotation, e.g., `"PANCAN"`). Only variants whose human gene symbol appeared in this list were retained for the oncoplot.

### 5.3 Oncoplot Construction

The oncoplot was constructed as follows:

1. **Impact prioritization**: When a gene harbored multiple variants in the same sample, the variant with the highest impact was retained (HIGH > MODERATE, priority values 1 and 2 respectively). Rows were sorted by `['Sample', 'gene_symbol', 'IMPACT_PRIORITY']` and the first entry per gene-sample pair was kept via `groupby().first()`.
2. **Matrix pivoting**: Data were pivoted to a gene-by-sample matrix with VEP consequence types as cell values (`pivot_table()` with `aggfunc='first'`).
3. **Gene ordering**: Genes were sorted by mutation frequency (number of samples with a non-null entry) in descending order.
4. **Sample ordering**: Samples were ordered with Young first (descending numeric ID), then Old (descending numeric ID).
5. **Visualization**: The consequence matrix was encoded numerically (unique consequence types mapped to integers 1..N, 0 for wild-type), rendered as a heatmap via `seaborn.heatmap()` with `tab20`-derived colormap (white for wild-type), 0.5-width gray gridlines, and a categorical legend. Figure height scaled dynamically as `max(num_genes * 0.4, 6)` with fixed width of 10. Output was saved as SVG and PNG (300 DPI).

X-axis labels displayed sample IDs with age suffixes (e.g., `"167_Y"`, `"102_O"`).

> **[WARNING]** The script uses `applymap()` (line 187), which was deprecated in pandas 2.1.0 in favor of `map()`. This will raise a `FutureWarning` in current pandas versions and will fail in future releases.

## 6. Pipeline Orchestration

The analysis pipeline was orchestrated via `run_analysis.sbatch`, which executed steps 01 through 03 sequentially under the `aging_wes` conda environment. The SLURM job requested 8 CPUs, 32 GB memory, and an 8-hour time limit. VEP annotation (`00_run_vep.sbatch`) was submitted as a separate job (4 CPUs, 16 GB, 2-hour limit) and was a prerequisite for the main pipeline. Upon completion, a marker file (`.pipeline_markers/03_rat_wes.complete`) was written with the job ID and UTC timestamp.

> **[WARNING]** No explicit dependency is enforced between `00_run_vep.sbatch` and `run_analysis.sbatch`. The pipeline relies on the user submitting them in the correct order and waiting for VEP completion, rather than using SLURM `--dependency=afterok:<jobid>`.

## 7. Software Versions

From `environment.yml` (environment name: `erp_brca_aging`):

| Software | Version | Source |
|----------|---------|--------|
| Python | 3.10 | conda-forge |
| pandas | (not pinned) | conda-forge |
| numpy | (not pinned) | conda-forge |
| matplotlib | (not pinned) | conda-forge |
| seaborn | (not pinned) | conda-forge |

Additional dependencies not in `environment.yml` but required by the scripts:

| Software | Version | Notes |
|----------|---------|-------|
| Ensembl VEP | 114.2 | Loaded via `module load ensembl-vep/114.2` |
| VEP cache | 95 | Rnor_6.0, offline cache |
| SigProfilerAssignment | (not pinned) | `cosmic_fit()` with COSMIC v3.4 |
| pybiomart | (not pinned) | BioMart API queries |

> **[WARNING]** Key Python dependencies (`SigProfilerAssignment`, `pybiomart`) are not listed in the project `environment.yml`. The pipeline uses a separate conda environment (`aging_wes`) whose specification is not tracked in this repository. This impedes reproducibility.

## 8. Reproducibility Notes

- **Random seeds**: No random seeds were set in any Python script, contrary to the project requirement of `np.random.seed(12345)`.
- **Caching**: BioMart homolog queries were cached to `outputs/homolog_cache.csv`, making subsequent runs independent of external API state. SigProfiler was skipped on re-runs if its output file already existed.
- **Determinism**: VEP annotation and VEP parsing are deterministic. SigProfilerAssignment's internal optimization may not be deterministic without explicit seed control. BioMart query results may vary across Ensembl releases.
- **Completion markers**: The pipeline wrote a `.pipeline_markers/03_rat_wes.complete` file upon successful completion, enabling downstream dependency checks.
- **Output formats**: Figures were saved in both SVG (vector) and PNG (300 DPI raster). Intermediate data were saved as CSV.

---

## Main Text Summary

Somatic variants from whole exome sequencing of six rat mammary tumors (three young, three old; matched spleen germline controls) were annotated using Ensembl VEP v114.2 (Rnor_6.0 assembly, cache v95) and filtered for HIGH or MODERATE impact consequences on Ensembl-annotated genes. Mutational signatures were decomposed against COSMIC v3.4 SBS references using SigProfilerAssignment with 96-trinucleotide context on the rn6 genome. For oncoplot visualization, rat genes were mapped to human orthologs via BioMart and filtered against a curated panel of 229 cancer-associated genes. When multiple variants affected the same gene in a sample, the highest-impact consequence was retained. Samples were grouped by age and ordered by descending numeric ID within each group.
