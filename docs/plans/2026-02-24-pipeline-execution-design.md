# Pipeline Execution for Validation Prerequisites Design

**Date:** 2026-02-24
**Status:** Approved
**Author:** Claude (with user approval)

## Overview

Run missing analysis pipelines to generate outputs required by the manuscript validation pipeline. Uses marker-file polling to orchestrate parallel SLURM jobs and auto-trigger validation.

## Problem

The validation pipeline (`scripts/validation/`) expects specific CSV outputs that don't yet exist:
- `03_rat_wes`: No outputs directory
- `01_human_bulk_rnaseq`: Has RDS but missing `correlation_results.csv`, `mica_results.csv`, `hsd17b7_correlations.csv`
- `04_human_scrnaseq`: Missing `tam_analysis.csv`

## Scope

### Analyses to Run

| Analysis | Current Status | Action |
|----------|---------------|--------|
| `03_rat_wes` | No outputs | Run full pipeline |
| `01_human_bulk_rnaseq` | Has RDS | Run full pipeline + add CSV exports |
| `04_human_scrnaseq` | Mostly complete | Run full pipeline + add TAM CSV export |
| `02_rat_snrnaseq` | Complete | No action |
| `05_rat_bulk_rnaseq` | Complete | No action |

## Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                     PIPELINE ORCHESTRATION                           │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  sbatch_run_missing_analyses.sh                                     │
│       │                                                              │
│       ├──▶ sbatch 03_rat_wes/run_analysis.sh                        │
│       │         └──▶ .complete marker on success                    │
│       │                                                              │
│       ├──▶ sbatch 01_human_bulk_rnaseq/run_analysis.sh              │
│       │         └──▶ exports CSVs, writes .complete marker          │
│       │                                                              │
│       ├──▶ sbatch 04_human_scrnaseq/run_analysis.sh                 │
│       │         └──▶ exports TAM CSV, writes .complete marker       │
│       │                                                              │
│       └──▶ sbatch poll_and_validate.sbatch                          │
│                 │                                                    │
│                 ├── Poll for 3 .complete markers (60s interval)     │
│                 └── Run validation pipeline when all present        │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### Marker File Location

```
analysis/.pipeline_markers/
├── 03_rat_wes.complete
├── 01_human_bulk_rnaseq.complete
└── 04_human_scrnaseq.complete
```

### Polling Job

- Partition: HTC (minimal resources)
- Poll interval: 60 seconds
- Timeout: 4 hours
- On success: runs `scripts/validation/run_validation.sh`

## Implementation Strategy

**Principle: Minimal changes to existing scripts**

### Modifications to Existing Scripts

Each `run_analysis.sh` gets appended with:

1. **CSV Export Block** (for 01 and 04):
   - R code to read existing RDS outputs
   - Write required CSV files for validation

2. **Completion Marker**:
   ```bash
   MARKER_DIR="$(dirname "$0")/../.pipeline_markers"
   mkdir -p "$MARKER_DIR"
   touch "$MARKER_DIR/$(basename $(dirname "$0")).complete"
   ```

### New Scripts

| Script | Purpose | Location |
|--------|---------|----------|
| `sbatch_run_missing_analyses.sh` | Master launcher | `analysis/` |
| `poll_and_validate.sbatch` | Polling watcher | `analysis/` |

### SLURM Requirements

All sbatch scripts include:
```bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu
```

## Expected Outputs

After successful execution:

### 01_human_bulk_rnaseq/outputs/
- `correlation_results.csv` (new)
- `mica_results.csv` (new)
- `hsd17b7_correlations.csv` (new)

### 03_rat_wes/outputs/
- `cosmic_signatures.csv`
- `oncoplot_data.csv`

### 04_human_scrnaseq/outputs/
- `tam_analysis.csv` (new)

## Error Handling

- Each job writes marker only on successful completion
- Poll job logs missing markers every interval
- Poll job exits with error after 4-hour timeout
- SLURM failure emails sent to alc376@pitt.edu

## Success Criteria

1. All 3 analyses complete successfully (markers present)
2. All expected CSV outputs exist
3. Validation pipeline runs automatically
4. Validation report generated at `validation/report/validation_report.md`

## Notes

- Cross-partition SLURM dependencies don't work; marker-file polling handles this
- Existing RDS outputs are preserved; CSV exports are additive
- Validation pipeline already tested and functional
