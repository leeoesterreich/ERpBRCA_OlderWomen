---
name: pipeline-checkpointing-design-2026-06-28
description: Split monolithic scRNA-seq pipeline into two checkpointed stages to avoid reloading 37GB Seurat object unnecessarily
metadata:
  type: project
  status: approved
---

# Pipeline Checkpointing: Two-Stage Split for scRNA-seq Analysis

## Problem

`run_analysis.sbatch` runs Steps 01-11 sequentially in a single 12h job. The 37GB Seurat object (`SeuratObj_Xu2024_HRpos_AfterQCSCT.rds`) is loaded multiple times across steps. When the job fails mid-pipeline (e.g., Step 06 GSVA), the entire job restarts from Step 01 — reloading the Seurat object unnecessarily, wasting ~2h of I/O, and risking another crash.

Job 10338426 demonstrated this: Steps 01-04 completed successfully, but the job disappeared before reaching Step 06. The preprocessed outputs (`seurat_annotated.rds` at 640MB) are perfectly valid but must be regenerated on re-run.

## Design

Split the monolithic pipeline into two independent sbatch scripts with checkpoint-based dependency:

### Stage 1: Preprocessing (run_preprocessing.sbatch)

- **Resource requests**: 12h wall-clock, 256G RAM, 8 CPUs, HTC cluster
- **Steps**: 01 (load & subset), 03 (cell type annotation), 04 (cell fractions), 05 (gene expression)
- **Key output**: `seurat_annotated.rds` (640MB) — this is the checkpoint
- **Marker**: `.pipeline_markers/04_human_scrnaseq_preprocessing.complete`

### Stage 2: Downstream (run_downstream.sbatch)

- **Resource requests**: 4h wall-clock, 16G RAM, 4 CPUs, HTC cluster
- **Steps**: 06 (GSVA pseudobulk), 07 (PROGENy), 08 (WCSEA), 09 (CellPhoneDB prep), 10 (TAM export), 11 (ensure PNG figures)
- **Dependency**: Checks that `seurat_annotated.rds` exists before starting
- **Marker**: `.pipeline_markers/04_human_scrnaseq.complete` (same as current — backwards compatible)

### Orchestrator (run_analysis.sbatch)

- **Kept for backwards compatibility** — runs both stages sequentially
- Checks preprocessing marker: if already done, skips Stage 1
- After Stage 1 completes, checks downstream marker: if already done, skips Stage 2

## Components

### File Changes

| File | Action | Purpose |
|------|--------|---------|
| `run_preprocessing.sbatch` | CREATE | Stage 1: Steps 00b-05 |
| `run_downstream.sbatch` | CREATE | Stage 2: Steps 06-11 |
| `run_analysis.sbatch` | MODIFY | Orchestrator with checkpoint logic |

### Data Flow

1. **Input**: `SeuratObj_Xu2024_HRpos_AfterQCSCT.rds` (37GB, shared storage) → Step 01
2. Step 01 output → `seurat_young_midage_elderly.rds` (32GB)
3. Step 03 output → `seurat_annotated.rds` (640MB) — **CHECKPOINT**
4. Step 06 reads `seurat_annotated.rds` for metadata + the 37GB file for raw counts
5. Steps 07-11 read intermediates from `outputs/` directory between stages

### Error Handling

- Each stage validates input files exist before starting its R scripts
- Stage 2 can be submitted independently — it checks for `seurat_annotated.rds` and exits with code 1 if missing
- The orchestrator checks markers before dispatching; existing markers skip completed stages
- If Stage 1 fails, the orchestrator exits immediately with code 1 (no cascading to Stage 2)
- Completion markers written only at the end of each stage — partial failures do not leave stale markers

## Testing Strategy

1. **Stage 2 standalone**: Submit `run_downstream.sbatch` directly. `seurat_annotated.rds` already exists from the previous run (Steps 01-04 completed before the crash). If it produces GSVA outputs, the split works.
2. **Full pipeline**: After Stage 2 is verified, test the orchestrator by removing completion markers and submitting `run_analysis.sbatch`.
3. **Resume scenario**: Let Stage 1 complete, then kill the downstream job mid-Step-06. Re-submit Stage 2 and verify it uses existing `seurat_annotated.rds`.

## Success Criteria

- Stage 2 produces identical GSVA outputs to the original monolithic pipeline
- Stage 2 completes in under 30 minutes (pseudobulk GSVA on 37 samples is fast)
- Orchestrator correctly skips completed stages
- Preexisting outputs from partial runs are preserved, not regenerated

## What's NOT Changing

- R scripts (01_load_subset_data.R, 03_cell_type_annotation.R, 06_run_gsva.R, etc.) remain untouched
- Output file paths remain the same
- Figure formats and naming remain the same
- `00b_preprocess_seurat.R` behavior (checkpoint on preprocessed data existence) is preserved
