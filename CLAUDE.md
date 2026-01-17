# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

processBAM is a BAM processing pipeline for TMSP and CEBPA/CEBNX panels. It automates QC statistics, coverage depth calculation with cross-validation, Pindel structural variant detection (FLT3/CALR), and generates PDF reports.

## Running the Pipeline

```bash
# Navigate to BAM directory and run
cd /path/to/analysis/bam
~/Shared/SCRIPTS/claude/processBAM/processBAM.sh

# Check status of completed stages
./processBAM.sh --status

# Run specific stages
./processBAM.sh --qc           # Picard/mosdepth metrics
./processBAM.sh --coverage     # Depth with cross-validation
./processBAM.sh --pindel       # FLT3/CALR breakpoints
./processBAM.sh --plots        # Coverage plot PDFs
./processBAM.sh --qcreport     # QC report PDFs

# Force re-run a completed stage
./processBAM.sh --qc --force

# Check dependencies
./processBAM.sh --check
```

### Docker Execution

```bash
# Build
docker build -t processbam:latest .

# Run via docker-compose
DATABASES_PATH=~/Databases ANALYSIS_PATH=/path/to/analysis \
  docker-compose run --rm processbam
```

## Architecture

### Pipeline Stages (in order)

1. **QC** - Picard metrics + mosdepth coverage thresholds → `QCdata.xlsx`
2. **Coverage** - Per-base depth with leave-one-out cross-validation → `CoverageData.xlsx`
3. **Pindel** - Structural variants for FLT3 (chr13) and CALR (chr19) → `*.pindel.vcf`
4. **Plots** - Multi-page coverage PDFs with gene shading → `*_Coverage.pdf`
5. **QCReport** - 2-page QC PDFs (Picard stats + gene coverage table) → `*_QC.pdf`

### File Responsibilities

| File | Purpose |
|------|---------|
| `processBAM.sh` | Main orchestrator with all stages, checkpoint detection, Excel writers |
| `generateCoveragePlots.sh` | Excel→CSV extraction + calls R for plot generation |
| `generateQCreport.sh` | Excel→CSV extraction + calls R for QC report generation |
| `plotCoverage.R` | Generates multi-page PDFs with gene-shaded coverage plots |
| `plotQC.R` | Generates 2-page PDFs with Picard stats and coverage tables |
| `transpose.c` | Text transpose utility (compiled into Docker) |

### Key Patterns

- **Checkpoint system**: Stages skip if output exists (e.g., `QCdata.xlsx`); use `--force` to re-run
- **Cross-validation**: For each sample, computes AVE/STD from all other samples at each position
- **Parallel processing**: Uses GNU `parallel --memsuspend 3G` for BAM operations
- **Excel writing**: Embedded Python/openpyxl in bash heredocs (no external Perl scripts)
- **Docker mode**: Auto-detected via `/.dockerenv`; skips conda activation

### Expected Directory Structure

```
/path/to/analysis/
├── bam/           ← Run script from here (TMSP BAM files)
├── vcf/           ← Required for QC (matching VCF files)
├── cebpa/bam/     ← Optional CEBPA BAM files
└── output/        ← Created by script
    ├── QCdata.xlsx
    ├── CoverageData.xlsx
    ├── pindel/
    ├── CoveragePlots/
    └── QCreports/
```

## Dependencies

**Tools**: samtools, picard, mosdepth, pindel, pindel2vcf, gatk3, parallel, transpose, R

**Python**: openpyxl

**Reference files** (via environment variables or `~/Databases/`):
- `HG19`: WholeGenomeFASTA/genome.fa
- `BED_TMSP`: TSMP.UCSCexons.NoZero.v2.bed
- `BED_TMSP_REFSEQ`: TMSP.RefSeqexons.bed
- `BED_CEBNX`: TSMP.UCSCexons.CEBPA.bed

## Conventions

- Sample names truncated to 20 characters via `rename` command
- BAM files auto-indexed if `.bai` missing
- Intermediate files cleaned up after each stage
- Coverage plot gene boundaries defined in `plotCoverage.R` (56 TMSP genes)
- QC report gene list in `plotQC.R` (54 TMSP genes)
