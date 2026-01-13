# processBAM Pipeline

A self-contained BAM processing pipeline for TMSP and CEBPA/CEBNX panels. This pipeline automates QC statistics, coverage depth calculation, and Pindel breakpoint analysis for FLT3 and CALR genes.

## Overview

This pipeline:
- Runs Picard and mosdepth QC metrics on BAM files
- Calculates coverage depth with cross-validation statistics
- Detects FLT3 (chr13) and CALR (chr19) breakpoints using Pindel
- Generates Excel reports for QC and coverage data

## Directory Structure

```
processBAM/
├── processBAM.sh              # Main pipeline script
├── write2WStoXLS.pl           # Excel writer (2 worksheets)
├── writeTMSPDepthtoXLS.pl     # TMSP coverage Excel writer
├── writeCEBNXDepthtoXLS.pl    # CEBNX coverage Excel writer
└── README.md                  # This file
```

## Prerequisites

### System Tools

The following tools must be installed (available via conda):

| Tool | Purpose |
|------|---------|
| samtools | BAM manipulation and indexing |
| picard | QC metrics collection |
| mosdepth | Coverage depth calculation |
| pindel | Structural variant detection |
| pindel2vcf | Pindel VCF conversion |
| gatk3 | Variant evaluation |
| parallel | GNU Parallel for job distribution |
| transpose | Text transpose utility |

### Perl Modules

```bash
cpan install Excel::Writer::XLSX
```

### Reference Files

These should be present in `$HOME/Databases/`:
- `WholeGenomeFASTA/genome.fa` - hg19 reference genome
- `TMSPvcf/BEDfiles/TSMP.UCSCexons.NoZero.v2.bed` - TMSP exons BED
- `TMSPvcf/BEDfiles/TMSP.RefSeqexons.bed` - TMSP RefSeq exons
- `TMSPvcf/BEDfiles/TSMP.UCSCexons.CEBPA.bed` - CEBPA exons BED

### Conda Environment

The pipeline requires conda with a `base` environment containing the bioinformatics tools:

```bash
# Install conda environment with required tools
conda create -n base samtools picard mosdepth pindel gatk
```

## Usage

### Expected Input Directory Structure

```
/path/to/analysis/
├── bam/                  <- Run script from here
│   └── *.bam            <- TMSP BAM files
├── vcf/                  <- Required for QC (variant stats)
│   └── *.vcf            <- TMSP VCF files (matching BAM names)
├── cebpa/
│   └── bam/
│       └── *.bam        <- CEBPA BAM files (optional)
└── output/               <- Created by script
```

### Running the Pipeline

1. Navigate to the BAM directory:
   ```bash
   cd /path/to/analysis/bam
   ```

2. Run the pipeline:
   ```bash
   ~/Shared/SCRIPTS/claude/processBAM/processBAM.sh
   ```

### Command Line Options

```
Usage: processBAM.sh [options]

Pipeline Stages:
  1. QC        - Run Picard and mosdepth QC metrics
  2. Coverage  - Calculate coverage depth with cross-validation
  3. Pindel    - Detect FLT3 and CALR breakpoints

Options:
  (no options)    Run all stages (skip completed stages)
  --status        Show pipeline completion status
  --qc            Run QC stage only
  --coverage      Run coverage stage only
  --pindel        Run Pindel stage only
  --force, -f     Force re-run even if stage is complete
  --check         Check dependencies
  --help          Show this help
```

### Examples

```bash
# Run full pipeline
./processBAM.sh

# Check pipeline status
./processBAM.sh --status

# Force regenerate QC data
./processBAM.sh --qc --force

# Run only Pindel analysis
./processBAM.sh --pindel

# Check dependencies
./processBAM.sh --check
```

## Output

### QC Data (`QCdata.xlsx`)

Excel file with two worksheets:

**Picard.ALL worksheet:**
- Sample name
- Total reads, R1/R2 reads
- PF reads, PF bases
- Q20 and Q30 bases
- Alignment statistics
- Variant counts (SNPs, insertions, deletions)
- Ti/Tv ratio

**Mosdepth.ALL worksheet:**
- Chromosome, start, end, gene, transcript, exon
- Per-sample mean coverage
- Coverage thresholds (1X, 100X, 250X, 1000X)

### Coverage Data (`CoverageData.xlsx`)

Excel file with two worksheets:

**Depth.TMSP worksheet:**
- Per-base coverage for all TMSP exons
- Per-sample depth with cross-validation statistics (AVE, STD)
- ~68,800 positions

**Depth.CEBNX worksheet:**
- Per-base coverage for CEBPA exons
- Per-sample depth with cross-validation statistics
- ~1,200 positions

### Pindel Output

Per-sample VCF files for structural variants:
- `*.FLT3.pindel.vcf` - FLT3 breakpoints (chr13)
- `*.CALR.pindel.vcf` - CALR breakpoints (chr19)
- `pindel/` - Intermediate Pindel files

## Processing Details

### QC Stage

1. Runs `picard SamFormatConverter` to convert BAM to SAM
2. Runs `picard CollectQualityYieldMetrics` for read quality
3. Runs `picard CollectAlignmentSummaryMetrics` for alignment stats
4. Runs `gatk3 VariantEval` for variant statistics
5. Runs `mosdepth` for coverage thresholds
6. Collates all samples into single Excel file

### Coverage Stage

1. Calculates per-base depth using `samtools depth`
2. Performs cross-validation: for each sample, computes average and stdev of all other samples at each position
3. Writes TMSP and CEBNX coverage to separate worksheets

### Pindel Stage

1. Creates config file for each sample
2. Runs `pindel` for FLT3 (chr13) and CALR (chr19)
3. Converts output to VCF using `pindel2vcf`
4. Stores intermediate files in `pindel/` directory

## Checkpoint System

The pipeline automatically detects completed stages:
- Checks for `QCdata.xlsx` for QC stage
- Checks for `CoverageData.xlsx` for coverage stage
- Checks for `*.pindel.vcf` files for Pindel stage

Use `--status` to view completion status, `--force` to re-run completed stages.

## Troubleshooting

### Check Dependencies
```bash
./processBAM.sh --check
```

### Common Issues

1. **Conda not found**: Ensure conda is installed at `$HOME/Software/miniconda3/`
2. **VCF not found**: QC stage requires matching VCF files in `../vcf/`
3. **Missing tools**: Install required tools via conda: `conda install samtools picard mosdepth pindel`
4. **Excel::Writer::XLSX missing**: Run `cpan install Excel::Writer::XLSX`

## Version History

- v1.0 (2025-01): Initial release
  - Picard and mosdepth QC metrics
  - Coverage depth with cross-validation
  - Pindel FLT3/CALR breakpoint detection
  - Checkpoint system with stage-based execution
  - Parallel processing using GNU parallel

## Author

Alvin Ng
