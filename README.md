# processBAM Pipeline

A fully self-contained BAM processing pipeline for TMSP and CEBPA/CEBNX panels. This pipeline automates QC statistics, coverage depth calculation, and Pindel breakpoint analysis for FLT3 and CALR genes.

## Installation

### Clone the Repository

```bash
# Clone via HTTPS
git clone https://github.com/alvin8-git/processBAM.git

# Or clone via SSH
git clone git@github.com:alvin8-git/processBAM.git

# Navigate to the directory
cd processBAM
```

## Overview

This pipeline:
- Runs Picard and mosdepth QC metrics on BAM files
- Calculates coverage depth with cross-validation statistics
- Detects FLT3 (chr13) and CALR (chr19) breakpoints using Pindel
- Generates Excel reports for QC and coverage data
- **Generates coverage plot PDFs** with alternating gene shading (using R)
- **Generates QC report PDFs** with Picard stats and gene coverage tables (using R)
- **Fully self-contained** - all Excel writing functions integrated (no external scripts)

## Directory Structure

```
processBAM/
├── processBAM.sh              # Main pipeline script (self-contained)
├── generateCoveragePlots.sh   # Coverage plot generation wrapper
├── plotCoverage.R             # R script for coverage plot generation
├── generateQCreport.sh        # QC report generation wrapper
├── plotQC.R                   # R script for QC report generation
├── Dockerfile                 # Docker container definition
├── docker-compose.yml         # Docker Compose configuration
└── README.md                  # This file
```

### Output Directory Structure

After running the pipeline, the output directory will contain:

```
output/
├── QCdata.xlsx                # QC metrics Excel file
├── CoverageData.xlsx          # Coverage data Excel file
├── pindel/                    # Pindel structural variant results
│   ├── *.FLT3.pindel.vcf     # FLT3 breakpoint VCFs
│   ├── *.CALR.pindel.vcf     # CALR breakpoint VCFs
│   └── intermediate/          # Pindel intermediate files
├── CoveragePlots/             # Coverage plot PDFs
│   └── *_Coverage.pdf
└── QCreports/                 # QC report PDFs
    └── *_QC.pdf
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

### Python Modules

```bash
pip install openpyxl
```

### Reference Files

These should be present in `$HOME/Databases/`:
- `WholeGenomeFASTA/genome.fa` - hg19 reference genome
- `TMSPvcf/BEDfiles/TSMP.UCSCexons.NoZero.v2.bed` - TMSP exons BED
- `TMSPvcf/BEDfiles/TMSP.RefSeqexons.bed` - TMSP RefSeq exons
- `TMSPvcf/BEDfiles/TSMP.UCSCexons.CEBPA.bed` - CEBPA exons BED

### Conda Environment (Native)

The pipeline requires conda with a `base` environment containing the bioinformatics tools:

```bash
# Install conda environment with required tools
conda create -n base samtools picard mosdepth pindel gatk
```

### Docker (Recommended)

A Docker container is provided with all tools pre-installed. This is the recommended way to run the pipeline.

```bash
# Build the Docker image
cd processBAM
docker build -t processbam:latest .

# Run with docker-compose (recommended)
cd /path/to/analysis
DATABASES_PATH=~/Databases ANALYSIS_PATH=. docker-compose \
    -f ~/path/to/processBAM/docker-compose.yml \
    run --rm processbam

# Or run directly with docker
docker run --rm \
    -v ~/Databases:/databases:ro \
    -v /path/to/analysis:/data \
    -w /data/bam \
    processbam:latest
```

**Enter the container interactively:**

```bash
# Start an interactive shell inside the container
docker run --rm -it \
    -v ~/Databases:/databases:ro \
    -v /path/to/analysis:/data \
    -w /data/bam \
    --entrypoint bash \
    processbam:latest

# Once inside the container, you can run commands directly:
samtools --version
picard --version
mosdepth --version
pindel
/scripts/processBAM.sh --help
```

**Docker volume mounts:**
- `/databases` - Reference files (read-only): genome.fa, BED files
- `/data` - Analysis directory (read-write): bam/, vcf/, output/

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
  4. Plots     - Generate coverage plots with gene shading
  5. QCReport  - Generate QC report PDFs (Picard + gene coverage)

Options:
  (no options)    Run all stages (skip completed stages)
  --status        Show pipeline completion status
  --qc            Run QC stage only
  --coverage      Run coverage stage only
  --pindel        Run Pindel stage only
  --plots         Run coverage plot generation only
  --qcreport      Run QC report generation only
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

# Generate coverage plots only
./processBAM.sh --plots

# Generate QC reports only
./processBAM.sh --qcreport

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

### Pindel Output (`../output/pindel/`)

Per-sample VCF files for structural variants:
- `*.FLT3.pindel.vcf` - FLT3 breakpoints (chr13)
- `*.CALR.pindel.vcf` - CALR breakpoints (chr19)
- `intermediate/` - Intermediate Pindel files

### Coverage Plots (`*_Coverage.pdf`)

Multi-page PDF files per sample with:
- **TMSP pages (5 pages)**: ~13,700 positions per page covering all 56 TMSP genes
- **CEBPA page (1 page)**: ~1,200 positions for CEBPA gene

Each plot includes:
- **Sample depth** (dark blue line) - actual coverage at each position
- **UCL/LCL lines** (light grey) - confidence interval (AVE +/- STD)
- **250x reference line** (red dashed) - minimum coverage threshold
- **Z-score subplot** (pink) - deviation from mean (-2 to 2 range)
- **Ratio subplot** (purple) - coverage relative to mean (0 to 2 range)
- **Alternating gene shading** (cyan/white) - to demarcate gene boundaries
- **Gene labels** - rotated 90 degrees at top of each region

### QC Reports (`*_QC.pdf`)

Two-page PDF files per sample with:

**Page 1 - Picard Statistics:**
- Sequencing statistics (total reads, R1/R2, read length)
- Base quality (total bases, Q20, Q30)
- Alignment statistics (reads aligned, bases aligned)
- Variant statistics (SNPs, insertions, deletions, Ti/Tv ratio)
- Zygosity (heterozygous, homozygous, Het/Hom ratio)

**Page 2 - Gene Coverage Table:**
- Per-gene coverage statistics for all 54 TMSP genes
- Exons targeted per gene
- Average coverage depth
- Proportion of bases covered at 1X, 100X, 250X, >1000X
- Color-coded cells (red=0% to white=100% coverage)
- Overall summary row

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

### Plots Stage

1. Extracts TMSP and CEBPA data from `CoverageData.xlsx` to CSV
2. Uses R base graphics to generate multi-page PDF plots
3. Splits 68,800 TMSP positions across 5 pages
4. Creates alternating gene shading for all 56 genes
5. Adds sample depth, UCL/LCL lines, Z-score, and Ratio subplots

### QC Reports Stage

1. Extracts Picard.ALL and Mosdepth.ALL from `QCdata.xlsx` to CSV
2. Uses R base graphics to generate 2-page PDF reports
3. Page 1: Formats Picard metrics into readable statistics sections
4. Page 2: Aggregates mosdepth data by gene and creates coverage table
5. Applies color scale to coverage proportions (red=0% to white=100%)

## Integrated Functions

The following functions are built into the script (no external dependencies):

### Excel Writing
- `write_2ws_to_xlsx` - Creates Excel file with 2 worksheets from TSV files
- `write_1ws_to_xlsx` - Creates Excel file with 1 worksheet from TSV file
- Uses Python with openpyxl (automatically handles numeric conversion)

### QC Processing
- `run_picard_qc` - Runs Picard metrics for a single sample
- `run_mosdepth` - Runs mosdepth coverage for a single sample
- `generate_qc_data` - Orchestrates QC for all samples

### Coverage Processing
- `calculate_tmsp_depth` - Calculates TMSP coverage with cross-validation
- `calculate_cebnx_depth` - Calculates CEBNX coverage with cross-validation
- `generate_coverage_data` - Orchestrates coverage for all samples

### Pindel Processing
- `run_pindel_flt3` - Runs Pindel for FLT3 gene
- `run_pindel_calr` - Runs Pindel for CALR gene
- `run_pindel_analysis` - Orchestrates Pindel for all samples

### Coverage Plot Generation
- `generate_coverage_plots` - Orchestrates plot generation for all samples
- `generateCoveragePlots.sh` - Wrapper script for Excel to CSV conversion and R plotting
- `plotCoverage.R` - R script using base graphics for PDF generation

### QC Report Generation
- `generate_qc_reports` - Orchestrates QC report generation for all samples
- `generateQCreport.sh` - Wrapper script for Excel to CSV conversion and R plotting
- `plotQC.R` - R script using base graphics for 2-page PDF reports

## Checkpoint System

The pipeline automatically detects completed stages:
- Checks for `QCdata.xlsx` for QC stage
- Checks for `CoverageData.xlsx` for coverage stage
- Checks for `*.pindel.vcf` files in `../output/pindel/` for Pindel stage
- Checks for `*_Coverage.pdf` files in `../output/CoveragePlots/` for Plots stage
- Checks for `*_QC.pdf` files in `../output/QCreports/` for QC Reports stage

Use `--status` to view completion status, `--force` to re-run completed stages.

## Troubleshooting

### Check Dependencies
```bash
./processBAM.sh --check
```

### Common Issues

**Native installation:**
1. **Conda not found**: Ensure conda is installed at `$HOME/Software/miniconda3/`
2. **VCF not found**: QC stage requires matching VCF files in `../vcf/`
3. **Missing tools**: Install required tools via conda: `conda install samtools picard mosdepth pindel`
4. **openpyxl missing**: Run `pip install openpyxl`

**Docker:**
1. **Volume mount errors**: Ensure paths exist and have correct permissions
2. **Memory issues**: Increase Docker memory limit (default 16GB recommended)
3. **Database not found**: Mount your Databases directory to `/databases`

## Version History

- v1.4 (2025-01): Docker support and output organization
  - Added Docker container with all dependencies pre-installed
  - Dockerfile and docker-compose.yml for containerized execution
  - Reorganized Pindel output to `../output/pindel/` directory
  - Environment variable support for configurable paths
  - Docker mode auto-detection
- v1.3 (2025-01): QC report generation
  - Added Stage 5: QC report PDFs (Picard + gene coverage)
  - R script for 2-page PDF generation with color-coded coverage table
  - Per-gene coverage aggregation for 54 TMSP genes
  - New `--qcreport` command line option
- v1.2 (2025-01): Coverage plot generation
  - Added Stage 4: Coverage plots with gene shading
  - R script for multi-page PDF generation using base graphics
  - Alternating cyan/white gene shading for 56 TMSP genes
  - Sample depth, UCL/LCL lines, Z-score and Ratio subplots
  - New `--plots` command line option
- v1.1 (2025-01): Self-contained Excel writers
  - Integrated Excel writing functions using Python/openpyxl
  - Removed external Perl script dependencies
  - Single-file deployment
- v1.0 (2025-01): Initial release
  - Picard and mosdepth QC metrics
  - Coverage depth with cross-validation
  - Pindel FLT3/CALR breakpoint detection
  - Checkpoint system with stage-based execution
  - Parallel processing using GNU parallel

## Author

Alvin Ng
