#!/bin/bash
# generateCoveragePlots.sh - Generate coverage plot PDFs from CoverageData.xlsx
#
# Usage: generateCoveragePlots.sh [CoverageData.xlsx] [output_dir]
#
# This script:
# 1. Extracts TMSP and CEBPA data from Excel to CSV (using Python)
# 2. Generates multi-page PDFs per sample with gene-shaded coverage plots (using R)
#
# Output: One PDF per sample with:
#   - 5 pages of TMSP coverage (~68k positions split across pages)
#   - 1 page of CEBPA coverage (~1.2k positions)
#   - Alternating shaded backgrounds for each gene
#   - Sample depth (black), AVE (blue), UCL/LCL bands, 250x reference line

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default paths
COVERAGE_FILE="${1:-CoverageData.xlsx}"
OUTPUT_DIR="${2:-.}"

# Activate conda for Python and R
CONDA_PATH="$HOME/Software/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_PATH" ]; then
    source "$CONDA_PATH"
    conda activate base 2>/dev/null || true
fi

# Check inputs
if [ ! -f "$COVERAGE_FILE" ]; then
    echo "Error: Coverage file not found: $COVERAGE_FILE"
    echo "Usage: $0 [CoverageData.xlsx] [output_dir]"
    exit 1
fi

echo "=== Coverage Plot Generation ==="
echo "Input: $COVERAGE_FILE"
echo "Output: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Temporary CSV files
TMSP_CSV="$OUTPUT_DIR/.tmsp_coverage.csv"
CEBPA_CSV="$OUTPUT_DIR/.cebpa_coverage.csv"

# Step 1: Extract Excel to CSV using Python
echo "Step 1: Extracting Excel data to CSV..."
python3 << PYTHON_EOF
import openpyxl
import csv
import sys

coverage_file = "${COVERAGE_FILE}"
tmsp_csv = "${TMSP_CSV}"
cebpa_csv = "${CEBPA_CSV}"

print(f"  Reading: {coverage_file}")
wb = openpyxl.load_workbook(coverage_file, read_only=True, data_only=True)

# Extract TMSP data
if 'Depth.TMSP' in wb.sheetnames:
    ws = wb['Depth.TMSP']
    print(f"  Extracting Depth.TMSP ({ws.max_row} rows)...")
    with open(tmsp_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in ws.iter_rows(values_only=True):
            writer.writerow(row)
    print(f"  Wrote: {tmsp_csv}")
else:
    print("  Warning: Depth.TMSP sheet not found")
    sys.exit(1)

# Extract CEBPA/CEBNX data
cebpa_sheet = None
for name in ['Depth.CEBNX', 'Depth.CEBPA']:
    if name in wb.sheetnames:
        cebpa_sheet = name
        break

if cebpa_sheet:
    ws = wb[cebpa_sheet]
    print(f"  Extracting {cebpa_sheet} ({ws.max_row} rows)...")
    with open(cebpa_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in ws.iter_rows(values_only=True):
            writer.writerow(row)
    print(f"  Wrote: {cebpa_csv}")
else:
    print("  Note: CEBPA/CEBNX sheet not found, skipping")
    # Create empty file
    open(cebpa_csv, 'w').close()

wb.close()
print("  Excel extraction complete")
PYTHON_EOF

# Step 2: Generate plots using R
echo ""
echo "Step 2: Generating coverage plots..."
Rscript "$SCRIPT_DIR/plotCoverage.R" "$TMSP_CSV" "$CEBPA_CSV" "$OUTPUT_DIR"

# Cleanup temporary files
echo ""
echo "Cleaning up temporary files..."
rm -f "$TMSP_CSV" "$CEBPA_CSV"

echo ""
echo "=== Coverage plot generation complete ==="
echo "PDF files created in: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"/*.pdf 2>/dev/null || echo "No PDF files found"
