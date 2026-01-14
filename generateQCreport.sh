#!/bin/bash
# generateQCreport.sh - Generate QC report PDFs from QCdata.xlsx
#
# Usage: generateQCreport.sh [QCdata.xlsx] [output_dir]
#
# This script:
# 1. Extracts Picard.ALL and Mosdepth.ALL data from Excel to CSV (using Python)
# 2. Generates 2-page PDFs per sample with QC statistics (using R)
#
# Output: One PDF per sample with:
#   - Page 1: Sequencing, Alignment, and Variant statistics (Picard)
#   - Page 2: Gene coverage table with color-coded proportions (Mosdepth)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default paths
QC_FILE="${1:-QCdata.xlsx}"
OUTPUT_DIR="${2:-.}"

# Activate conda for Python and R
CONDA_PATH="$HOME/Software/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_PATH" ]; then
    source "$CONDA_PATH"
    conda activate base 2>/dev/null || true
fi

# Check inputs
if [ ! -f "$QC_FILE" ]; then
    echo "Error: QC data file not found: $QC_FILE"
    echo "Usage: $0 [QCdata.xlsx] [output_dir]"
    exit 1
fi

echo "=== QC Report Generation ==="
echo "Input: $QC_FILE"
echo "Output: $OUTPUT_DIR"
echo ""

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Temporary CSV files
PICARD_CSV="$OUTPUT_DIR/.picard_data.csv"
MOSDEPTH_CSV="$OUTPUT_DIR/.mosdepth_data.csv"

# Step 1: Extract Excel to CSV using Python
echo "Step 1: Extracting Excel data to CSV..."
python3 << PYTHON_EOF
import openpyxl
import csv
import sys

qc_file = "${QC_FILE}"
picard_csv = "${PICARD_CSV}"
mosdepth_csv = "${MOSDEPTH_CSV}"

print(f"  Reading: {qc_file}")
wb = openpyxl.load_workbook(qc_file, read_only=True, data_only=True)

# Extract Picard.ALL data
if 'Picard.ALL' in wb.sheetnames:
    ws = wb['Picard.ALL']
    print(f"  Extracting Picard.ALL ({ws.max_row} rows, {ws.max_column} cols)...")
    with open(picard_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in ws.iter_rows(values_only=True):
            writer.writerow(row)
    print(f"  Wrote: {picard_csv}")
else:
    print("  Error: Picard.ALL sheet not found")
    sys.exit(1)

# Extract Mosdepth.ALL data
if 'Mosdepth.ALL' in wb.sheetnames:
    ws = wb['Mosdepth.ALL']
    print(f"  Extracting Mosdepth.ALL ({ws.max_row} rows, {ws.max_column} cols)...")
    with open(mosdepth_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        for row in ws.iter_rows(values_only=True):
            writer.writerow(row)
    print(f"  Wrote: {mosdepth_csv}")
else:
    print("  Error: Mosdepth.ALL sheet not found")
    sys.exit(1)

wb.close()
print("  Excel extraction complete")
PYTHON_EOF

# Step 2: Generate QC reports using R
echo ""
echo "Step 2: Generating QC report PDFs..."
Rscript "$SCRIPT_DIR/plotQC.R" "$PICARD_CSV" "$MOSDEPTH_CSV" "$OUTPUT_DIR"

# Cleanup temporary files
echo ""
echo "Cleaning up temporary files..."
rm -f "$PICARD_CSV" "$MOSDEPTH_CSV"

echo ""
echo "=== QC report generation complete ==="
echo "PDF files created in: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"/*_QC.pdf 2>/dev/null || echo "No QC PDF files found"
