#!/bin/bash
# =============================================================================
# processBAM.sh - BAM Processing Pipeline for TMSP and CEBPA
# =============================================================================
# This script automates BAM processing for TMSP and CEBNX panels.
# It generates QC statistics, coverage depth, and Pindel breakpoint analysis.
#
# Usage: Run from the bam directory containing TMSP BAM files
#   cd /path/to/analysis/bam
#   ~/Shared/SCRIPTS/claude/processBAM/processBAM.sh
#
# Expected directory structure:
#   /path/to/analysis/
#   ├── bam/                  <- Run script from here
#   │   └── *.bam            <- TMSP BAM files
#   ├── vcf/                  <- Required for Picard QC (variant stats)
#   │   └── *.vcf            <- TMSP VCF files
#   ├── cebpa/
#   │   └── bam/
#   │       └── *.bam        <- CEBPA BAM files (optional)
#   └── output/               <- Final output directory
#
# =============================================================================

set -e  # Exit on error

# =============================================================================
# CONFIGURATION
# =============================================================================

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Reference files
HG19="$HOME/Databases/WholeGenomeFASTA/genome.fa"
BED_TMSP="$HOME/Databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.NoZero.v2.bed"
BED_TMSP_REFSEQ="$HOME/Databases/TMSPvcf/BEDfiles/TMSP.RefSeqexons.bed"
BED_CEBNX="$HOME/Databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.CEBPA.bed"

# Conda environment
CONDA_PATH="$HOME/Software/miniconda3/etc/profile.d/conda.sh"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

activate_conda() {
    if [ -f "$CONDA_PATH" ]; then
        source "$CONDA_PATH"
        conda activate base
    else
        log_error "Conda not found at $CONDA_PATH"
        exit 1
    fi
}

# =============================================================================
# EXCEL WRITING FUNCTIONS (Integrated - no external scripts needed)
# =============================================================================
# These functions use Python with openpyxl to write Excel files

# Write two worksheets to Excel
# Usage: write_2ws_to_xlsx <output.xlsx> <ws1_name> <ws1_file> <ws2_name> <ws2_file>
write_2ws_to_xlsx() {
    local output_file="$1"
    local ws1_name="$2"
    local ws1_file="$3"
    local ws2_name="$4"
    local ws2_file="$5"

    python3 << PYTHON_EOF
import openpyxl
from openpyxl import Workbook

wb = Workbook()

# First worksheet
ws1 = wb.active
ws1.title = "${ws1_name}"
with open("${ws1_file}", 'r') as f:
    for row_idx, line in enumerate(f, 1):
        fields = line.rstrip('\n').split('\t')
        for col_idx, value in enumerate(fields, 1):
            # Try to convert to number if possible
            try:
                if '.' in value:
                    ws1.cell(row=row_idx, column=col_idx, value=float(value))
                else:
                    ws1.cell(row=row_idx, column=col_idx, value=int(value))
            except ValueError:
                ws1.cell(row=row_idx, column=col_idx, value=value)

# Second worksheet
ws2 = wb.create_sheet(title="${ws2_name}")
with open("${ws2_file}", 'r') as f:
    for row_idx, line in enumerate(f, 1):
        fields = line.rstrip('\n').split('\t')
        for col_idx, value in enumerate(fields, 1):
            try:
                if '.' in value:
                    ws2.cell(row=row_idx, column=col_idx, value=float(value))
                else:
                    ws2.cell(row=row_idx, column=col_idx, value=int(value))
            except ValueError:
                ws2.cell(row=row_idx, column=col_idx, value=value)

wb.save("${output_file}")
print(f"Created ${output_file}")
PYTHON_EOF
}

# Write single worksheet to Excel
# Usage: write_1ws_to_xlsx <output.xlsx> <ws_name> <data_file>
write_1ws_to_xlsx() {
    local output_file="$1"
    local ws_name="$2"
    local data_file="$3"

    python3 << PYTHON_EOF
import openpyxl
from openpyxl import Workbook

wb = Workbook()
ws = wb.active
ws.title = "${ws_name}"

with open("${data_file}", 'r') as f:
    for row_idx, line in enumerate(f, 1):
        fields = line.rstrip('\n').split('\t')
        for col_idx, value in enumerate(fields, 1):
            try:
                if '.' in value:
                    ws.cell(row=row_idx, column=col_idx, value=float(value))
                else:
                    ws.cell(row=row_idx, column=col_idx, value=int(value))
            except ValueError:
                ws.cell(row=row_idx, column=col_idx, value=value)

wb.save("${output_file}")
print(f"Created ${output_file}")
PYTHON_EOF
}

# =============================================================================
# PICARD QC FUNCTION
# =============================================================================
# Runs Picard QC metrics and mosdepth for all BAM files
# Output: QCdata.xlsx with Picard and Mosdepth worksheets

run_picard_qc() {
    local sample="$1"
    local sample_name=$(basename "$sample" .bam)

    # Check that corresponding VCF exists
    if [ ! -f "../vcf/$sample_name.vcf" ]; then
        log_error "VCF not found: ../vcf/$sample_name.vcf"
        return 1
    fi

    # Convert BAM to SAM
    [ -f "$sample_name.sam" ] || \
    picard SamFormatConverter \
        I="$sample_name.bam" \
        O="$sample_name.sam" \
        VALIDATION_STRINGENCY=SILENT

    # Collect quality metrics
    [ -f "$sample_name.quality.txt" ] || \
    picard CollectQualityYieldMetrics \
        I="$sample_name.sam" \
        O="$sample_name.quality.txt" \
        VALIDATION_STRINGENCY=SILENT

    # Collect alignment metrics
    [ -f "$sample_name.align.txt" ] || \
    picard CollectAlignmentSummaryMetrics \
        I="$sample_name.sam" \
        O="$sample_name.align.txt" \
        R="$HG19" \
        VALIDATION_STRINGENCY=SILENT

    # Combine metrics
    [ -f "$sample_name.metrics.txt" ] || \
    (
        echo "$sample_name" > "$sample_name.metrics.txt"
        awk '{if($1!~/^\#/)print}' "$sample_name.quality.txt" | transpose -t >> "$sample_name.metrics.txt"
        awk '{if($1!~/^\#/)print}' "$sample_name.align.txt" | transpose -t >> "$sample_name.metrics.txt"
        gatk3 -T VariantEval \
            -R "$HG19" \
            --eval "../vcf/$sample_name.vcf" \
            | awk -F"\t" '{if($1~/^(TiTv|Variant|Count)/)print}' >> "$sample_name.metrics.txt"
    )

    # Select required fields
    [ -f "$sample_name.SelectMetrics.txt" ] || \
    (
        cat \
            <(awk -F"\t" '{if(NR==1)print "Sample\t"$1}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==2)print $2"\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==14)print $2"_R1\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==14)print $2"_R2\t"$4}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==4)print $2"\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==5)print $2"\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==7)print $2"\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==9)print $2"\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==18)print $2"\t"$5}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==18)print $2"_R1\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==18)print $2"_R2\t"$4}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==20)print $2"\t"$5}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==20)print $2"_R1\t"$3}' "$sample_name.metrics.txt") \
            <(awk -F"\t" '{if(NR==20)print $2"_R2\t"$4}' "$sample_name.metrics.txt") \
            <(paste <(echo "Variants") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f7)) \
            <(paste <(echo "SNPs") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f12)) \
            <(paste <(echo "Insertions") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f14)) \
            <(paste <(echo "Deletions") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f15)) \
            <(paste <(echo "Heterozygous") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f20)) \
            <(paste <(echo "Homozygous") <(awk -F"\t" '{if(NR==49)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f22)) \
            <(paste <(echo "Transitions") <(awk -F"\t" '{if(NR==53)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f6)) \
            <(paste <(echo "Transversions") <(awk -F"\t" '{if(NR==53)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f7)) \
            <(paste <(echo "Ti/Tv") <(awk -F"\t" '{if(NR==53)print $1}' "$sample_name.metrics.txt" | sed 's/ \+ /\t/g' | cut -f8)) \
            > "$sample_name.SelectMetrics.txt"
    )

    # Get field 2 for collation
    [ -f "$sample_name.SelectMetrics.f2.txt" ] || \
    cut -f2 "$sample_name.SelectMetrics.txt" > "$sample_name.SelectMetrics.f2.txt"

    # Cleanup
    rm -f "$sample_name.quality.txt" "$sample_name.align.txt" "$sample_name.metrics.txt" "$sample_name.sam"
}

run_mosdepth() {
    local bam="$1"
    local sample=$(basename "$bam" .bam)

    # Run mosdepth
    mosdepth \
        --threshold 1,100,250,1000 \
        --fast-mode \
        --by "$BED_TMSP_REFSEQ" \
        --no-per-base \
        --threads 4 \
        "$sample" \
        "$bam"

    # Modify output
    gunzip "$sample".*.bed.gz 2>/dev/null || true

    paste <(cat "$BED_TMSP_REFSEQ") \
          <(cut -f5 "$sample.regions.bed" | sed "1i $sample") \
          <(cut -f5-8 "$sample.thresholds.bed") \
          > "$sample.mosdepthTemp.txt"

    cat <(awk '{if(NR==1)print}' "$sample.mosdepthTemp.txt") \
        <(awk '{if(NR!=1)print}' "$sample.mosdepthTemp.txt" | sort -k4,4) \
        > "$sample.mosdepth.txt"

    # Get fields for multisample paste
    cut -f7-11 "$sample.mosdepth.txt" > "$sample.mosdepth.f.txt"

    # Cleanup
    rm -f "$sample.mosdepth.global.dist.txt" \
          "$sample.mosdepth.region.dist.txt" \
          "$sample.mosdepth.summary.txt" \
          "$sample.regions.bed" \
          "$sample.regions.bed.gz.csi" \
          "$sample.thresholds.bed" \
          "$sample.thresholds.bed.gz.csi" \
          "$sample.mosdepthTemp.txt"
}

# Export functions for parallel
export -f run_picard_qc run_mosdepth
export HG19 BED_TMSP_REFSEQ

generate_qc_data() {
    log_info "######################## RUNNING PICARD QC #########################"

    [ -f QCdata.xlsx ] && { log_info "QCdata.xlsx already exists, skipping..."; return 0; }

    # Run Picard QC in parallel
    log_info "Running Picard metrics..."
    parallel --memfree 3G "run_picard_qc {}" ::: *.bam

    # Run mosdepth in parallel
    log_info "Running mosdepth..."
    parallel --memfree 3G "run_mosdepth {}" ::: *.bam

    # Get first sample name for collation
    local SAMPLE=$(basename $(ls *.bam | head -1) .bam)

    # Collate Picard files
    log_info "Collating Picard results..."
    paste \
        <(cut -f1 "$SAMPLE.SelectMetrics.txt") \
        <(paste *.SelectMetrics.f2.txt) \
        > Picard.ALL.txt

    # Collate mosdepth files
    log_info "Collating mosdepth results..."
    paste \
        <(cut -f1-6 "$SAMPLE.mosdepth.txt") \
        <(paste *.mosdepth.f.txt) \
        > Mosdepth.ALL.txt

    # Write to Excel (using integrated Python function)
    log_info "Writing QCdata.xlsx..."
    write_2ws_to_xlsx "QCdata.xlsx" "Picard.ALL" "Picard.ALL.txt" "Mosdepth.ALL" "Mosdepth.ALL.txt"

    # Cleanup
    rm -f *.SelectMetrics.txt *.SelectMetrics.f2.txt *.mosdepth.txt *.mosdepth.f.txt Picard.ALL.txt Mosdepth.ALL.txt

    log_info "QCdata.xlsx generated"
}

# =============================================================================
# COVERAGE DEPTH CALCULATION
# =============================================================================
# Calculates per-base depth and cross-validation statistics
# Output: CoverageData.xlsx with TMSP and CEBNX worksheets

calculate_tmsp_depth() {
    log_info "################# CALCULATING TMSP DEPTH #########################"

    [ -f Depth.TMSP.txt ] && { log_info "Depth.TMSP.txt already exists, skipping..."; return 0; }

    log_info "Calculating depth for each BAM..."
    local start=$(date +%s)
    parallel "[ -f {/.}.depth ] || samtools depth -aa -m 0 -b $BED_TMSP {} > {/.}.depth" ::: *.bam
    local end=$(date +%s)
    log_info "Depth calculation took $((end-start))s"

    # Add sample name to individual depths
    parallel "[ -f 1.{.}.field ] || cat <(echo {.}) <(cut -f3 {}) > 1.{.}.field" ::: *.depth
    # Get chr and pos
    parallel "[ -f 1.Chr ] || cat <(echo 'Chr') <(cut -f1 {}) > 1.Chr" ::: *.depth
    parallel "[ -f 1.Pos ] || cat <(echo 'Pos') <(cut -f2 {}) > 1.Pos" ::: *.depth

    # Concatenate together
    [ -f 2.Depth.txt ] || \
    paste <(cat 1.Chr) <(cat 1.Pos) <(paste *.field) > 2.Depth.txt

    # Cross-validation for each sample
    [ -f 2.Data.txt ] || cut -f3- 2.Depth.txt > 2.Data.txt
    [ -f 2.Data.content ] || awk '{if(NR!=1) print}' 2.Data.txt > 2.Data.content

    local num_fields=$(awk '{print NF}' 2.Data.txt | uniq)
    for ((field=1; field<=num_fields; field++)); do
        local sample_name=$(cut -f$field 2.Data.txt | awk '{if(NR==1)print}')
        log_info "Processing cross-validation for $sample_name..."

        # Get depth from all other samples
        [ -f "3.$sample_name.others.txt" ] || \
        cut -f$field --complement 2.Data.content > "3.$sample_name.others.txt"

        # Get depth for current sample
        [ -f "3.$sample_name.txt" ] || \
        cat <(echo "$sample_name") <(cut -f$field 2.Data.content) > "3.$sample_name.txt"

        # Calculate row average
        [ -f "3.$sample_name.ave" ] || \
        cat <(echo "$sample_name.AVE") \
            <(awk '{sum=cnt=0; for (i=1;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print (cnt ? sum/cnt : "NA") }' "3.$sample_name.others.txt") \
            > "3.$sample_name.ave"

        # Calculate row stdev
        [ -f "3.$sample_name.stdev" ] || \
        cat <(echo "$sample_name.STD") \
            <(awk '{ A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ;
                for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }' "3.$sample_name.others.txt") \
            > "3.$sample_name.stdev"

        # Paste values together
        [ -f "3.$sample_name.results" ] || \
        paste "3.$sample_name.txt" "3.$sample_name.ave" "3.$sample_name.stdev" > "3.$sample_name.results"
    done

    # Combine all samples
    [ -f Depth.TMSP.txt ] || \
    paste <(cut -f1-2 2.Depth.txt) <(paste *.results) > Depth.TMSP.txt

    # Write to Excel (using integrated Python function)
    [ -f CoverageTMSP.xlsx ] || \
    write_1ws_to_xlsx "CoverageTMSP.xlsx" "Depth.TMSP" "Depth.TMSP.txt"

    # Cleanup
    rm -f [1-3].* *.depth

    log_info "TMSP depth calculation complete"
}

calculate_cebnx_depth() {
    log_info "################# CALCULATING CEBNX DEPTH #########################"

    [ -f Depth.CEBNX.txt ] && { log_info "Depth.CEBNX.txt already exists, skipping..."; return 0; }

    # Process depth for each sample
    parallel "[ -f 1.{.}.depth ] || samtools depth -aa -m 0 -b $BED_CEBNX {} > 1.{.}.depth" ::: *.bam
    parallel "[ -f 1.{.}.field ] || cat <(echo {.}) <(cut -f3 {}) > 1.{.}.field" ::: *.depth
    parallel "[ -f 1.Chr ] || cat <(echo 'Chr') <(cut -f1 {}) > 1.Chr" ::: *.depth
    parallel "[ -f 1.Pos ] || cat <(echo 'Pos') <(cut -f2 {}) > 1.Pos" ::: *.depth

    # Concatenate together
    [ -f 2.Depth.txt ] || \
    paste <(cat 1.Chr) <(cat 1.Pos) <(paste *.field) > 2.Depth.txt

    # Cross-validation
    [ -f 2.Data.txt ] || cut -f3- 2.Depth.txt > 2.Data.txt
    [ -f 2.Data.content ] || awk '{if(NR!=1) print}' 2.Data.txt > 2.Data.content

    local num_fields=$(awk '{print NF}' 2.Data.txt | uniq)
    for ((field=1; field<=num_fields; field++)); do
        local sample_name=$(cut -f$field 2.Data.txt | awk '{if(NR==1)print}')
        log_info "Processing cross-validation for $sample_name..."

        [ -f "3.$sample_name.others.txt" ] || \
        cut -f$field --complement 2.Data.content > "3.$sample_name.others.txt"

        [ -f "3.$sample_name.txt" ] || \
        cat <(echo "$sample_name") <(cut -f$field 2.Data.content) > "3.$sample_name.txt"

        [ -f "3.$sample_name.ave" ] || \
        cat <(echo "$sample_name.AVE") \
            <(awk '{sum=cnt=0; for (i=1;i<=NF;i++) if ($i != "NA") { sum+=$i; cnt++ } print (cnt ? sum/cnt : "NA") }' "3.$sample_name.others.txt") \
            > "3.$sample_name.ave"

        [ -f "3.$sample_name.stdev" ] || \
        cat <(echo "$sample_name.STD") \
            <(awk '{ A=0; V=0; for(N=1; N<=NF; N++) A+=$N ; A/=NF ;
                for(N=1; N<=NF; N++) V+=(($N-A)*($N-A))/(NF-1); print sqrt(V) }' "3.$sample_name.others.txt") \
            > "3.$sample_name.stdev"

        [ -f "3.$sample_name.results" ] || \
        paste "3.$sample_name.txt" "3.$sample_name.ave" "3.$sample_name.stdev" > "3.$sample_name.results"
    done

    # Combine all samples
    [ -f Depth.CEBNX.txt ] || \
    paste <(cut -f1-2 2.Depth.txt) <(paste *.results) > Depth.CEBNX.txt

    # Write to Excel (using integrated Python function)
    [ -f CoverageCEBNX.xlsx ] || \
    write_1ws_to_xlsx "CoverageCEBNX.xlsx" "Depth.CEBNX" "Depth.CEBNX.txt"

    # Cleanup
    rm -f 1.* 2.* 3.*

    log_info "CEBNX depth calculation complete"
}

generate_coverage_data() {
    log_info "######################## GENERATING COVERAGE DATA #########################"

    [ -f CoverageData.xlsx ] && { log_info "CoverageData.xlsx already exists, skipping..."; return 0; }

    local start_dir=$(pwd)

    # Calculate TMSP depth
    calculate_tmsp_depth

    # Calculate CEBNX depth if directory exists
    if [ -d "../cebpa/bam" ]; then
        cd ../cebpa/bam
        calculate_cebnx_depth
        cp Depth.CEBNX.txt "$start_dir/"
        cd "$start_dir"

        # Write combined coverage data (using integrated Python function)
        log_info "Writing CoverageData.xlsx..."
        write_2ws_to_xlsx "CoverageData.xlsx" "Depth.TMSP" "Depth.TMSP.txt" "Depth.CEBNX" "Depth.CEBNX.txt"

        # Cleanup
        rm -f Depth.TMSP.txt Depth.CEBNX.txt CoverageTMSP.xlsx
    else
        log_info "CEBPA directory not found, skipping CEBNX depth calculation"
        # Just rename the TMSP coverage file
        mv CoverageTMSP.xlsx CoverageData.xlsx 2>/dev/null || true
    fi

    log_info "Coverage data generation complete"
}

# =============================================================================
# PINDEL BREAKPOINT DETECTION
# =============================================================================
# Runs Pindel for FLT3 (chr13) and CALR (chr19) genes
# Output: *.FLT3.pindel.vcf and *.CALR.pindel.vcf

run_pindel_flt3() {
    local sample="$1"
    local sample_name=$(basename "$sample" .bam)

    [ -f "$sample_name.FLT3.pindel.vcf" ] && return 0

    log_info "Running Pindel FLT3 for $sample_name..."

    # Create pindel directory
    [ ! -d ./pindel ] && mkdir -p ./pindel

    # Create config file
    echo "$sample_name.bam 400 $sample_name-FLT3" > "$sample_name-FLT3.bam_config.txt"

    # Run pindel
    pindel -f "$HG19" \
           -i "$sample_name-FLT3.bam_config.txt" \
           -o "$sample_name-FLT3" \
           -c chr13 \
           -T 2

    # Convert to VCF
    pindel2vcf -P "$sample_name-FLT3" \
               -r "$HG19" \
               -d 20190531 \
               -v "$sample_name.FLT3.pindel.vcf" \
               -R NCBI36 \
               -is 10 \
               -e 25

    # Move intermediate files
    mv "$sample_name-FLT3"* pindel/
}

run_pindel_calr() {
    local sample="$1"
    local sample_name=$(basename "$sample" .bam)

    [ -f "$sample_name.CALR.pindel.vcf" ] && return 0

    log_info "Running Pindel CALR for $sample_name..."

    # Create pindel directory
    [ ! -d ./pindel ] && mkdir -p ./pindel

    # Create config file
    echo "$sample_name.bam 400 $sample_name-CALR" > "$sample_name-CALR.bam_config.txt"

    # Run pindel
    pindel -f "$HG19" \
           -i "$sample_name-CALR.bam_config.txt" \
           -o "$sample_name-CALR" \
           -c chr19 \
           -T 8

    # Convert to VCF
    pindel2vcf -P "$sample_name-CALR" \
               -r "$HG19" \
               -d 20190531 \
               -v "$sample_name.CALR.pindel.vcf" \
               -R NCBI36 \
               -is 25 \
               -e 25

    # Move intermediate files
    mv "$sample_name-CALR"* pindel/
}

# Export functions for parallel
export -f run_pindel_flt3 run_pindel_calr
export HG19

run_pindel_analysis() {
    log_info "######################## RUNNING PINDEL ANALYSIS #########################"

    # Run FLT3 Pindel in parallel
    log_info "################# PINDEL FOR FLT3 #########################"
    parallel --memsuspend 3G "run_pindel_flt3 {}" ::: *.bam

    # Run CALR Pindel in parallel
    log_info "################# PINDEL FOR CALR #########################"
    parallel --memsuspend 3G "run_pindel_calr {}" ::: *.bam

    log_info "Pindel analysis complete"
}

# =============================================================================
# COVERAGE PLOT GENERATION
# =============================================================================
# Generates multi-page PDF coverage plots with gene shading
# Output: *_Coverage.pdf files (one per sample)

generate_coverage_plots() {
    log_info "######################## GENERATING COVERAGE PLOTS #########################"

    local plot_count=$(ls *_Coverage.pdf 2>/dev/null | wc -l)
    local bam_count=$(ls *.bam 2>/dev/null | wc -l)

    if [ "$plot_count" -ge "$bam_count" ] && [ "$bam_count" -gt 0 ]; then
        log_info "Coverage plots already exist, skipping..."
        return 0
    fi

    if [ ! -f "CoverageData.xlsx" ]; then
        log_error "CoverageData.xlsx not found. Run --coverage first."
        return 1
    fi

    log_info "Generating coverage plots from CoverageData.xlsx..."
    "$SCRIPT_DIR/generateCoveragePlots.sh" "CoverageData.xlsx" "."

    log_info "Coverage plot generation complete"
}

# =============================================================================
# CHECKPOINT DETECTION
# =============================================================================

check_qc_complete() {
    [ -f "QCdata.xlsx" ] && return 0
    return 1
}

check_coverage_complete() {
    [ -f "CoverageData.xlsx" ] && return 0
    return 1
}

check_pindel_complete() {
    local flt3_count=$(ls *.FLT3.pindel.vcf 2>/dev/null | wc -l)
    local calr_count=$(ls *.CALR.pindel.vcf 2>/dev/null | wc -l)
    local bam_count=$(ls *.bam 2>/dev/null | wc -l)

    [ "$flt3_count" -eq "$bam_count" ] && [ "$calr_count" -eq "$bam_count" ] && return 0
    return 1
}

check_plots_complete() {
    local plot_count=$(ls *_Coverage.pdf 2>/dev/null | wc -l)
    local bam_count=$(ls *.bam 2>/dev/null | wc -l)

    [ "$plot_count" -ge "$bam_count" ] && [ "$bam_count" -gt 0 ] && return 0
    return 1
}

show_status() {
    echo ""
    echo "Pipeline Status:"
    echo "================"

    if check_qc_complete; then
        echo "  [✓] QC data complete (QCdata.xlsx)"
    else
        echo "  [ ] QC data not complete"
    fi

    if check_coverage_complete; then
        echo "  [✓] Coverage data complete (CoverageData.xlsx)"
    else
        echo "  [ ] Coverage data not complete"
    fi

    if check_pindel_complete; then
        local vcf_count=$(ls *.pindel.vcf 2>/dev/null | wc -l)
        echo "  [✓] Pindel analysis complete ($vcf_count VCF files)"
    else
        echo "  [ ] Pindel analysis not complete"
    fi

    if check_plots_complete; then
        local pdf_count=$(ls *_Coverage.pdf 2>/dev/null | wc -l)
        echo "  [✓] Coverage plots complete ($pdf_count PDF files)"
    else
        echo "  [ ] Coverage plots not complete"
    fi
    echo ""
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

main() {
    local start_time=$(date +%s)

    # Stage flags
    local run_qc=false
    local run_coverage=false
    local run_pindel=false
    local run_plots=false
    local run_all=true
    local force=false
    local show_status_only=false

    log_info "=========================================="
    log_info "BAM Processing Pipeline (TMSP + CEBPA)"
    log_info "=========================================="
    log_info "Working directory: $(pwd)"
    log_info "Script directory: $SCRIPT_DIR"

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --qc)
                run_qc=true
                run_all=false
                ;;
            --coverage)
                run_coverage=true
                run_all=false
                ;;
            --pindel)
                run_pindel=true
                run_all=false
                ;;
            --plots)
                run_plots=true
                run_all=false
                ;;
            --force|-f)
                force=true
                ;;
            --status)
                show_status_only=true
                ;;
            --check)
                log_info "Checking dependencies..."
                # Activate conda for tool checking
                if [ -f "$CONDA_PATH" ]; then
                    source "$CONDA_PATH"
                    conda activate base 2>/dev/null || true
                fi
                echo ""
                echo "Required tools:"
                for cmd in samtools picard mosdepth pindel pindel2vcf gatk3 parallel transpose; do
                    if command -v "$cmd" &> /dev/null; then
                        echo "  $cmd: OK ($(which $cmd))"
                    else
                        echo "  $cmd: MISSING"
                    fi
                done
                echo ""
                echo "Reference files:"
                for ref in "$HG19" "$BED_TMSP" "$BED_TMSP_REFSEQ" "$BED_CEBNX"; do
                    if [ -f "$ref" ]; then
                        echo "  $ref: OK"
                    else
                        echo "  $ref: MISSING"
                    fi
                done
                echo ""
                echo "Python modules:"
                if python3 -c "import openpyxl" 2>/dev/null; then
                    echo "  openpyxl: OK"
                else
                    echo "  openpyxl: MISSING (pip install openpyxl)"
                fi
                exit 0
                ;;
            --help|-h)
                echo "Usage: $0 [options]"
                echo "  Run from the bam directory containing TMSP BAM files"
                echo ""
                echo "Pipeline Stages:"
                echo "  1. QC        - Run Picard and mosdepth QC metrics"
                echo "  2. Coverage  - Calculate coverage depth with cross-validation"
                echo "  3. Pindel    - Detect FLT3 and CALR breakpoints"
                echo "  4. Plots     - Generate coverage plots with gene shading"
                echo ""
                echo "Options:"
                echo "  (no options)    Run all stages (skip completed stages)"
                echo "  --status        Show pipeline completion status"
                echo "  --qc            Run QC stage only"
                echo "  --coverage      Run coverage stage only"
                echo "  --pindel        Run Pindel stage only"
                echo "  --plots         Run coverage plot generation only"
                echo "  --force, -f     Force re-run even if stage is complete"
                echo "  --check         Check dependencies"
                echo "  --help          Show this help"
                echo ""
                echo "Examples:"
                echo "  $0                  # Run full pipeline"
                echo "  $0 --status         # Check what stages are complete"
                echo "  $0 --qc --force     # Force regenerate QC data"
                echo "  $0 --plots          # Generate coverage plots only"
                echo ""
                echo "Expected directory structure:"
                echo "  /path/to/analysis/"
                echo "  ├── bam/                  <- Run script from here"
                echo "  ├── vcf/                  <- Required for QC (variant stats)"
                echo "  ├── cebpa/bam/            <- CEBPA BAM files (optional)"
                echo "  └── output/               <- Created by script"
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
        shift
    done

    # Show status only
    if [ "$show_status_only" = true ]; then
        show_status
        exit 0
    fi

    # Check for BAM files
    local bam_count=$(ls *.bam 2>/dev/null | wc -l)
    if [ "$bam_count" -eq 0 ]; then
        log_error "No BAM files found in current directory"
        exit 1
    fi
    log_info "Found $bam_count BAM files"

    # Activate conda environment
    log_info "Activating conda environment..."
    activate_conda

    # Rename files if >20 characters
    rename 's/^(.{20}).*.bam/$1\.bam/' *.bam 2>/dev/null || true

    # Index BAM files if needed
    log_info "Checking BAM indexes..."
    for bam in *.bam; do
        if [ ! -f "${bam}.bai" ] && [ ! -f "${bam%.*}.bai" ]; then
            log_info "Indexing $bam..."
            samtools index "$bam"
        fi
    done

    # Determine which stages to run
    if [ "$run_all" = true ]; then
        run_qc=true
        run_coverage=true
        run_pindel=true
        run_plots=true
    fi

    # =========================================================================
    # STAGE 1: QC DATA
    # =========================================================================
    if [ "$run_qc" = true ]; then
        if [ "$force" = true ] || ! check_qc_complete; then
            log_info ">>> STAGE 1: QC DATA"
            generate_qc_data
        else
            log_info ">>> STAGE 1: QC DATA [SKIPPED - already complete]"
        fi
    fi

    # =========================================================================
    # STAGE 2: COVERAGE DATA
    # =========================================================================
    if [ "$run_coverage" = true ]; then
        if [ "$force" = true ] || ! check_coverage_complete; then
            log_info ">>> STAGE 2: COVERAGE DATA"
            generate_coverage_data
        else
            log_info ">>> STAGE 2: COVERAGE DATA [SKIPPED - already complete]"
        fi
    fi

    # =========================================================================
    # STAGE 3: PINDEL ANALYSIS
    # =========================================================================
    if [ "$run_pindel" = true ]; then
        if [ "$force" = true ] || ! check_pindel_complete; then
            log_info ">>> STAGE 3: PINDEL ANALYSIS"
            run_pindel_analysis
        else
            log_info ">>> STAGE 3: PINDEL ANALYSIS [SKIPPED - already complete]"
        fi
    fi

    # =========================================================================
    # STAGE 4: COVERAGE PLOTS
    # =========================================================================
    if [ "$run_plots" = true ]; then
        if [ "$force" = true ] || ! check_plots_complete; then
            log_info ">>> STAGE 4: COVERAGE PLOTS"
            generate_coverage_plots
        else
            log_info ">>> STAGE 4: COVERAGE PLOTS [SKIPPED - already complete]"
        fi
    fi

    # =========================================================================
    # CREATE OUTPUT DIRECTORY
    # =========================================================================
    log_info "######################## CREATING OUTPUT #########################"
    [ ! -d ../output ] && mkdir -p ../output
    cp *.xlsx ../output/ 2>/dev/null || true
    cp *_Coverage.pdf ../output/ 2>/dev/null || true

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))

    log_info "=========================================="
    log_info "######################## DONE #########################"
    log_info "Pipeline completed in ${duration}s"

    # Show final status
    show_status

    log_info "=========================================="
}

# Run main function
main "$@"
