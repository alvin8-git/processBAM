#!/usr/bin/env Rscript
# plotCoverage.R - Generate coverage plots from CSV files using base R graphics
# Creates multi-page PDF with alternating gene shading for TMSP and CEBPA panels
#
# Usage: Rscript plotCoverage.R tmsp.csv cebpa.csv output_dir
#
# Input: CSV files with columns: Chr, Pos, Sample1, Sample1.AVE, Sample1.STD, ...
# Output: One PDF per sample with TMSP (5 pages) and CEBPA (1 page) coverage plots
#
# Each page contains:
#   - Coverage plot (log scale Y-axis) with UCL/LCL bands
#   - Z-score subplot (pink)
#   - Ratio subplot (purple)

# =============================================================================
# Gene boundary definitions (from TMSP panel)
# =============================================================================
tmsp_genes <- data.frame(
  gene = c("CSF3R", "MPL", "NRAS", "DNMT3A", "SF3B1", "IDH1", "MYD88", "CBLB",
           "GATA2", "GATA2-AS1", "PDGFRA", "KIT", "TET2", "FBXW7", "NPM1",
           "IKZF1", "CUX1", "BRAF", "EZH2", "RAD21", "JAK2", "CDKN2A",
           "CDKN2B-AS1", "ABL1", "NOTCH1", "PTEN", "SMC3", "HRAS", "WT1",
           "MLL", "CBL", "ETV6", "KRAS", "PTPN11", "FLT3", "IDH2", "TP53",
           "SRSF2", "SETBP1", "CALR", "JAK3", "CEBPA", "CBLC", "ASXL1",
           "GNAS", "RUNX1", "U2AF1", "ZRSR2", "BCOR", "KDM6A", "GATA1",
           "SMC1A", "ATRX", "STAG2", "BCORL1", "PHF6"),
  start = c(1, 1132, 1229, 1536, 4374, 5025, 5317, 5832, 6168, 6742, 7658,
            8025, 9072, 15141, 15760, 15920, 17578, 22039, 22158, 24437,
            26356, 26572, 27123, 27305, 27841, 30139, 30545, 31317, 31660,
            31904, 33008, 33344, 34798, 35099, 35446, 35884, 36045, 37266,
            37656, 37844, 38047, 38132, 38898, 39072, 42066, 42237, 43805,
            44027, 45502, 50812, 55118, 55357, 56014, 58439, 62373, 67647),
  stringsAsFactors = FALSE
)
# Calculate end positions
tmsp_genes$end <- c(tmsp_genes$start[-1] - 1, 68801)

# CEBPA gene boundary
cebpa_genes <- data.frame(
  gene = "CEBPA",
  start = 1,
  end = 1217,
  stringsAsFactors = FALSE
)

# =============================================================================
# Color definitions
# =============================================================================
COL_GENE_SHADE <- rgb(0, 1, 1, 0.3)       # Cyan with transparency
COL_GENE_SHADE_ALT <- "white"              # White for alternating
COL_COVERAGE <- "darkblue"                 # Bold dark blue for coverage
COL_AVERAGE <- "blue"                      # Blue for average line
COL_UCL_LCL <- rgb(1, 1, 0.6, 0.5)        # Light yellow for UCL/LCL
COL_ZSCORE <- rgb(1, 0.4, 0.7, 0.8)       # Pink for Z-score
COL_RATIO <- rgb(0.6, 0.2, 0.8, 0.8)      # Purple for Ratio
COL_REF_LINE <- "red"                      # Red for 250x reference

# =============================================================================
# Helper functions
# =============================================================================

# Parse sample names from column headers
get_sample_names <- function(df) {
  cols <- colnames(df)[-(1:2)]  # Skip Chr, Pos
  # Get base sample names (every 3rd column starting from 1)
  sample_cols <- cols[seq(1, length(cols), by = 3)]
  # Clean up names
  sample_names <- gsub("\\.(AVE|STD)$", "", sample_cols)
  sample_names <- gsub("_S[0-9]+$", "", sample_names)
  return(unique(sample_names))
}

# Extract data for one sample
get_sample_data <- function(df, sample_idx) {
  col_start <- 3 + (sample_idx - 1) * 3
  data.frame(
    idx = 1:nrow(df),
    chr = df[, 1],
    pos = df[, 2],
    depth = as.numeric(df[, col_start]),
    ave = as.numeric(df[, col_start + 1]),
    std = as.numeric(df[, col_start + 2]),
    stringsAsFactors = FALSE
  )
}

# Draw alternating gene shading (cyan/white)
draw_gene_shading <- function(gene_bounds, x_range, y_min, y_max) {
  # Filter genes to visible range
  visible <- gene_bounds[gene_bounds$end >= x_range[1] & gene_bounds$start <= x_range[2], ]
  if (nrow(visible) == 0) return()

  visible$start <- pmax(visible$start, x_range[1])
  visible$end <- pmin(visible$end, x_range[2])

  # Draw alternating rectangles (cyan/white)
  for (i in 1:nrow(visible)) {
    col <- if (i %% 2 == 1) COL_GENE_SHADE else COL_GENE_SHADE_ALT
    rect(visible$start[i], y_min, visible$end[i], y_max, col = col, border = NA)
  }

  # Draw gene labels at top
  label_y <- 10^(log10(y_max) * 0.95)
  for (i in 1:nrow(visible)) {
    gene_width <- visible$end[i] - visible$start[i]
    if (gene_width > 200) {  # Only label genes with enough space
      label_x <- (visible$start[i] + visible$end[i]) / 2
      text(label_x, label_y, visible$gene[i], cex = 0.5, srt = 45, adj = c(0, 0.5))
    }
  }
}

# Draw gene shading for subplots (no labels, linear scale)
draw_gene_shading_linear <- function(gene_bounds, x_range, y_min, y_max) {
  visible <- gene_bounds[gene_bounds$end >= x_range[1] & gene_bounds$start <= x_range[2], ]
  if (nrow(visible) == 0) return()

  visible$start <- pmax(visible$start, x_range[1])
  visible$end <- pmin(visible$end, x_range[2])

  for (i in 1:nrow(visible)) {
    col <- if (i %% 2 == 1) COL_GENE_SHADE else COL_GENE_SHADE_ALT
    rect(visible$start[i], y_min, visible$end[i], y_max, col = col, border = NA)
  }
}

# Create coverage plot with subplots for a single page
create_coverage_plot <- function(data, gene_bounds, title, x_range = NULL) {

  # Calculate derived values
  data$ucl <- data$ave + data$std
  data$lcl <- pmax(1, data$ave - data$std)  # Min 1 for log scale
  data$zscore <- ifelse(data$std > 0, (data$depth - data$ave) / data$std, 0)
  data$ratio <- ifelse(data$ave > 0, data$depth / data$ave, 1)

  # Ensure depth is at least 1 for log scale
  data$depth <- pmax(1, data$depth)
  data$ave <- pmax(1, data$ave)
  data$ucl <- pmax(1, data$ucl)
  data$lcl <- pmax(1, data$lcl)

  # Filter to x_range if specified
  if (!is.null(x_range)) {
    data <- data[data$idx >= x_range[1] & data$idx <= x_range[2], ]
  } else {
    x_range <- c(min(data$idx), max(data$idx))
  }

  if (nrow(data) == 0) return()

  # Calculate y limits for coverage (log scale)
  y_min <- 1
  y_max <- max(c(data$depth, data$ucl, 500), na.rm = TRUE) * 2
  y_max <- min(y_max, 100000)  # Cap at 100000

  # Set up layout: Coverage (large), Z-score (small), Ratio (small)
  layout(matrix(c(1, 2, 3), nrow = 3), heights = c(3, 1, 1))

  # =========================================================================
  # COVERAGE PLOT (Top - main plot with log scale)
  # =========================================================================
  par(mar = c(0.5, 4, 3, 1))  # Minimal bottom margin, no x-axis labels

  # Create empty plot with log scale
  plot(NULL, xlim = x_range, ylim = c(y_min, y_max),
       xlab = "", ylab = "Coverage Depth",
       main = title, type = "n", xaxs = "i", yaxs = "i",
       log = "y", xaxt = "n")  # No x-axis

  # Draw gene shading (background)
  draw_gene_shading(gene_bounds, x_range, y_min, y_max)

  # Draw UCL/LCL polygon (light yellow confidence band)
  polygon(c(data$idx, rev(data$idx)),
          c(data$lcl, rev(data$ucl)),
          col = COL_UCL_LCL, border = NA)

  # Draw reference line at 250x
  abline(h = 250, col = COL_REF_LINE, lty = 2, lwd = 1)

  # Draw average line
  lines(data$idx, data$ave, col = COL_AVERAGE, lwd = 0.5)

  # Draw sample depth line (bold dark blue)
  lines(data$idx, data$depth, col = COL_COVERAGE, lwd = 1.5)

  # Add legend
  legend("topright",
         legend = c("Coverage", "Average", "UCL/LCL", "250x"),
         col = c(COL_COVERAGE, COL_AVERAGE, COL_UCL_LCL, COL_REF_LINE),
         lty = c(1, 1, NA, 2),
         lwd = c(1.5, 0.5, NA, 1),
         pch = c(NA, NA, 15, NA),
         pt.cex = c(NA, NA, 2, NA),
         cex = 0.6, bg = "white")

  box()

  # =========================================================================
  # Z-SCORE SUBPLOT (Middle - pink)
  # =========================================================================
  par(mar = c(0.5, 4, 0.5, 1))  # Minimal margins

  # Calculate Z-score limits (symmetric around 0)
  z_max <- max(abs(data$zscore), na.rm = TRUE)
  z_max <- min(z_max, 10)  # Cap at +/- 10
  z_max <- max(z_max, 3)   # At least +/- 3

  plot(NULL, xlim = x_range, ylim = c(-z_max, z_max),
       xlab = "", ylab = "",
       type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")

  # Draw gene shading
  draw_gene_shading_linear(gene_bounds, x_range, -z_max, z_max)

  # Reference line at 0
  abline(h = 0, col = "gray50", lty = 1, lwd = 0.5)

  # Draw Z-score line (pink)
  lines(data$idx, data$zscore, col = COL_ZSCORE, lwd = 1)

  # Label
  mtext("Z-score", side = 2, line = 1, cex = 0.7, col = COL_ZSCORE)

  box()

  # =========================================================================
  # RATIO SUBPLOT (Bottom - purple)
  # =========================================================================
  par(mar = c(3, 4, 0.5, 1))  # Bottom margin for x-axis

  # Calculate Ratio limits (centered around 1)
  r_max <- max(abs(data$ratio - 1), na.rm = TRUE) + 1
  r_max <- min(r_max, 3)  # Cap at 0-3
  r_max <- max(r_max, 2)  # At least 0-2

  plot(NULL, xlim = x_range, ylim = c(0, r_max),
       xlab = "Position Index", ylab = "",
       type = "n", xaxs = "i", yaxs = "i", yaxt = "n")

  # Draw gene shading
  draw_gene_shading_linear(gene_bounds, x_range, 0, r_max)

  # Reference line at 1
  abline(h = 1, col = "gray50", lty = 1, lwd = 0.5)

  # Draw Ratio line (purple)
  lines(data$idx, data$ratio, col = COL_RATIO, lwd = 1)

  # Label
  mtext("Ratio", side = 2, line = 1, cex = 0.7, col = COL_RATIO)

  box()

  # Reset layout
  layout(1)
}

# Generate PDF for one sample
generate_sample_pdf <- function(tmsp_df, cebpa_df, sample_idx, sample_name, output_dir) {
  pdf_file <- file.path(output_dir, sprintf("%s_Coverage.pdf", sample_name))

  cat(sprintf("  Generating %s...\n", basename(pdf_file)))

  # Get sample data
  tmsp_data <- get_sample_data(tmsp_df, sample_idx)

  # Create TMSP plots (5 pages)
  n_pos <- nrow(tmsp_data)
  page_size <- ceiling(n_pos / 5)

  pdf(pdf_file, width = 11, height = 8.5, paper = "USr")

  for (i in 1:5) {
    x_start <- (i - 1) * page_size + 1
    x_end <- min(i * page_size, n_pos)
    title <- sprintf("%s - TMSP Coverage %d/5", sample_name, i)

    create_coverage_plot(tmsp_data, tmsp_genes, title, c(x_start, x_end))
  }

  # Create CEBPA plot if data exists
  if (!is.null(cebpa_df) && ncol(cebpa_df) >= (2 + sample_idx * 3)) {
    cebpa_data <- get_sample_data(cebpa_df, sample_idx)
    if (!all(is.na(cebpa_data$depth))) {
      title <- sprintf("%s - CEBPA Coverage", sample_name)
      create_coverage_plot(cebpa_data, cebpa_genes, title)
    }
  }

  dev.off()

  cat(sprintf("  Created: %s\n", pdf_file))
  return(pdf_file)
}

# =============================================================================
# Main function
# =============================================================================
main <- function(args) {
  if (length(args) < 3) {
    cat("Usage: Rscript plotCoverage.R tmsp.csv cebpa.csv output_dir\n")
    cat("\nGenerates coverage plot PDFs from CSV files\n")
    cat("Output: One PDF per sample with TMSP (5 pages) and CEBPA coverage\n")
    cat("\nEach plot includes:\n")
    cat("  - Coverage (log scale, dark blue) with UCL/LCL (yellow)\n")
    cat("  - Z-score subplot (pink)\n")
    cat("  - Ratio subplot (purple)\n")
    quit(status = 1)
  }

  tmsp_file <- args[1]
  cebpa_file <- args[2]
  output_dir <- args[3]

  if (!file.exists(tmsp_file)) {
    stop(sprintf("File not found: %s", tmsp_file))
  }

  cat(sprintf("Reading TMSP data from: %s\n", tmsp_file))
  tmsp_df <- read.csv(tmsp_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  TMSP: %d positions, %d columns\n", nrow(tmsp_df), ncol(tmsp_df)))

  # Read CEBPA data
  cebpa_df <- NULL
  if (file.exists(cebpa_file) && file.info(cebpa_file)$size > 0) {
    cat(sprintf("Reading CEBPA data from: %s\n", cebpa_file))
    cebpa_df <- read.csv(cebpa_file, stringsAsFactors = FALSE, check.names = FALSE)
    cat(sprintf("  CEBPA: %d positions\n", nrow(cebpa_df)))
  } else {
    cat("  CEBPA: No data file found\n")
  }

  # Get sample names
  sample_names <- get_sample_names(tmsp_df)
  n_samples <- length(sample_names)

  cat(sprintf("\nFound %d samples: %s\n", n_samples, paste(sample_names, collapse = ", ")))
  cat(sprintf("Output directory: %s\n\n", output_dir))

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate PDF for each sample
  pdf_files <- c()
  for (i in 1:n_samples) {
    cat(sprintf("Processing sample %d/%d: %s\n", i, n_samples, sample_names[i]))
    pdf_file <- generate_sample_pdf(tmsp_df, cebpa_df, i, sample_names[i], output_dir)
    pdf_files <- c(pdf_files, pdf_file)
  }

  cat(sprintf("\nGenerated %d PDF files\n", length(pdf_files)))
  return(invisible(pdf_files))
}

# Run if called as script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  main(args)
}
