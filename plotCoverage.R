#!/usr/bin/env Rscript
# plotCoverage.R - Generate coverage plots from CSV files using base R graphics
# Creates multi-page PDF with alternating gene shading for TMSP and CEBPA panels
#
# Usage: Rscript plotCoverage.R tmsp.csv cebpa.csv output_dir
#
# Input: CSV files with columns: Chr, Pos, Sample1, Sample1.AVE, Sample1.STD, ...
# Output: One PDF per sample with TMSP (5 pages) and CEBPA (1 page) coverage plots

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

# Draw alternating gene shading
draw_gene_shading <- function(gene_bounds, x_range, y_max) {
  # Filter genes to visible range
  visible <- gene_bounds[gene_bounds$end >= x_range[1] & gene_bounds$start <= x_range[2], ]
  if (nrow(visible) == 0) return()

  visible$start <- pmax(visible$start, x_range[1])
  visible$end <- pmin(visible$end, x_range[2])

  # Draw alternating rectangles
  for (i in 1:nrow(visible)) {
    col <- if (i %% 2 == 1) rgb(0.9, 0.9, 0.9) else rgb(1, 1, 1)
    rect(visible$start[i], 0, visible$end[i], y_max, col = col, border = NA)
  }

  # Draw gene labels
  label_y <- y_max * 0.95
  for (i in 1:nrow(visible)) {
    gene_width <- visible$end[i] - visible$start[i]
    if (gene_width > 200) {  # Only label genes with enough space
      label_x <- (visible$start[i] + visible$end[i]) / 2
      text(label_x, label_y, visible$gene[i], cex = 0.5, srt = 45, adj = c(0, 0.5))
    }
  }
}

# Create coverage plot for a single page
create_coverage_plot <- function(data, gene_bounds, title, x_range = NULL) {

  # Add UCL/LCL
  data$ucl <- data$ave + data$std
  data$lcl <- pmax(0, data$ave - data$std)

  # Filter to x_range if specified
  if (!is.null(x_range)) {
    data <- data[data$idx >= x_range[1] & data$idx <= x_range[2], ]
  } else {
    x_range <- c(min(data$idx), max(data$idx))
  }

  if (nrow(data) == 0) return()

  # Calculate y limits
  y_max <- max(c(data$depth, data$ucl, 500), na.rm = TRUE) * 1.15
  y_max <- min(y_max, 10000)  # Cap at 10000

  # Set up plot area with margins for labels
  par(mar = c(4, 4, 3, 1))

  # Create empty plot
  plot(NULL, xlim = x_range, ylim = c(0, y_max),
       xlab = "Position Index", ylab = "Coverage Depth",
       main = title, type = "n", xaxs = "i", yaxs = "i")

  # Draw gene shading first (background)
  draw_gene_shading(gene_bounds, x_range, y_max)

  # Draw UCL/LCL polygon (confidence band)
  polygon(c(data$idx, rev(data$idx)),
          c(data$lcl, rev(data$ucl)),
          col = rgb(0.68, 0.85, 0.9, 0.3), border = NA)

  # Draw reference line at 250x
  abline(h = 250, col = "red", lty = 2, lwd = 1)

  # Draw average line
  lines(data$idx, data$ave, col = "blue", lwd = 0.5)

  # Draw sample depth line
  lines(data$idx, data$depth, col = "black", lwd = 0.7)

  # Add legend
  legend("topright",
         legend = c("Sample Depth", "Average", "UCL/LCL", "250x"),
         col = c("black", "blue", rgb(0.68, 0.85, 0.9), "red"),
         lty = c(1, 1, NA, 2),
         lwd = c(1, 0.5, NA, 1),
         pch = c(NA, NA, 15, NA),
         pt.cex = c(NA, NA, 2, NA),
         cex = 0.6, bg = "white")

  # Re-draw box
  box()
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
