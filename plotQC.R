#!/usr/bin/env Rscript
# plotQC.R - Generate QC report PDFs from QCdata.xlsx
# Creates a 2-page PDF per sample with Picard and Mosdepth statistics
#
# Usage: Rscript plotQC.R picard.csv mosdepth.csv output_dir
#
# Input:
#   - picard.csv: Picard metrics (rows: metrics, cols: samples)
#   - mosdepth.csv: Mosdepth coverage per exon (rows: exons, cols: samples)
#
# Output: One PDF per sample with:
#   - Page 1: Sequencing, Alignment, and Variant statistics
#   - Page 2: Gene coverage table with color-coded proportions

# =============================================================================
# Gene definitions for aggregation (54 genes in TMSP panel)
# =============================================================================
TMSP_GENES <- c(
  "ABL1", "ASXL1", "ATRX", "BCOR", "BCORL1", "BRAF", "CALR", "CBL",
  "CBLB", "CBLC", "CDKN2A", "CEBPA", "CSF3R", "CUX1", "DNMT3A", "ETV6",
  "EZH2", "FBXW7", "FLT3", "GATA1", "GATA2", "GNAS", "HRAS", "IDH1",
  "IDH2", "IKZF1", "JAK2", "JAK3", "KDM6A", "KIT", "KMT2A", "KRAS",
  "MPL", "MYD88", "NOTCH1", "NPM1", "NRAS", "PDGFRA", "PHF6", "PTEN",
  "PTPN11", "RAD21", "RUNX1", "SETBP1", "SF3B1", "SMC1A", "SMC3", "SRSF2",
  "STAG2", "TET2", "TP53", "U2AF1", "WT1", "ZRSR2"
)

# Exons targeted per gene (from TMSP panel design)
GENE_EXONS <- list(
  ABL1 = "4-6", ASXL1 = "12", ATRX = "8,10-11,17-31", BCOR = "2-15",
  BCORL1 = "2-13", BRAF = "15", CALR = "9", CBL = "8-9",
  CBLB = "9-10", CBLC = "9,11", CDKN2A = "1-3", CEBPA = "1",
  CSF3R = "14-17", CUX1 = "1-24", DNMT3A = "1-22", ETV6 = "2-8",
  EZH2 = "2-20", FBXW7 = "9-11", FLT3 = "13-15,20", GATA1 = "2",
  GATA2 = "2-6", GNAS = "8-10", HRAS = "2-3", IDH1 = "4",
  IDH2 = "4", IKZF1 = "2-8", JAK2 = "12,14", JAK3 = "13",
  KDM6A = "1-30", KIT = "2,8-11,13,17", KMT2A = "1,3,5-8,27", KRAS = "2-3",
  MPL = "10", MYD88 = "3-5", NOTCH1 = "26-28,34", NPM1 = "11",
  NRAS = "2-3", PDGFRA = "12,14,18", PHF6 = "2-10", PTEN = "5,7",
  PTPN11 = "3,13", RAD21 = "2-13", RUNX1 = "2-9", SETBP1 = "4",
  SF3B1 = "13-16", SMC1A = "2,11,16-17", SMC3 = "10,13,19,23,25,28", SRSF2 = "1",
  STAG2 = "3-35", TET2 = "3-11", TP53 = "2-11", U2AF1 = "2,6",
  WT1 = "7,9", ZRSR2 = "1-10"
)

# =============================================================================
# Color functions for coverage proportion
# =============================================================================
# Color scale: dark red (0) -> medium red (0.5) -> white (1)
# 0.0 = dark red (bad coverage)
# 0.5 = medium red
# 1.0 = white (good coverage)
coverage_color <- function(prop) {
  prop <- max(0, min(1, prop))  # Clamp to 0-1 (single value)
  # Use a gradient from dark red (139,0,0) to salmon (248,105,107) to white (255,255,255)
  if (prop <= 0.5) {
    # 0 to 0.5: dark red to salmon
    t <- prop * 2  # Scale 0-0.5 to 0-1
    r <- 139 + (248 - 139) * t
    g <- 0 + (105 - 0) * t
    b <- 0 + (107 - 0) * t
  } else {
    # 0.5 to 1: salmon to white
    t <- (prop - 0.5) * 2  # Scale 0.5-1 to 0-1
    r <- 248 + (255 - 248) * t
    g <- 105 + (255 - 105) * t
    b <- 107 + (255 - 107) * t
  }
  rgb(r/255, g/255, b/255)
}

# =============================================================================
# Format large numbers with commas
# =============================================================================
format_number <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE)
}

# Format percentage with brackets
format_pct <- function(x) {
  sprintf("(%.1f%%)", x * 100)
}

# =============================================================================
# Parse Picard data
# =============================================================================
parse_picard <- function(picard_df, sample_idx) {
  # Column 1 is metric names, columns 2+ are samples
  col <- sample_idx + 1

  # Create named list of metrics
  metrics <- list()
  for (i in 1:nrow(picard_df)) {
    name <- as.character(picard_df[i, 1])
    value <- picard_df[i, col]
    metrics[[name]] <- value
  }
  metrics
}

# =============================================================================
# Aggregate mosdepth data by gene
# =============================================================================
aggregate_mosdepth <- function(mosdepth_df, sample_idx) {
  # Columns: Chr, Start, End, Gene, Transcript, Exon, then for each sample:
  # Coverage, 1X, 100X, 250X, 1000X (5 columns per sample)
  #
  # Proportions are calculated as per TSMPQCtemplate:
  # - T (Actual Bases) = 1X value (bases with at least 1X coverage)
  # - AA (1X prop) = 1X / 1X = 1.0 (always)
  # - AB (100X prop) = 100X / 1X
  # - AC (250X prop) = 250X / 1X
  # - AD (1000X prop) = 1000X / 1X

  col_start <- 7 + (sample_idx - 1) * 5

  # Extract sample data
  data <- data.frame(
    gene = mosdepth_df[, 4],
    exon_bases = mosdepth_df[, 3] - mosdepth_df[, 2],  # End - Start (for coverage calc)
    coverage = as.numeric(mosdepth_df[, col_start]),
    x1 = as.numeric(mosdepth_df[, col_start + 1]),
    x100 = as.numeric(mosdepth_df[, col_start + 2]),
    x250 = as.numeric(mosdepth_df[, col_start + 3]),
    x1000 = as.numeric(mosdepth_df[, col_start + 4]),
    stringsAsFactors = FALSE
  )

  # Aggregate by gene
  gene_summary <- data.frame(
    gene = character(),
    exons = character(),
    bases = numeric(),
    coverage = numeric(),
    x1 = numeric(),
    x100 = numeric(),
    x250 = numeric(),
    x1000 = numeric(),
    stringsAsFactors = FALSE
  )

  for (g in TMSP_GENES) {
    gene_data <- data[data$gene == g, ]
    if (nrow(gene_data) > 0) {
      total_exon_bases <- sum(gene_data$exon_bases)
      total_cov <- sum(gene_data$coverage * gene_data$exon_bases)
      total_x1 <- sum(gene_data$x1)      # Actual Bases (T)
      total_x100 <- sum(gene_data$x100)
      total_x250 <- sum(gene_data$x250)
      total_x1000 <- sum(gene_data$x1000)

      # Bases column = exon length (S column in TSMPQCtemplate)
      # Proportions use 1X (Actual Bases / T column) as denominator
      gene_summary <- rbind(gene_summary, data.frame(
        gene = g,
        exons = GENE_EXONS[[g]],
        bases = total_exon_bases,  # Exon length (S column), not 1X
        coverage = total_cov / total_exon_bases,  # Weighted average coverage
        x1 = if (total_x1 > 0) total_x1 / total_x1 else 0,      # Always 1.0
        x100 = if (total_x1 > 0) total_x100 / total_x1 else 0,
        x250 = if (total_x1 > 0) total_x250 / total_x1 else 0,
        x1000 = if (total_x1 > 0) total_x1000 / total_x1 else 0,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Add overall row
  # Per TSMPQCtemplate:
  # - Overall Coverage (V56) = U56/T56 = sum(coverage*bases) / sum(1X)
  # - Overall Bases (PRINT2) = T56 = sum(1X), NOT S56 (exon length)
  total_all_x1 <- sum(data$x1)  # T56 equivalent (Actual Bases)
  overall <- data.frame(
    gene = "Overall",
    exons = "",
    bases = total_all_x1,  # T56 = Total Actual Bases (1X sum)
    coverage = sum(data$coverage * data$exon_bases) / total_all_x1,  # U56/T56
    x1 = if (total_all_x1 > 0) total_all_x1 / total_all_x1 else 0,  # 1.0
    x100 = if (total_all_x1 > 0) sum(data$x100) / total_all_x1 else 0,
    x250 = if (total_all_x1 > 0) sum(data$x250) / total_all_x1 else 0,
    x1000 = if (total_all_x1 > 0) sum(data$x1000) / total_all_x1 else 0,
    stringsAsFactors = FALSE
  )

  rbind(gene_summary, overall)
}

# =============================================================================
# Draw Page 1: Picard Statistics
# =============================================================================
draw_page1 <- function(metrics, sample_name, total_targeted_bases) {
  par(mar = c(1, 1, 2, 1))
  plot.new()

  # Title
  title(main = paste(sample_name, "- QC Report"), cex.main = 1.5, font.main = 2)

  # Calculate derived values
  total_reads <- as.numeric(metrics$TOTAL_READS)
  r1_reads <- as.numeric(metrics$TOTAL_READS_R1)
  r2_reads <- as.numeric(metrics$TOTAL_READS_R2)
  read_length <- as.numeric(metrics$READ_LENGTH)
  total_bases <- as.numeric(metrics$TOTAL_BASES)
  q20_bases <- as.numeric(metrics$Q20_BASES)
  q30_bases <- as.numeric(metrics$Q30_BASES)
  aligned_reads <- as.numeric(metrics$PF_READS_ALIGNED)
  aligned_r1 <- as.numeric(metrics$PF_READS_ALIGNED_R1)
  aligned_r2 <- as.numeric(metrics$PF_READS_ALIGNED_R2)
  aligned_bases <- as.numeric(metrics$PF_ALIGNED_BASES)
  aligned_bases_r1 <- as.numeric(metrics$PF_ALIGNED_BASES_R1)
  aligned_bases_r2 <- as.numeric(metrics$PF_ALIGNED_BASES_R2)
  variants <- as.numeric(metrics$Variants)
  snps <- as.numeric(metrics$SNPs)
  insertions <- as.numeric(metrics$Insertions)
  deletions <- as.numeric(metrics$Deletions)
  het <- as.numeric(metrics$Heterozygous)
  hom <- as.numeric(metrics$Homozygous)
  ti <- as.numeric(metrics$Transitions)
  tv <- as.numeric(metrics$Transversions)
  titv <- as.numeric(metrics[["Ti/Tv"]])

  # Calculate Mean Coverage = PF_ALIGNED_BASES / total_targeted_bases
  mean_coverage <- aligned_bases / total_targeted_bases

  # Text positions - values right-aligned closer to units
  y_start <- 0.96  # Moved higher to reduce space above title
  y_step <- 0.025    # Reduced from 0.028 to fit DNA substitutions section
  x_label <- 0.08
  x_value <- 0.48    # Moved right (was 0.35) - right edge of value
  x_unit <- 0.50     # Unit position (was 0.52)
  x_pct <- 0.65      # Percentage position (was 0.62)

  y <- y_start

  # Section: Sequencing Statistics
  text(0.05, y, "Sequencing Statistics", cex = 1.2, font = 2, adj = 0)
  y <- y - y_step * 1.5

  # Reads subsection
  text(x_label, y, "Reads", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Total", adj = 0)
  text(x_value, y, format_number(total_reads), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  y <- y - y_step
  text(x_label, y, "R1", adj = 0)
  text(x_value, y, format_number(r1_reads), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  text(x_pct, y, format_pct(r1_reads / total_reads), adj = 0)
  y <- y - y_step
  text(x_label, y, "R2", adj = 0)
  text(x_value, y, format_number(r2_reads), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  text(x_pct, y, format_pct(r2_reads / total_reads), adj = 0)
  y <- y - y_step
  text(x_label, y, "Read Length", adj = 0)
  text(x_value, y, format_number(read_length), adj = 1)
  text(x_unit, y, "bp", adj = 0)
  y <- y - y_step * 1.5

  # Bases subsection
  text(x_label, y, "Bases", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Total", adj = 0)
  text(x_value, y, format_number(total_bases), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  y <- y - y_step
  text(x_label, y, "Q20", adj = 0)
  text(x_value, y, format_number(q20_bases), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  text(x_pct, y, format_pct(q20_bases / total_bases), adj = 0)
  y <- y - y_step
  text(x_label, y, "Q30", adj = 0)
  text(x_value, y, format_number(q30_bases), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  text(x_pct, y, format_pct(q30_bases / total_bases), adj = 0)
  y <- y - y_step * 2

  # Section: Alignment Statistics
  text(0.05, y, "Alignment Statistics", cex = 1.2, font = 2, adj = 0)
  y <- y - y_step * 1.5

  # Mean coverage = PF_ALIGNED_BASES / total_targeted_bases
  text(x_label, y, "Mean Coverage", adj = 0)
  text(x_value, y, sprintf("%.0f", mean_coverage), adj = 1)
  text(x_unit, y, "X", adj = 0)
  y <- y - y_step * 1.5

  # Reads Aligned subsection
  text(x_label, y, "Reads Aligned", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Total", adj = 0)
  text(x_value, y, format_number(aligned_reads), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  text(x_pct, y, format_pct(aligned_reads / total_reads), adj = 0)
  y <- y - y_step
  text(x_label, y, "R1", adj = 0)
  text(x_value, y, format_number(aligned_r1), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  text(x_pct, y, format_pct(aligned_r1 / r1_reads), adj = 0)
  y <- y - y_step
  text(x_label, y, "R2", adj = 0)
  text(x_value, y, format_number(aligned_r2), adj = 1)
  text(x_unit, y, "reads", adj = 0)
  text(x_pct, y, format_pct(aligned_r2 / r2_reads), adj = 0)
  y <- y - y_step * 1.5

  # Bases Aligned subsection
  text(x_label, y, "Bases Aligned", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Total", adj = 0)
  text(x_value, y, format_number(aligned_bases), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  text(x_pct, y, format_pct(aligned_bases / total_bases), adj = 0)
  y <- y - y_step
  text(x_label, y, "R1", adj = 0)
  text(x_value, y, format_number(aligned_bases_r1), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  text(x_pct, y, format_pct(aligned_bases_r1 / (r1_reads * read_length)), adj = 0)
  y <- y - y_step
  text(x_label, y, "R2", adj = 0)
  text(x_value, y, format_number(aligned_bases_r2), adj = 1)
  text(x_unit, y, "bases", adj = 0)
  text(x_pct, y, format_pct(aligned_bases_r2 / (r2_reads * read_length)), adj = 0)
  y <- y - y_step * 2

  # Section: Variant Statistics
  text(0.05, y, "Variant Statistics", cex = 1.2, font = 2, adj = 0)
  y <- y - y_step * 1.5

  text(x_label, y, "Total Variants", adj = 0)
  text(x_value, y, format_number(variants), adj = 1)
  y <- y - y_step
  text(x_label, y, "SNPs", adj = 0)
  text(x_value, y, format_number(snps), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(snps / variants), adj = 0)
  y <- y - y_step
  text(x_label, y, "Insertions", adj = 0)
  text(x_value, y, format_number(insertions), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(insertions / variants), adj = 0)
  y <- y - y_step
  text(x_label, y, "Deletions", adj = 0)
  text(x_value, y, format_number(deletions), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(deletions / variants), adj = 0)
  y <- y - y_step * 1.5

  # Zygosity
  text(x_label, y, "Zygosity", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Heterozygous", adj = 0)
  text(x_value, y, format_number(het), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(het / variants), adj = 0)
  y <- y - y_step
  text(x_label, y, "Homozygous", adj = 0)
  text(x_value, y, format_number(hom), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(hom / variants), adj = 0)
  y <- y - y_step
  text(x_label, y, "Het/Hom ratio", adj = 0)
  text(x_value, y, sprintf("%.2f", het / hom), adj = 1)
  y <- y - y_step * 1.5

  # DNA substitutions
  text(x_label, y, "DNA substitutions", cex = 1.0, font = 2, adj = 0)
  y <- y - y_step
  text(x_label, y, "Transitions", adj = 0)
  text(x_value, y, format_number(ti), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(ti / snps), adj = 0)
  y <- y - y_step
  text(x_label, y, "Transversions", adj = 0)
  text(x_value, y, format_number(tv), adj = 1)
  text(x_unit, y, "variants", adj = 0)
  text(x_pct, y, format_pct(tv / snps), adj = 0)
  y <- y - y_step
  text(x_label, y, "Ti/Tv ratio", adj = 0)
  text(x_value, y, sprintf("%.2f", titv), adj = 1)
}

# =============================================================================
# Draw Page 2: Mosdepth Gene Coverage Table
# =============================================================================
draw_page2 <- function(gene_summary, sample_name) {
  par(mar = c(1, 1, 2, 1))
  plot.new()

  # Title
  title(main = paste("Gene Coverage Statistics -", sample_name), cex.main = 1.2, font.main = 2)

  # Column headers
  cols <- c(0.05, 0.18, 0.30, 0.42, 0.54, 0.66, 0.78, 0.90)
  headers <- c("Gene", "Exons", "Coverage", "Bases", "1X", "100X", "250X", ">1000X")

  y_start <- 0.94
  y_header <- y_start
  y_step <- 0.0155  # Smaller step for 55 rows + header

  # Draw header - "Proportion of Bases Covered @" centered above 1X, 100X, 250X, >1000X columns
  # Columns 5-8 are at 0.54, 0.66, 0.78, 0.90 - center is (0.54 + 0.90)/2 = 0.72
  text(0.72, y_header + 0.02, "Proportion of Bases Covered @", cex = 0.8, adj = 0.5)
  for (i in 1:length(headers)) {
    text(cols[i], y_header, headers[i], cex = 0.7, font = 2, adj = 0)
  }

  # Draw horizontal line under header
  segments(0.03, y_header - 0.008, 0.97, y_header - 0.008, lwd = 0.5)

  y <- y_header - y_step

  # Draw data rows
  for (row in 1:nrow(gene_summary)) {
    g <- gene_summary[row, ]

    # Gene name (bold for Overall)
    font_weight <- if (g$gene == "Overall") 2 else 1
    text(cols[1], y, g$gene, cex = 0.6, font = font_weight, adj = 0)

    # Exons
    text(cols[2], y, g$exons, cex = 0.6, adj = 0)

    # Coverage (formatted)
    text(cols[3], y, sprintf("%.0f", g$coverage), cex = 0.6, adj = 0)

    # Bases
    text(cols[4], y, format_number(g$bases), cex = 0.6, adj = 0)

    # Coverage proportions with colored background
    prop_cols <- c(5, 6, 7, 8)
    prop_vals <- c(g$x1, g$x100, g$x250, g$x1000)

    for (j in 1:4) {
      # Draw colored rectangle
      rect(cols[prop_cols[j]] - 0.01, y - 0.006,
           cols[prop_cols[j]] + 0.09, y + 0.01,
           col = coverage_color(prop_vals[j]), border = NA)

      # Draw proportion as decimal (0.00 to 1.00) with 2 decimal places
      text(cols[prop_cols[j]], y, sprintf("%.2f", prop_vals[j]), cex = 0.6, adj = 0)
    }

    # Draw horizontal line after ZRSR2 (last gene before Overall)
    if (g$gene == "ZRSR2") {
      segments(0.03, y - 0.008, 0.97, y - 0.008, lwd = 0.5)
    }

    y <- y - y_step
  }
}

# =============================================================================
# Generate QC report PDF for one sample
# =============================================================================
generate_qc_pdf <- function(picard_df, mosdepth_df, sample_idx, sample_name, output_dir) {
  pdf_file <- file.path(output_dir, sprintf("%s_QC.pdf", sample_name))

  cat(sprintf("  Generating %s...\n", basename(pdf_file)))

  # Parse Picard metrics
  metrics <- parse_picard(picard_df, sample_idx)

  # Aggregate mosdepth by gene
  gene_summary <- aggregate_mosdepth(mosdepth_df, sample_idx)

  # Calculate total targeted bases (S56 = sum of exon lengths) for Mean Coverage
  # Note: gene_summary$bases for Overall row is T56 (1X sum), not S56
  # We need S56 for Mean Coverage calculation
  total_targeted_bases <- sum(mosdepth_df[, 3] - mosdepth_df[, 2])  # Sum of (End - Start)

  # Create PDF
  pdf(pdf_file, width = 8.5, height = 11, paper = "letter")

  # Page 1: Picard statistics
  draw_page1(metrics, sample_name, total_targeted_bases)

  # Page 2: Gene coverage table
  draw_page2(gene_summary, sample_name)

  dev.off()

  cat(sprintf("  Created: %s\n", pdf_file))
  return(pdf_file)
}

# =============================================================================
# Get sample names from Picard data
# =============================================================================
get_sample_names <- function(picard_df) {
  # Column 1 is metric names, rest are sample names
  colnames(picard_df)[-1]
}

# =============================================================================
# Main function
# =============================================================================
main <- function(args) {
  if (length(args) < 3) {
    cat("Usage: Rscript plotQC.R picard.csv mosdepth.csv output_dir\n")
    cat("\nGenerates QC report PDFs from CSV files\n")
    cat("Output: One PDF per sample with:\n")
    cat("  - Page 1: Sequencing, Alignment, and Variant statistics\n")
    cat("  - Page 2: Gene coverage table with color-coded proportions\n")
    quit(status = 1)
  }

  picard_file <- args[1]
  mosdepth_file <- args[2]
  output_dir <- args[3]

  if (!file.exists(picard_file)) {
    stop(sprintf("Picard file not found: %s", picard_file))
  }
  if (!file.exists(mosdepth_file)) {
    stop(sprintf("Mosdepth file not found: %s", mosdepth_file))
  }

  cat(sprintf("Reading Picard data from: %s\n", picard_file))
  picard_df <- read.csv(picard_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  Picard: %d metrics, %d samples\n", nrow(picard_df), ncol(picard_df) - 1))

  cat(sprintf("Reading Mosdepth data from: %s\n", mosdepth_file))
  mosdepth_df <- read.csv(mosdepth_file, stringsAsFactors = FALSE, check.names = FALSE)
  cat(sprintf("  Mosdepth: %d exons, %d columns\n", nrow(mosdepth_df), ncol(mosdepth_df)))

  # Get sample names
  sample_names <- get_sample_names(picard_df)
  # Clean up sample names (remove _S[0-9]+ suffix)
  sample_names_clean <- gsub("_S[0-9]+$", "", sample_names)
  n_samples <- length(sample_names)

  cat(sprintf("\nFound %d samples: %s\n", n_samples, paste(sample_names_clean, collapse = ", ")))
  cat(sprintf("Output directory: %s\n\n", output_dir))

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate PDF for each sample
  pdf_files <- c()
  for (i in 1:n_samples) {
    cat(sprintf("Processing sample %d/%d: %s\n", i, n_samples, sample_names_clean[i]))
    pdf_file <- generate_qc_pdf(picard_df, mosdepth_df, i, sample_names_clean[i], output_dir)
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
