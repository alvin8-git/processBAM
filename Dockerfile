# =============================================================================
# Dockerfile for processBAM Pipeline
# =============================================================================
# Lightweight container with all bioinformatics tools via Mambaforge
#
# Build:
#   docker build -t processbam:latest .
#
# Run:
#   docker run -v /path/to/Databases:/databases \
#              -v /path/to/analysis:/data \
#              processbam:latest /data/bam
#
# The container expects:
#   - /databases volume mounted with reference files
#   - /data volume mounted with analysis directory (containing bam/, vcf/, etc.)
# =============================================================================

FROM condaforge/mambaforge:latest

LABEL maintainer="Alvin Ng"
LABEL description="processBAM pipeline for TMSP/CEBPA BAM processing"
LABEL version="1.5"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# =============================================================================
# Install system dependencies
# =============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    ca-certificates \
    mawk \
    rename \
    procps \
    gcc \
    make \
    libc6-dev \
    locales \
    && rm -rf /var/lib/apt/lists/* \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen \
    && update-alternatives --set awk /usr/bin/mawk

# Install GNU Parallel (newer version with --memsuspend support)
RUN wget -q https://ftpmirror.gnu.org/parallel/parallel-20211222.tar.bz2 -O /tmp/parallel.tar.bz2 && \
    cd /tmp && tar -xjf parallel.tar.bz2 && \
    cd parallel-20211222 && ./configure && make && make install && \
    cd / && rm -rf /tmp/parallel*

# Set locale environment
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# =============================================================================
# Set pipeline environment variables (for both entrypoint and interactive use)
# =============================================================================
ENV DOCKER_MODE=true
ENV SCRIPT_DIR=/scripts
ENV HG19=/databases/WholeGenomeFASTA/genome.fa
ENV BED_TMSP=/databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.NoZero.v2.bed
ENV BED_TMSP_REFSEQ=/databases/TMSPvcf/BEDfiles/TMSP.RefSeqexons.bed
ENV BED_CEBNX=/databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.CEBPA.bed

# =============================================================================
# Install bioinformatics tools via conda/mamba
# =============================================================================
# Install specific versions matching local conda environment
RUN mamba install -y -c bioconda -c conda-forge \
    samtools=1.14 \
    mosdepth=0.3.3 \
    pindel \
    python \
    openpyxl \
    r-base \
    && mamba clean -afy

# Install Java 8 (required for Picard 2.x and GATK3)
RUN mamba install -y -c conda-forge \
    openjdk=8 \
    && mamba clean -afy

# Install Picard 2.26.10 (matching local version - Picard 3.x has different output format)
RUN mamba install -y -c bioconda -c conda-forge \
    picard=2.26.10 \
    && mamba clean -afy

# Download and setup GATK3 manually (required for VariantEval)
RUN mkdir -p /opt/gatk3 && \
    wget -q https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
    -O /tmp/gatk3.tar.bz2 && \
    tar -xjf /tmp/gatk3.tar.bz2 -C /opt/gatk3 --strip-components=1 --no-same-owner && \
    rm /tmp/gatk3.tar.bz2 && \
    echo '#!/bin/bash' > /usr/local/bin/gatk3 && \
    echo 'java -jar /opt/gatk3/GenomeAnalysisTK.jar "$@"' >> /usr/local/bin/gatk3 && \
    chmod +x /usr/local/bin/gatk3

# =============================================================================
# Install transpose utility (use pre-compiled binary from local system)
# =============================================================================
# Note: transpose.bin is the "Version 2.0 - Dr. Alex Sheppard" binary that
# correctly handles empty lines in Picard output files
COPY transpose.bin /usr/local/bin/transpose
RUN chmod +x /usr/local/bin/transpose

# =============================================================================
# Create directories for mounted volumes
# =============================================================================
RUN mkdir -p /databases /data /scripts

# =============================================================================
# Copy pipeline scripts
# =============================================================================
COPY processBAM.sh /scripts/
COPY generateCoveragePlots.sh /scripts/
COPY generateQCreport.sh /scripts/
COPY plotCoverage.R /scripts/
COPY plotQC.R /scripts/
RUN chmod +x /scripts/*.sh

# =============================================================================
# Create wrapper script for Docker environment
# =============================================================================
COPY processbam-wrapper.sh /usr/local/bin/processbam
RUN chmod +x /usr/local/bin/processbam

# Set working directory
WORKDIR /data

# Default command
ENTRYPOINT ["/usr/local/bin/processbam"]
CMD ["--help"]
