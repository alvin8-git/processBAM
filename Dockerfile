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
LABEL version="1.4"

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
    parallel \
    gawk \
    rename \
    procps \
    gcc \
    libc6-dev \
    locales \
    && rm -rf /var/lib/apt/lists/* \
    && sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen \
    && locale-gen

# Set locale environment
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# =============================================================================
# Install bioinformatics tools via conda/mamba
# =============================================================================
# Install tools in stages to avoid version conflicts
RUN mamba install -y -c bioconda -c conda-forge \
    samtools \
    mosdepth \
    pindel \
    python \
    openpyxl \
    r-base \
    && mamba clean -afy

# Install Java tools with compatible openjdk
RUN mamba install -y -c conda-forge \
    openjdk=17 \
    && mamba clean -afy

# Install picard (uses openjdk from above)
RUN mamba install -y -c bioconda -c conda-forge \
    picard \
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
# Install transpose utility
# =============================================================================
COPY transpose.c /tmp/transpose.c
RUN gcc -O2 -o /usr/local/bin/transpose /tmp/transpose.c && rm /tmp/transpose.c

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
