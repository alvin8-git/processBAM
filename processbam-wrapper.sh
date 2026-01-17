#!/bin/bash
# Wrapper for processBAM in Docker environment
export HOME=/root
export DOCKER_MODE=true
export HG19=/databases/WholeGenomeFASTA/genome.fa
export BED_TMSP=/databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.NoZero.v2.bed
export BED_TMSP_REFSEQ=/databases/TMSPvcf/BEDfiles/TMSP.RefSeqexons.bed
export BED_CEBNX=/databases/TMSPvcf/BEDfiles/TSMP.UCSCexons.CEBPA.bed
export SCRIPT_DIR=/scripts

# Run from provided directory or current directory
if [ -n "$1" ] && [ -d "$1" ]; then
    cd "$1"
    shift
fi
/scripts/processBAM.sh "$@"
