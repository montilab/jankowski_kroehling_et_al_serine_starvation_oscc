#!/bin/bash

# Usage: ./calculate_frip.sh <bam_file> <peak_file> <min_overlap>

# Input parameters
BAM_FILE=$1
PEAK_FILE=$2
MIN_OVERLAP=$3 #we used 0.1

# Count total fragments
TOTAL_FRAGMENTS=$(samtools view -c -F 260 ${BAM_FILE})

# Count fragments overlapping peaks
FRAGMENTS_IN_PEAKS=$(bedtools intersect -a ${BAM_FILE} -b ${PEAK_FILE} -bed -f ${MIN_OVERLAP} -c | awk '$NF>0' | wc -l)

# Calculate FRiP score
FRIP=$(echo "scale=4; ${FRAGMENTS_IN_PEAKS} / ${TOTAL_FRAGMENTS}" | bc)

# Output results
echo "Total Fragments: ${TOTAL_FRAGMENTS}"
echo "Fragments in Peaks: ${FRAGMENTS_IN_PEAKS}"
echo "FRiP Score: ${FRIP}"
