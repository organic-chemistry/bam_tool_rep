#!/bin/bash

# Usage: ./run_bam_filters.sh /path/to/folder

# folder containing BAM files
FOLDER="$1"

if [ -z "$FOLDER" ]; then
  echo "Usage: $0 /path/to/folder"
  exit 1
fi

# loop through all BAM files in the folder
for BAM in "$FOLDER"/*.bam; do
  if [ -f "$BAM" ]; then
    PICKLE="${BAM}.pickle"
    if [ -f "$PICKLE" ]; then
      echo "Skipping $BAM (pickle already exists: $PICKLE)"
      continue
    fi

    echo "Processing $BAM ..."
    python bam_filter_cli.py "$BAM" \
      --remove-shorter-than 5000 \
      --brdu-threshold 0.05 \
      --threads 3 \
      --check-size
  fi
done