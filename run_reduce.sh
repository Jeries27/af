#!/usr/bin/env bash
# File name: run_reduce.sh

MAIN_DIR="/nfs/scistore20/bronsgrp/jsaleh/data/organized_missing_structures"

shopt -s nullglob

for folder in "$MAIN_DIR"/*; do
  if [ -d "$folder" ]; then
    echo "Processing folder: $folder"

    # Look for PDB files matching the unrelaxed_rank_00 pattern
    for pdb_file in "$folder"/*_unrelaxed_rank_00*.pdb; do
      if [ -f "$pdb_file" ]; then
        out_file="${pdb_file%.pdb}_withH.pdb"

        # Check if the output file already exists
        if [ -f "$out_file" ]; then
          echo "  Skipping $pdb_file, output file $out_file already exists."
          continue
        fi

        echo "  Found PDB: $pdb_file"
        echo "  Running reduce on $pdb_file..."
        reduce -BUILD "$pdb_file" > "$out_file"
        echo "  Created: $out_file"
      fi
    done
  fi
done

echo "Done."
