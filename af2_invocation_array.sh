#!/bin/bash
#SBATCH --partition=gpu                      # Partition to run in
#SBATCH --gres=gpu:1                         # GPU resources requested
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=01:00:00                      # Runtime in D-HH:MM format
#SBATCH --mem=128GB                          # Requested Memory
#SBATCH --output=output-%x-%A_%a.log         # Changed to include array task ID
#SBATCH --error=error-%x-%A_%a.log           # Changed to include array task ID
#SBATCH --array=1-3200                       # Array job based on number of lines in data file

echo $HOSTNAME

module purge
module load colabfold

export OPENMM_CUDA_COMPILER=$(which nvcc)
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"

# Read the specific line from the data file based on array task ID
# Assuming each line in af2_fix_runs.dat is formatted as: input_path;output_path
#line=$(sed "${SLURM_ARRAY_TASK_ID}q;d" af2_fix_runs.dat)

# Split the line into input and output paths
#MY_PROTEIN_PATH=$(echo $line | cut -d';' -f1)
#MY_OUTPUT_DIR=$(echo $line | cut -d';' -f2)

#echo "Processing: $MY_PROTEIN_PATH -> $MY_OUTPUT_DIR"
# Run colabfold_search to generate MSAs and structural templates
#colabfold_search \
#    --mmseqs mmseqs \
#    --fasta-path $MY_PROTEIN_PATH \
#    --output-path $MY_OUTPUT_DIR \
#    --multimer \
#    --custom-template $custom_template_path \

# Run colabfold_batch to generate protein structure predictions
colabfold_batch \
    --num-recycle 5 \
    --num-seeds 1 \
    --random-seed 0 \
    --num-models 5 \
    --rank "ptm" \
    --use-gpu-relax \
    /nfs/scistore20/bronsgrp/jsaleh/data/missing_seq4.fasta \
    /nfs/scistore20/bronsgrp/jsaleh/data/missing_structures

