#!/bin/bash

#SBATCH --partition=gpu                      # Partition to run in
#SBATCH --gres=gpu:1                         # GPU resources requested
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=1-00:00:00                        # Runtime in D-HH:MM format
#SBATCH --mem=128GB                          # Requested Memory
#SBATCH --output=output-%x-%j.log
#SBATCH --error=error-%x-%j.log

echo $HOSTNAME

module purge
module load colabfold

export OPENMM_CUDA_COMPILER=$(which nvcc)
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION="4.0"

# Set the path to the input protein
MY_PROTEIN_PATH="/nfs/scistore14/schangrp/pschanda/Colabfold/Nb21/Nb_CA18921_MP1914_G10.fasta"

# Replace this path with the path to your custom template file
#custom_template_path=""



# Set the output directory for AlphaFold
MY_OUTPUT_DIR="/nfs/scistore14/schangrp/pschanda/Colabfold/Nb21"

# Run colabfold_search to generate MSAs and structural templates
#colabfold_search \
#    --mmseqs mmseqs \
#    --fasta-path $MY_PROTEIN_PATH \
#    --output-path $MY_OUTPUT_DIR \
#    --multimer \
#    --custom-template $custom_template_path \

# Run colabfold_batch to generate protein structure predictions
colabfold_batch \
    --model-type "alphafold2_multimer_v3" \
    --num-recycle 5 \
    --recycle-early-stop-tolerance 0.01 \
    --num-ensemble 1 \
    --num-seeds 1 \
    --random-seed 0 \
    --num-models 5 \
    --recompile-padding 20 \
    --num-relax 5 \
    --rank "ptm" \
    --use-gpu-relax \
    $MY_PROTEIN_PATH \
    $MY_OUTPUT_DIR
#    --num-recycle 5 \
#    --recycle-early-stop-tolerance 0.01 \
#    --num-ensemble 1 \
#    --num-seeds 1 \
#    --random-seed 0 \
#    --num-models 5 \
#    --recompile-padding 20 \
#    --model-order "model_1,model_2,model_3,model_4,model_5" \
#    --data "casp14" \
#    --msa-mode "mmseqs2_uniref_env" \
#    --model-type "alphafold2_multimer_v3" \
#    --num-relax 3 \
#    --rank "ptm" \
#    --pair-mode "unpaired" \
#    --sort-queries-by "length" \
#    --zip \
#    --use-gpu-relax \
#    $MY_PROTEIN_PATH \
#    $MY_OUTPUT_DIR
#    --custom-template-path $custom_template_path \

