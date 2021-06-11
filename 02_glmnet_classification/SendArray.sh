#!/bin/bash


# Set files to perform actions on

FILES=($(cat /home/mainciburu/scRNA/pipeline_2.0/02_glmnet_classification/celltypes.txt))
NUM=${#FILES[@]} # get size of array
ZBNUM=$(($NUM - 1)) 

# 1) Build individual models

if [ $ZBNUM -ge 0 ]; then
  sbatch --wait --array=0-$ZBNUM%5 /home/mainciburu/scRNA/pipeline_2.0/02_glmnet_classification/SlurmBatch_binary_models_R.sbs 
fi

wait

# 2) Final classification

sbatch --wait /home/mainciburu/scRNA/pipeline_2.0/02_glmnet_classification/SlurmBatch_final_classification_R.sbs
