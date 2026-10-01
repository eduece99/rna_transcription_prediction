#!/bin/bash

conda activate venv_scb

python -V
echo "Using Anaconda environment: $CONDA_PREFIX"

echo "Running 1_create_models.py"
python src/1_create_models.py

echo "Running 2_create_predictions.py"
python src/2_create_predictions.py

echo "Running 3_create_constructs.py"
python src/3_create_constructs.py

echo "Finished running all scripts."
