#!/bin/bash

# create a conda environment from the env.yaml file
mamba env create -f ./env.yml

mamba activate venv_scb  # or failing that, conda activate venv

pip install pre-commit black ruff
pre-commit install
pre-commit run --all-files  # execute the pre-commit hooks on all files in the repository

