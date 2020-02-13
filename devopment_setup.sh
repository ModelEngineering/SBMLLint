# Setup for code development
conda update conda
conda env create -f environment.yml
conda activate sbmllint
pip install python-libsbml
pip install tellurium
