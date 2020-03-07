#!/bin/bash
# Script that installs SBMLLint in a virtual environment.
DIR="test-sbmllint"
echo ""
echo "**Will delete $DIR. Press any key to continue."
read -p "$*"
if [ -d "$DIR" ]; then
    rm -rf $DIR
fi
python3 -m venv $DIR
source $DIR/bin/activate
pip install Cython numpy
python3 setup.py install
echo "Success."
echo "Do: 'source $DIR/bin/activate' before running installed codes."
