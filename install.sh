#!/bin/bash
python3 -m venv test-sbmllint
source test-sbmllint/bin/activate
python3 setup.py install
echo "Success."
echo "Do: test-sbmllint/bin/activate before using installed codes."
