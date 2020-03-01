#!/bin/bash
# Randomly selects a different setup file
COUNT=$(($RANDOM*100/32767))
echo $COUNT
if [ "$COUNT" -gt 50 ]
then
  echo "Running without Tellurium setup"
  python setup.py install
else
  echo "Running with Tellurium setup"
  python setup_tellurium.py install
fi
