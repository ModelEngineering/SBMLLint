# Creates a conda environment in which to run the tests
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a

if [ ! -e /home/ubuntu/miniconda3/envs/test-environment ]; then
  echo "***Create environment***"
  conda create -q -n test-environment python=3.4 numpy pandas matplotli
  sudo $HOME/miniconda3/bin/conda install -c SBMLTeam python-libsbml
fi
conda info --envs
echo "***Done!."
echo "***Next: 'source activate test-environment'"
