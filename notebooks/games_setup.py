import collections
import itertools  
import matplotlib.pyplot as plt
import networkx as nx
import os
import re
import os, sys

import init
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games import print_model as pm
# from SBMLLint.games.som import SOM
# from SBMLLint.games.mesgraph import MESGraph


cwd = os.getcwd()
print("Current Directory:", cwd)

def load_file_from_games(num):

  format_num = format(num, '03d')
  file = os.path.join(os.getcwd(), os.pardir, 'data/biomodels/BIOMD0000000' + format_num + '_url.xml')
  simple = SimpleSBML()
  simple.initialize(file)
  return simple

def load_file_from_curated_data(num):

  format_num = format(num, '03d')
  file = os.path.join(os.getcwd(), os.pardir, os.pardir, 'curated_data/curated_' + format_num + '.xml')
  simple = SimpleSBML()
  simple.initialize(file)
  return simple
