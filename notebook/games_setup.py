import collections
import itertools  
import matplotlib.pyplot as plt
import networkx as nx
import os
import tellurium as te
import tesbml
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

def load_file(num):

  format_num = format(num, '03d')
  file = os.path.join(os.getcwd(), os.pardir, 'SBMLLint/games/data/curated_' + format_num + '.xml')
  simple = SimpleSBML()
  simple.initialize(file)
  return simple

