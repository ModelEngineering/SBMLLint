"""
Performs Linear Programming Analysis to detection stoichiometric
inconsistencies.
"""

from SBMLLint.common import constants as cn
from SBMLLint.common import simple_sbml
from SBMLLint.common import stoichiometry_matrix

import argparse

def LPAnalysis(simple):
  """
  Does LP analysis for a simple model.
  :param bool: True if model is stoichiometric consistent.
  """
  sm_matrix = stoichiometry_matrix.StoichiometryMatrix(
      simple=simple)
  return sm_matrix.isConsistent()

def main():
  parser = argparse.ArgumentParser(
      description='LP Analysis of SBML XML file.')
  parser.add_argument('xml_fid', type=open, help='SBML file')
  args = parser.parse_args()
  simple = simple_sbml.SimpleSBML()
  simple.initialize(args.xml_fid)
  is_consistent = LPAnalysis(simple)
  if is_consistent:
    print "Model is consistent."
  else:
    print "Model is NOT consistent."
