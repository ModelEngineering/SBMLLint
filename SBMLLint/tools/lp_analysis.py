#!/usr/bin/env python
"""
Performs Linear Programming Analysis to detection stoichiometric
inconsistencies.
"""

from SBMLLint.common import constants as cn
from SBMLLint.common import simple_sbml
from SBMLLint.common import stoichiometry_matrix

import argparse

import sys


def LPAnalysis(simple, is_report=False):
  """
  Does LP analysis for a simple model.
  :param bool: True if model is stoichiometric consistent.
  """
  sm_matrix = stoichiometry_matrix.StoichiometryMatrix(
      simple=simple)
  return sm_matrix.isConsistent(is_report_warning=is_report)

def main():
  def str2Bool(stg):
    if "T" in stg.upper():
      return True
    if "F" in stg.upper():
      return False
    return False
  #
  parser = argparse.ArgumentParser(
      description='LP Analysis of SBML XML file.')
  parser.add_argument('xml_fid', type=open, help='SBML file')
  parser.add_argument('--report_warnings', nargs=1,
      type=str2Bool,
      help="Print warnings if ill-formed matrix True or False",
      default = ['True'])
  args = parser.parse_args()
  simple = simple_sbml.SimpleSBML()
  simple.initialize(args.xml_fid)
  is_consistent = LPAnalysis(simple, 
      is_report=args.report_warnings[0])
  if is_consistent:
    print("Model is consistent.")
  else:
    print("Model is NOT consistent!")


if __name__ == '__main__':
  main()
