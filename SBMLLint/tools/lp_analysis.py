#!/usr/bin/env python
"""
Performs Linear Programming Analysis to detection stoichiometric
inconsistencies.
"""

from SBMLLint.common import constants as cn
from SBMLLint.common import simple_sbml
from SBMLLint.common import stoichiometry_matrix
from SBMLLint.common import util

import argparse

import sys


def LPAnalysis(fid, is_report=False):
  """
  Does LP analysis for a simple model.
  :param IOStream fid: XML file
  :param bool: True if model is stoichiometric consistent.
  """
  simple = simple_sbml.SimpleSBML()
  simple.initialize(fid)
  sm_matrix = stoichiometry_matrix.StoichiometryMatrix(
      simple=simple)
  is_consistent = sm_matrix.isConsistent(is_report_warning=is_report)
  if is_consistent:
    print("Model is consistent.")
  else:
    print("Model is NOT consistent!")
  return is_consistent

def main():
  def str2Bool(stg):
    if "T" in stg.upper():
      return True
    if "F" in stg.upper():
      return False
    return False
  #
  parser = argparse.ArgumentParser(
      description='LP Analysis of SBML file(s).')
  parser.add_argument('xml_fid', type=open,
      help='SBML file or zipfile')
  parser.add_argument('--report_warnings', nargs=1,
      type=str2Bool,
      help="Print warnings if ill-formed matrix True or False",
      default = ['True'])
  args = parser.parse_args()
  for fid in util.getNextFid(args.xml_fid):
    try:
      LPAnalysis(fid, is_report=args.report_warnings[0])
    except ValueError:
      print ("  *** Bad SBML file.")
    except Exception as e:
      print(e)


if __name__ == '__main__':
  main()
