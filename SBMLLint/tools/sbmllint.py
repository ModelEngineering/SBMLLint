"""Checks for static errors in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util
from SBMLLint.structured_names.moiety_comparator import MoietyComparator

import argparse
import sys


def lint(model_reference, file_out=sys.stdout,
    mass_balance_check="structured_names"):
  """
  Reports on errors found in a model
  :param str model_reference: file, antimony string, xml string
  :param TextIOWrapper file_out:
  :param str mass_balance_check: how check for mass balance
  :return int, int: total reactions, number non-compliant
  """
  document = util.getSBMLDocument(model_reference)
  model = document.getModel()
  num_bad, report = MoietyComparator.analyzeReactions(model)
  num_reactions = len(Reaction.reactions)
  file_out.write("%d/%d reactions are impbalanced." 
      % (num_bad, num_reactions))
  for line in report.split('\n'):
    file_out.write("%s\n" % line)
  return num_reactions, num_bad
    

def main():
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('filename', type=str, help='SBML file')
  args = parser.parse_args()
  lint(args.filename)


if __name__ == '__main__':
  main()
