"""Prints the reactions in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util

import argparse


def pprint(model_reference, **kwargs):
  """
  Prints the reactions in a model.
  :param str model_reference: file, xml string, antimony string
  :param dict kwargs: arguments to SimpleSBML.getReactionString
  """
  document = util.getSBMLDocument(model_reference)
  model = document.getModel()
  simple = SimpleSBML(model)
  for reaction in simple.reactions:
    stg = SimpleSBML.getReactionString(reaction, **kwargs)
    print("%s" % stg)

def main():
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('filename', type=str, help='SBML file')
  parser.add_argument('--kinetics', nargs=1,
      help="Print kinetics formula True or False",
      default = ['True'])
  args = parser.parse_args()
  is_include_kinetics = True if args.kinetics[0][0] == 'T' else False
  pprint(args.filename, is_include_kinetics=is_include_kinetics)


if __name__ == '__main__':
  main()
