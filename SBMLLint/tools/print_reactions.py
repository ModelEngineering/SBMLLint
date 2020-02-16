#!/usr/bin/env python
"""Prints the reactions in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util

import argparse
import sys
import libsbml


def prettyPrint(model_reference, file_out=sys.stdout, **kwargs):
  """
  Prints the reactions in a model.
  :param str model_reference: file, xml string, antimony string
  :param dict kwargs: arguments to Reaction.getId
  """
  xml = util.getXML(model_reference)
  reader = libsbml.SBMLReader()
  document = reader.readSBMLFromString(xml)
  util.checkSBMLDocument(document)
  model = document.getModel()
  simple = SimpleSBML()
  simple.initialize(model)
  stgs = []
  for reaction in simple.reactions:
    stg = reaction.getId(**kwargs)
    stgs.append(stg)
    file_out.write("%s\n" % stg)
  return stgs

def main():
  def str2Bool(stg):
    if "T" in stg.upper():
      return True
    if "F" in stg.upper():
      return False
    return False
  #
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('filename', type=open, help='SBML file')
  parser.add_argument('--kinetics', nargs=1, type=str2Bool,
      help="Print kinetics formula True or False",
      default = ['True'])
  args = parser.parse_args()
  for fid in util.getNextFid(args.filename):
    try:
      prettyPrint(fid, is_include_kinetics=args.kinetics[0])
    except ValueError:
      print ("  *** Bad SBML file.")


if __name__ == '__main__':
  main()
