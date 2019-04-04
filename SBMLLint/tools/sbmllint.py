"""Checks for static errors in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util
from SBMLLint.games.mesgraph import MESGraph
from SBMLLint.structured_names.moiety_comparator import MoietyComparator

import argparse
import os
import sys
import tesbml


def lint(model_reference, file_out=sys.stdout,
    mass_balance_check="structured_names",
    is_report=True):
  """
  Reports on errors found in a model
  :param str model_reference: 
      libsbml_model, file, antimony string, xml string
  :param TextIOWrapper file_out:
  :param str mass_balance_check: how check for mass balance
  :param bool is_report: print result
  :return MoietyComparatorResult/bull/None:
  """
  if util.isSBMLModel(model_reference):
    model = model_reference
  else:
    xml = util.getXML(model_reference)
    reader = tesbml.SBMLReader()
    document = reader.readSBMLFromString(xml)
    util.checkSBMLDocument(document)
    model = document.getModel()
  simple = SimpleSBML()
  simple.initialize(model)
  if mass_balance_check == "structured_names":
    result = MoietyComparator.analyzeReactions(simple)
    if is_report:
      for line in result.report.split('\n'):
          file_out.write("%s\n" % line)
    return result
  elif mass_balance_check == "games":
    m = MESGraph(simple)
    report = m.analyze(simple.reactions)
    if is_report:
      print(report)
    return True
  else:
    print ("Specified method doesn't exist")
    return None
    

def main():
  parser = argparse.ArgumentParser(description='SBML XML file.')
  parser.add_argument('filename', type=str, help='SBML file')
  args = parser.parse_args()
  lint(args.filename)


if __name__ == '__main__':
  main()
