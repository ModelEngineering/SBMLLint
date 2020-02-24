"""Script for running mass balance checking tools."""


from SBMLLint.common import config
from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util
from SBMLLint.games.games_pp import GAMES_PP
from SBMLLint.games.games_report import GAMESReport
from SBMLLint.moiety_analysis.moiety_comparator import MoietyComparator

import os
import sys
import libsbml

TYPE_I = "type1"
TYPE_II = "type2"
TYPE_III = "type3"
CANCELING = "canceling"
ECHELON = "echelon"
GAMES = "games"


def lint(model_reference=None, 
    file_out=sys.stdout,
    mass_balance_check=GAMES,
    config_fid=None,
    is_report=True,
    implicit_games=False):
  """
  Reports on errors found in a model
  :param str model_reference: 
      libsbml_model file in
      file, antimony string, xml string
  :param TextIOWrapper model_fid: fid for an XML file
  :param TextIOWrapper file_out:
  :param str mass_balance_check: how check for mass balance
  :param TextIOWrapper config_fid: readable stream
  :param bool is_report: print result
  :return MoietyComparatorResult/null/None:
  """
  config.setConfiguration(fid=config_fid)
  config_dct = config.getConfiguration()
  if util.isSBMLModel(model_reference):
    model = model_reference
  else:
    xml = util.getXML(model_reference)
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromString(xml)
    util.checkSBMLDocument(document)
    model = document.getModel()
  #
  simple = SimpleSBML()
  simple.initialize(model)
  if mass_balance_check==cn.MOIETY_ANALYSIS:
    result = MoietyComparator.analyzeReactions(simple)
    if is_report:
      for line in result.report.split('\n'):
          file_out.write("%s\n" % line)
    return result
  elif mass_balance_check == GAMES:
    if implicit_games:
      for ignored in config_dct[cn.CFG_IGNORED_MOLECULES]:
        simple = removeIgnored(simple, ignored)
    m = GAMES_PP(simple)
    games_result = m.analyze(simple.reactions)
    if games_result and is_report:
      gr = GAMESReport(m, explain_threshold=config_dct[cn.CFG_GAMES_THRESHOLD])
      errortype_dic = {TYPE_I: gr.reportTypeOneError,
                       TYPE_II: gr.reportTypeTwoError,
                       TYPE_III: gr.reportTypeThreeError,
                       CANCELING: gr.reportCancelingError,
                       ECHELON: gr.reportEchelonError
                      }
      for errors in m.error_summary:
        for category in errortype_dic.keys():
          if errors.type == category:
            func = errortype_dic[category]            
            report, _ = func(errors.errors, explain_details=True)
            print(report)
    return games_result
  else:
    print ("Specified method doesn't exist")
    return None

def removeIgnored(simple, ignored):
  """
  Remove an ignored molecule
  from all reactions in a simpleSBML model.
  :param SimpleSBML simple:
  :param str ignored: a molecule name
  :return SimpleSBML:
  """
  modified_reactions = []
  for r in simple.reactions:
    modified = False
    reactant_names = [reactant.molecule.name for reactant in r.reactants]
    product_names = [product.molecule.name for product in r.products]
    if ignored in reactant_names:
      r.reactants = [ms for ms in r.reactants 
          if ms.molecule.name != ignored]
      modified = True
    if ignored in product_names:
      r.products = [ms for ms in r.products
           if ms.molecule.name != ignored]
      modified = True
    r.identifier = r.makeIdentifier()
    r.category = r.getCategory()
    if modified:
      modified_reactions.append(r)
  return simple
