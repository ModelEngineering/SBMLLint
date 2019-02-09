"""Checks for static errors in a model."""

from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.reaction import Reaction
from SBMLLint.common import util
from SBMLLint.structured_names.moiety_comparator import MoietyComparator

NUM_S1 = 2
NUM_S2 = 3
MOLECULE1 = "S1"
MOLECULE2 = "S2"
ANTIMONY_STG = '''
%d%s-> %d%s; 1
S1 = 0
S2 = 0
''' % (NUM_S1, MOLECULE1, NUM_S2, MOLECULE2)

def lint(model_reference,
    mass_balance_check="structured_names"):
  """
  Reports on errors found in a model
  :param str model_reference: file, antimony, xml
  :param str mass_balance_check: how check for mass balance
  """
  document = util.getSBMLDocument(model_reference)
  model = document.getModel()
  simple = SimpleSBML(model)
  Reaction.initialize(simple)
  print(MoietyComparator.analyzeReactions())


if __name__ == '__main__':
  lint(ANTIMONY_STG)
