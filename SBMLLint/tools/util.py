""" Utilities for tools."""

import SBMLLint.common.constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.reaction import Reaction
import SBMLLint.common.util as util

def initialize(model_reference):
  """
  :param str model_reference: source of the SBML model
     file, antimony string, xml string
  """
  document = util.getSBMLDocument(model_reference)
  model = document.getModel()
  simple = SimpleSBML(model)
  Reaction.initialize(simple)
