"""
Utilities that use Tellurium.
These utilities are isolated because of problems with using
tellurium in combination with libsbml.
"""

from SBMLLint.common import constants as cn

import tellurium as te
import tesbml

NUM_S1 = 2
NUM_S2 = 3
ANTIMONY_STG = '''
%dS1 -> %dS2; 1
S1 = 0
S2 = 0
''' % (NUM_S1, NUM_S2)
XML_HEADER = '<?xml version="1.0" encoding="UTF-8"?>'


def getSBMLFromAntimony(antimony_stg):
  """
  Constructs an SBML model from the antimony string.
  :param str antimony_stg:
  :return str: SBML model in xml format
  """
  rr = te.loada(antimony_stg)
  model_sbml = rr.getSBML()
  reader = tesbml.libsbml.SBMLReader()
  document = reader.readSBMLFromString(model_sbml)
  return document.toSBML()

def makeSBMLFile(stg=ANTIMONY_STG, path=cn.TEST_FILE2):
  document = getSBMLFromAntimony(stg)
  document = "%s\n%s" % (XML_HEADER, document)
  with open(path, 'w') as fd:
    fd.writelines(document)
