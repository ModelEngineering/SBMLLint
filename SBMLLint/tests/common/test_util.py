from SBMLLint.common import constants as cn
from SBMLLint.common import util
import tesbml


import numpy as np
import os
import unittest


NUM_S1 = 2
NUM_S2 = 3
IGNORE_TEST = False
ANTIMONY_STG = '''
%dS1 -> %dS2; 1
S1 = 0
S2 = 0
''' % (NUM_S1, NUM_S2)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testGetModel(self):
    document = util.getSBMLDocument(ANTIMONY_STG)
    model = document.getModel()
    self.assertTrue('Reaction' in str(type(model.getReaction(0))))

  def testGetModelFromAntimony(self):
    document = util.getSBMLStringFromAntimony(ANTIMONY_STG)
    self.assertTrue(isinstance(document, str))
    reader = tesbml.libsbml.SBMLReader()
    libsbml_document = reader.readSBMLFromString(document)
    if (libsbml_document.getNumErrors() > 0):
      raise IOError("Errors in SBML document\n%s" 
          % libsbml_document.printErrors())
    model = libsbml_document.getModel()
    self.assertTrue('Reaction' in 
       str(type(model.getReaction(0))))

  def testMakeSBMLFile(self):
    util.makeSBMLFile(ANTIMONY_STG)
    self.assertTrue(os.path.isfile(cn.TEST_FILE2))


if __name__ == '__main__':
  unittest.main()
