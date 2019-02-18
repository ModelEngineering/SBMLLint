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
    def test(document):
      model = document.getModel()
      self.assertTrue('Reaction' in str(type(model.getReaction(0))))
    #
    test(util.getSBMLDocument(cn.TEST_FILE2))
    test(util.getSBMLDocument(ANTIMONY_STG))
    with open(cn.TEST_FILE2, 'r') as fd:
      lines = '\n'.join(fd.readlines())
    test(util.getSBMLDocument(lines))

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

  def testIsInt(self):
    self.assertTrue(util.isInt(1))
    self.assertFalse(util.isInt(1.5))
    self.assertFalse(util.isInt('ab'))

  def testIsFloat(self):
    self.assertTrue(util.isFloat(1))
    self.assertTrue(util.isFloat(1.5))
    self.assertTrue(util.isFloat('1.5'))
    self.assertFalse(util.isFloat('ab'))


if __name__ == '__main__':
  unittest.main()
