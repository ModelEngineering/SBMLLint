"""
Notes
  1. libsbml (tesbml) is fragile; will segmentation fault if return
     a model object. Each test must separately construct libsbml
     objects.
"""

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

  def testGetXMLString(self):
    def test(xml):
      reader = tesbml.SBMLReader()
      document = reader.readSBMLFromString(xml)
      if (document.getNumErrors() > 0):
        raise ValueError("Errors in SBML document\n%s" 
            % model_reference)
      model = document.getModel()
      self.assertTrue('Reaction' in str(type(model.getReaction(0))))
    def getString(path):
      with open(path, 'r') as fd:
        lines = '\n'.join(fd.readlines())
      return lines
    #
    for path in [cn.TEST_FILE2, cn.TEST_FILE3]:
      test(util.getXML(path))
      test(util.getXML(getString(path)))

  def testGetXMLFromAntimony(self):
    xml = util.getXMLFromAntimony(ANTIMONY_STG)
    self.assertTrue(isinstance(xml, str))
    reader = tesbml.libsbml.SBMLReader()
    libsbml_document = reader.readSBMLFromString(xml)
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

  def testIsSBMLModel(self):
    return
    self.assertFalse(util.isSBMLModel("dummy"))
    xml = util.getXML(cn.TEST_FILE2)
    reader = tesbml.SBMLReader()
    document = reader.readSBMLFromString(xml)
    if (document.getNumErrors() > 0):
      raise ValueError("Errors in SBML document\n%s" 
          % model_reference)
    model = document.getModel()
    self.assertTrue(util.isSBMLModel(model))

  def testUniqueify(self):
    class Tester():
    
      def __init__(self, name):
        self.name = name

      def __repr__(self):
        return self.name
 
      def isEqual(self, other):
        return self.name == other.name
    #
    STRING = 'abc'
    REPEATED_STRING = STRING + STRING   
    collection = [Tester(s) for s in REPEATED_STRING]
    result = util.uniqueify(collection)
    self.assertEqual(len(result), len(STRING))
    

if __name__ == '__main__':
  unittest.main()
