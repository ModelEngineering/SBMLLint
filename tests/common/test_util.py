from SBMLLint.common import constants as cn
from SBMLLint.common import exceptions
from SBMLLint.common import util
import libsbml


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
ZIP_PATH = os.path.join(cn.BIOMODELS_DIR, cn.BIOMODELS_ZIP_FILENAME)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testGetXMLString(self):
    def test(xml):
      reader = libsbml.SBMLReader()
      document = reader.readSBMLFromString(xml)
      util.checkSBMLDocument(document)
      model = document.getModel()
      self.assertTrue('Reaction' in str(type(model.getReaction(0))))
    def getString(path):
      with open(path, 'r') as fd:
        lines = '\n'.join(fd.readlines())
      return lines
    #
    for path in [cn.TEST_FILE2, cn.TEST_FILE3]:
      try:
        test(util.getXML(path))
        test(util.getXML(getString(path)))
      except exceptions.MissingTelluriumError:
        pass

  def testGetXMLFromAntimony(self):
    try:
      xml = util.getXMLFromAntimony(ANTIMONY_STG)
    except exceptions.MissingTelluriumError:
      return
    self.assertTrue(isinstance(xml, str))
    reader = libsbml.SBMLReader()
    libsbml_document = reader.readSBMLFromString(xml)
    util.checkSBMLDocument(libsbml_document)
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
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromString(xml)
    util.checkSBMLDocument(document)
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

  def testGetNextFid(self):
    fid = open(ZIP_PATH, "r")
    count = 0
    for zip_fid in util.getNextFid(fid, is_print=False):
      lines = zip_fid.read()
      count += 1
      self.assertGreater(len(lines), 0)
    self.assertGreater(count, 0)
    

if __name__ == '__main__':
  unittest.main()
