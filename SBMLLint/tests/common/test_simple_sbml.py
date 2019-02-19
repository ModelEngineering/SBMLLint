"""
Tests for simple_sbml
"""
from SBMLLint.common import constants as cn
from SBMLLint.common import simple_sbml
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
ANTIMONY_STG = '''
2S1 -> 3S2; 1
S1 = 0
S2 = 0
'''


#############################
# Tests
#############################
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()

  def testInitialize(self):
    if IGNORE_TEST:
      return
    self.simple.initialize(cn.TEST_FILE)
    self.assertEqual(len(self.simple.reactions), cn.NUM_REACTIONS)
    self.assertEqual(len(self.simple.parameters), cn.NUM_PARAMETERS)
    import pdb; pdb.set_trace()


class TestFunctions(unittest.TestCase):

  def testReadURL(self):
    pass

  def testModelIterator(self):
    itr = simple_sbml.modelIterator(final=1)
    for item in itr:
      model = item.model
      self.assertTrue(isinstance(model.getSpecies(0),
          tesbml.libsbml.Species))
    COUNT = 20
    itr = simple_sbml.modelIterator(final=COUNT)
    item_number = -1
    for item in itr:
      self.assertTrue(isinstance(item.filename, str))
      self.assertTrue('Model' in  str(type(item.model)))
      item_number = item.number
    self.assertEqual(item_number, COUNT - 1)
    

if __name__ == '__main__':
  unittest.main()
