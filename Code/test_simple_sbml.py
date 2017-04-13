"""
Tests for simple_sbml
"""
import unittest
import numpy as np
from simple_sbml import SimpleSBML
import libsbml


IGNORE_TEST = False
TEST_FILE = "chemotaxis.xml"
NUM_REACTIONS = 111
NUM_PARAMETERS = 27


#############################
# Tests
#############################
# pylint: disable=W0212,C0111,R0904
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.simple.getReactions()), NUM_REACTIONS)
    for reaction in self.simple.getReactions():
      self.assertTrue(isinstance(reaction, libsbml.Reaction))
    self.assertEqual(len(self.simple.getParameters()), NUM_PARAMETERS)
    for reaction in self.simple.getParameters():
      self.assertTrue(isinstance(reaction, libsbml.Parameter))


if __name__ == '__main__':
  unittest.main()
