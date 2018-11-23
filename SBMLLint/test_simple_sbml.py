"""
Tests for simple_sbml
"""
import unittest
import numpy as np
from simple_sbml import SimpleSBML
import tesbml


IGNORE_TEST = False
TEST_FILE = "chemotaxis.xml"
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10


#############################
# Tests
#############################
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.simple.getReactions()), NUM_REACTIONS)
    for reaction in self.simple.getReactions():
      self.assertTrue(isinstance(reaction, tesbml.Reaction))
      self.assertLessEqual(reaction.getNumReactants(), MAX_REACTANTS)
    self.assertEqual(len(self.simple.getParameters()), NUM_PARAMETERS)


if __name__ == '__main__':
  unittest.main()
