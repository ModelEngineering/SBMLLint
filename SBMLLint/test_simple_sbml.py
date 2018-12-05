"""
Tests for simple_sbml
"""
import unittest
import numpy as np
import tesbml
import simple_sbml
from simple_sbml import SimpleSBML
#import tesedml


IGNORE_TEST = True
TEST_FILE = "chemotaxis.xml"
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10


#############################
# Tests
#############################
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    pass
    #self.simple = SimpleSBML(TEST_FILE)

  def testConstructor(self):
    generator = simple_sbml.biomodelIterator(final=1)
    idx, model = [x for x in generator][0]
    simple = SimpleSBML(model)
    #
    self.assertEqual(len(simple.getReactions()), NUM_REACTIONS)
    self.assertEqual(len(self.simple.getReactions()), NUM_REACTIONS)
    for reaction in self.simple.getReactions():
      self.assertTrue(isinstance(reaction, tesbml.Reaction))
      self.assertLessEqual(reaction.getNumReactants(), MAX_REACTANTS)
    self.assertEqual(len(self.simple.getParameters()), NUM_PARAMETERS)



class TestFunctions(unittest.TestCase):

  def testBiomodelIterator(self):
    if IGNORE_TEST:
      return
    FINAL = 2
    generator = simple_sbml.biomodelIterator(final=FINAL)
    for idx, model in generator:
      self.assertGreaterEqual(FINAL, idx)
      #self.assertTrue(isinstance(model, tesedml.libsedml.Model))

if __name__ == '__main__':
  unittest.main()
