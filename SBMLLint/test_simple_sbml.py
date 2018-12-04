"""
Tests for simple_sbml
"""
import unittest
import numpy as np
import tesbml
import simple_sbml
from simple_sbml import SimpleSBML
import tellurium as te


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



class TestFunctions(unittest.TestCase):

  def testBiomodelIterator(self):
    FINAL = 2
    generator = simple_sbml.biomodelIterator(final=FINAL)
    import pdb; pdb.set_trace()
    for idx, model in generator:
      import pdb; pdb.set_trace()
      self.assertGreater(FINAL, idx)
      self.assertTrue(isinstance(model, te.libsbml.model))

if __name__ == '__main__':
  unittest.main()
