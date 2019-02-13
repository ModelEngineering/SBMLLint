"""
Tests for simple_sbml
"""
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml
from SBMLLint.common import constants as cn

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
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.simple.getReactions()), cn.NUM_REACTIONS)
    for reaction in self.simple.getReactions():
      self.assertTrue(isinstance(reaction, tesbml.Reaction))
      self.assertLessEqual(reaction.getNumReactants(), cn.MAX_REACTANTS)
    self.assertEqual(len(self.simple.getParameters()), cn.NUM_PARAMETERS)
    for param in self.simple.getParameters():
      self.assertTrue(isinstance(param, str))

  def testGetSpecies(self):
    species = self.simple._getSpecies()
    trues = [isinstance(s, tesbml.libsbml.Species)
        for s in species.values()]
    self.assertTrue(all(trues))

  def testGetReactions(self):
    reactions = self.simple._getReactions()
    trues = [isinstance(r, tesbml.libsbml.Reaction)
        for r in reactions]
    self.assertTrue(all(trues))

  def testGetParameters(self):
    parameters = self.simple._getParameters()
    trues = [isinstance(s, tesbml.libsbml.Parameter)
        for s in parameters.values()]
    self.assertTrue(all(trues))

  def testGetReactants(self):
    for reaction in self.simple.getReactions():
      reactants = self.simple.getReactants(reaction)
      trues = [isinstance(r, tesbml.libsbml.SpeciesReference)
          for r in reactants]
      self.assertTrue(all(trues))

  def testGetProducts(self):
    for reaction in self.simple.getReactions():
      products = self.simple.getProducts(reaction)
      trues = [isinstance(r, tesbml.libsbml.SpeciesReference)
          for r in products]
      self.assertTrue(all(trues))

  def testGetReactionKineticsTerms(self):
    for reaction in self.simple.getReactions():
      stgs = SimpleSBML.getReactionKineticsTerms(reaction)
      trues = [isinstance(s, str) for s in stgs]
      self.assertTrue(trues)

  def testGetReactionString(self):
    for reaction in self.simple.getReactions():
      stg = SimpleSBML.getReactionString(reaction)
      parts = stg.split('->')
      self.assertTrue(";" in parts[-1])  # Kinetics is last

  def testIsSpecies(self):
    species = list(self.simple._getSpecies().keys())
    self.assertTrue(self.simple.isSpecies(species[0]))
    self.assertFalse(self.simple.isSpecies("dummy"))

  def testIsParameter(self):
    parameters = list(self.simple._getParameters().keys())
    self.assertTrue(self.simple.isParameter(parameters[0]))
    self.assertFalse(self.simple.isParameter("dummy"))


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
    for item in itr:
      self.assertTrue(isinstance(item.filename, str))
      self.assertTrue(isinstance(item.model, tesbml.libsbml.Model))
      self.assertEqual(item.number, COUNT - 1)
      #self.assertTrue('Model' in  str(type(item.model)))
      #self.assertEqual(item.number, COUNT - 1)
    

if __name__ == '__main__':
  unittest.main()
