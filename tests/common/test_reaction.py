"""
Tests for Reactions
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.reaction import Reaction, REACTION_SEPARATOR

import numpy as np
import os
import unittest


IGNORE_TEST = False
NUM_S1 = 2
NUM_S2 = 3


#############################
# Tests
#############################
class TestReaction(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE)
    self.reactions = self.simple.reactions
    self.reaction = self.reactions[2]

  def testGetId(self):
    identifier1 = self.reaction.getId()
    self.assertTrue(REACTION_SEPARATOR in identifier1)
    identifier2 = self.reaction.getId(is_include_kinetics=False)
    self.assertGreater(len(identifier1), len(identifier2))
    self.assertFalse(cn.KINETICS_SEPARATOR in identifier2)
    identifier3 = self.reaction.getId(is_include_kinetics=False,
        is_include_label=False)
    self.assertGreater(len(identifier2), len(identifier3))
    self.assertFalse(cn.LABEL_SEPARATOR in identifier3)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.reaction.reactants[0],
        MoleculeStoichiometry))
    self.assertTrue(isinstance(self.reaction.products[0],
        MoleculeStoichiometry))
    self.assertEqual(self.reaction.category, cn.REACTION_1_n)
    self.assertGreater(len(self.simple.molecules), 0)
    count = len(self.simple.reactions)
    reaction = self.simple.reactions[3]
    self.assertEqual(len(self.simple.reactions), count)

  # Test the stoichiometry
  def testConstructor2(self):
    if IGNORE_TEST:
      return
    simple = SimpleSBML()
    simple.initialize(cn.TEST_FILE2)
    reaction = simple.reactions[0]
    self.assertEqual(reaction.reactants[0].stoichiometry, NUM_S1)
    self.assertEqual(reaction.products[0].stoichiometry, NUM_S2)

  def testIsEqual(self):
    if IGNORE_TEST:
      return
    for reaction in self.simple.reactions[1:]:
      self.assertTrue(reaction.isEqual(reaction))
      self.assertFalse(reaction.isEqual(self.simple.reactions[0]))

  def testMakeIdentifier(self):
    if IGNORE_TEST:
      return
    for reaction in self.simple.reactions:
      parts = reaction.identifier.split('->')
      self.assertTrue(";" in parts[-1])  # Kinetics is last

  def testFindReactions(self):
    if IGNORE_TEST:
      return
    reactions = Reaction.find(self.reactions,
        category=cn.REACTION_1_n)
    trues = [r.category == cn.REACTION_1_n for r in reactions]
    self.assertTrue(all(trues))


if __name__ == '__main__':
  unittest.main()
