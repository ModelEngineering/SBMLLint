"""
Tests for Reactions
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction, REACTION_SEPARATOR
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False

def strLen(a_list):
  return sum([len(x.name) for x in a_list])

#############################
# Tests
#############################
class TestMolecule(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)
    self.reactions = self.simple.reactions
    Reaction.reactions = []

  def testMakeId(self):
    self.reaction = Reaction(self.reactions[2])
    identifier = self.reaction.makeId()
    self.assertTrue(REACTION_SEPARATOR in identifier)
    self.assertGreater(len(identifier), len(REACTION_SEPARATOR))

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.reaction = Reaction(self.reactions[2])
    self.assertEqual(strLen(self.reaction.reactants),
        strLen(self.reaction.products))
    self.assertEqual(self.reaction.category, cn.REACTION_1_n)
    self.assertGreater(len(Molecule.molecules), 0)
    count = len(Reaction.reactions)
    reaction = Reaction(self.reactions[3])
    self.assertEqual(len(Reaction.reactions), count + 1)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(Reaction.reactions), 0)
    Reaction.initialize(self.simple)
    self.assertEqual(len(Reaction.reactions), cn.NUM_REACTIONS)


if __name__ == '__main__':
  unittest.main()
