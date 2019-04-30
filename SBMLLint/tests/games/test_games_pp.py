"""
Test for GAMES Plus (Reduced) Row Echelon Form (GAMES_PP)
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.games_pp import SOMStoichiometry, SOMReaction, GAMES_PP
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import os
import re
import tesbml
import unittest

IGNORE_TEST = False

# BIOMD0000000383
# Molecule names
PGA = "PGA"
RUBP = "RuBP"
# Reaction names
PGA_CONS = "PGA_cons"
PGA_CONS_SOMREACTION_IDENTIFIER = "PGA_cons: {PGA} -> {RuBP}"
PGA_PROD_VC = "PGA_prod_Vc"
# SOMStoichiometry identifier
RUBP_ONE = "{RuBP} * 1.00"
# StoichiometryMatrix value
PGA_CONS_WITH_PGA = -1.0
PGA_PROD_VC_WITH_RUBP = -1.0
#
NUM_REACTIONS = 4
NUM_MOLECULES = 5


#############################
# Tests
#############################
class TestSOMStoichiometry(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE_GAMES_PP1)
    self.reaction = self.simple.getReaction(PGA_CONS)
    self.rubp = self.reaction.products[0]
    self.rubl_ss = SOMStoichiometry(
        som = SOM([self.rubp.molecule]),
        stoichiometry = self.rubp.stoichiometry
        )

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(isinstance(self.rubl_ss, SOMStoichiometry))
    self.assertTrue(isinstance(self.rubl_ss.som, SOM))
    self.assertTrue(isinstance(self.rubl_ss.stoichiometry, float))

  def testMakeId(self):
  	self.assertEqual(self.rubl_ss.identifier, RUBP_ONE)


class TestSOMReaction(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE_GAMES_PP1)
    self.reaction = self.simple.getReaction(PGA_CONS)
    self.pga = self.reaction.reactants[0]
    self.pga_ss = SOMStoichiometry(
        som = SOM([self.pga.molecule]),
        stoichiometry = self.pga.stoichiometry
        )
    self.rubp = self.reaction.products[0]
    self.rubp_ss = SOMStoichiometry(
        som = SOM([self.rubp.molecule]),
        stoichiometry = self.rubp.stoichiometry
        )
    self.som_reaction = SOMReaction(
        reactants=[self.pga_ss],
        products=[self.rubp_ss],
        label=self.reaction.label)

  def testConstructor(self):
    self.assertTrue(isinstance(self.reaction, Reaction))
    self.assertEqual(self.pga.molecule.name, PGA)
    self.assertEqual(self.rubp.molecule.name, RUBP)    
    self.assertTrue(isinstance(self.pga_ss, SOMStoichiometry))
    self.assertTrue(isinstance(self.rubp_ss, SOMStoichiometry))
    self.assertTrue(isinstance(self.som_reaction, SOMReaction))      

  def testMakeId(self):
  	self.assertEqual(
  	    self.som_reaction.makeId(),
  	    PGA_CONS_SOMREACTION_IDENTIFIER
  	    )

  def testGetCategory(self):
  	self.assertEqual(
  	    self.som_reaction.category,
  	    cn.REACTION_1_1)


class TestGAMES_PP(unittest.TestCase):

  def setUp(self):
    self.simple1 = SimpleSBML()
    self.simple1.initialize(cn.TEST_FILE_GAMES_PP1)
    self.games_pp = GAMES_PP(self.simple1)

  def testConstructor(self):
    self.assertEqual(len(self.games_pp.reactions), NUM_REACTIONS)
    for r in self.games_pp.reactions:
      print(r.category)
      self.assertFalse(r.category==cn.REACTION_BOUNDARY)
    self.assertEqual(len(self.games_pp.molecules), NUM_MOLECULES)
    for m in self.games_pp.molecules:
      self.assertTrue(isinstance(m, Molecule))
    self.assertEqual(len(self.games_pp.soms), NUM_MOLECULES)
    self.assertTrue(isinstance(self.games_pp.soms[0], SOM))
    init_identifier = ""
    for som in self.games_pp.soms:
      init_identifier = init_identifier + som.identifier
      if som != self.games_pp.soms[-1]:
        init_identifier = init_identifier + ";"
    self.assertEqual(self.games_pp.identifier, init_identifier)

  def testGetStoichiometryMatrix(self):
    mat = self.games_pp.getStoichiometryMatrix(
        self.games_pp.reactions,
        self.games_pp.molecules,
        som=False)
    self.assertEqual(mat.shape, (NUM_MOLECULES, NUM_REACTIONS))
    self.assertEqual(mat[PGA_CONS][PGA], PGA_CONS_WITH_PGA)
    self.assertEqual(mat[PGA_PROD_VC][RUBP], PGA_PROD_VC_WITH_RUBP)

if __name__ == '__main__':
  unittest.main()