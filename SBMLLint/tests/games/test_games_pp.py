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
import pandas as pd
import os
import re
import tesbml
import unittest

IGNORE_TEST = False

ZERO = 0.0
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
# SOMStoichiometry NADPH identifier and stoichioemtry
SOM_NADPH = "{NADPH}"
SOM_PGA = "{PGA}"
SOM_RUBP = "{RuBP}"
NADPH_STOICHIOMETRY = 2.0
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

  def testConvertReactionToSOMReaction(self):
    reaction1 = self.games_pp.simple.getReaction(PGA_CONS)
    reaction2 = self.games_pp.simple.getReaction(PGA_PROD_VC)
    som_reaction1 = self.games_pp.convertReactionToSOMReaction(reaction1)
    som_reaction2 = self.games_pp.convertReactionToSOMReaction(reaction2)
    self.assertTrue(isinstance(som_reaction1, SOMReaction))
    self.assertTrue(isinstance(som_reaction2, SOMReaction))
    ms_rubp = som_reaction1.products[0]
    self.assertTrue(isinstance(ms_rubp, SOMStoichiometry))
    self.assertEqual(list(ms_rubp.som.molecules)[0], 
        self.games_pp.simple.getMolecule(RUBP))
    ms_nadph = None
    for ms in som_reaction2.reactants:
      if ms.som.identifier == SOM_NADPH:
        ms_nadph = ms
        break
    self.assertTrue(ms_nadph.som.identifier==SOM_NADPH)
    self.assertEqual(ms_nadph.stoichiometry, NADPH_STOICHIOMETRY)

  def testGetStoichiometryMatrix(self):
    # For regular stoichiometry matrix
    mat = self.games_pp.getStoichiometryMatrix(
        self.games_pp.reactions,
        self.games_pp.molecules,
        som=False)
    self.assertTrue(isinstance(mat, pd.DataFrame))
    self.assertEqual(mat.shape, (NUM_MOLECULES, NUM_REACTIONS))
    self.assertEqual(mat[PGA_CONS][PGA], PGA_CONS_WITH_PGA)
    self.assertEqual(mat[PGA_PROD_VC][RUBP], PGA_PROD_VC_WITH_RUBP)
    # For SOM stoichiometry matrix
    som_reactions = []
    for r in self.games_pp.reactions:
      som_reactions.append(self.games_pp.convertReactionToSOMReaction(r))
    som_mat = self.games_pp.getStoichiometryMatrix(
        som_reactions,
        self.games_pp.nodes,
        som=True)
    self.assertTrue(isinstance(som_mat, pd.DataFrame))
    self.assertEqual(som_mat[PGA_CONS][SOM_PGA], PGA_CONS_WITH_PGA)
    self.assertEqual(som_mat[PGA_PROD_VC][SOM_RUBP], PGA_PROD_VC_WITH_RUBP)

  def testDecomposeMatrix(self):
    # should be all None before decomposition
    self.assertTrue(self.games_pp.perm_inverse is None)
    self.assertTrue(self.games_pp.permuted_matrix is None)
    self.assertTrue(self.games_pp.lower_inverse is None)
    self.assertTrue(self.games_pp.echelon_df is None)
    som_reactions = []
    for r in self.games_pp.reactions:
      som_reactions.append(self.games_pp.convertReactionToSOMReaction(r))
    som_mat = self.games_pp.getStoichiometryMatrix(
        som_reactions,
        self.games_pp.nodes,
        som=True)
    echelon = self.games_pp.decomposeMatrix(som_mat).T
    self.assertFalse(self.games_pp.perm_inverse is None)
    self.assertFalse(self.games_pp.permuted_matrix is None)
    self.assertFalse(self.games_pp.lower_inverse is None)
    self.assertFalse(self.games_pp.echelon_df is None)
    self.assertTrue(isinstance(self.games_pp.perm_inverse, np.ndarray))    
    self.assertTrue(isinstance(self.games_pp.permuted_matrix, pd.DataFrame))
    self.assertTrue(isinstance(self.games_pp.lower_inverse, pd.DataFrame))
    self.assertTrue(isinstance(self.games_pp.echelon_df, pd.DataFrame))
    # Checking lower echelon - lower left element should be 0.0 (zero)
    self.assertTrue(echelon.iloc[echelon.shape[0]-1][0]==ZERO)
    # On the ohter hand, upper left should be nonzero
    self.assertFalse(echelon.iloc[0][0]==ZERO)




if __name__ == '__main__':
  unittest.main()





























