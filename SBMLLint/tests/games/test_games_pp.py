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
import unittest

IGNORE_TEST = False

ZERO = 0
ONE = 1
ZERO_F = 0.0
ONE_F = 1.0
# BIOMD0000000383
# Molecule names
PGA = "PGA"
RUBP = "RuBP"
CO2 = "CO2"
# Reaction names
PGA_CONS = "PGA_cons"
PGA_CONS_SOMREACTION_IDENTIFIER = "PGA_cons: {PGA} -> {RuBP}"
PGA_PROD_VC = "PGA_prod_Vc"
PGA_PROD_VO = "PGA_prod_Vo"
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
#"
REACTION = "reaction"
NUM_REACTIONS = 4
NUM_MOLECULES = 5
#
# BIOMD0000000018
# Reactions
CH2FH4toHCHO = "CH2FH4toHCHO"
HCHOtoCH2FH4 = "HCHOtoCH2FH4"
# Molecules
CH2FH4 = "CH2FH4"
FH4 = "FH4"
HCHO = "HCHO"

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
    # BIOMD0000000383
    self.simple1 = SimpleSBML()
    self.simple1.initialize(cn.TEST_FILE_GAMES_PP1)
    self.games_pp = GAMES_PP(self.simple1)
    # BIOMD0000000018
    self.simple2 = SimpleSBML()
    self.simple2.initialize(cn.TEST_FILE_GAMES_PP2)

  def testConstructor(self):
    if IGNORE_TEST:
      return
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
    if IGNORE_TEST:
      return
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
    if IGNORE_TEST:
      return
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
    if IGNORE_TEST:
      return
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
    self.assertTrue(echelon.iloc[echelon.shape[0]-1][0]==ZERO_F)
    # On the ohter hand, upper left should be nonzero
    self.assertFalse(echelon.iloc[0][0]==ZERO_F)

  def testGetRREFMatrix(self):
    if IGNORE_TEST:
      return
    self.assertTrue(self.games_pp.rref_operation is None) 
    self.assertTrue(self.games_pp.rref_df is None)
    som_reactions = []
    for r in self.games_pp.reactions:
      som_reactions.append(self.games_pp.convertReactionToSOMReaction(r))
    som_mat = self.games_pp.getStoichiometryMatrix(
        som_reactions,
        self.games_pp.nodes,
        som=True)
    echelon = self.games_pp.decomposeMatrix(som_mat).T
    rref = self.games_pp.getRREFMatrix(echelon)
    self.assertTrue(isinstance(self.games_pp.rref_operation, pd.DataFrame))
    self.assertTrue(isinstance(self.games_pp.rref_df, pd.DataFrame))
    # Choose the last row of rref.T 
    last_row = rref.T.iloc[-1:]
    # Get the first nonzero index
    nonzero_species = (last_row != 0).idxmax(axis=1)[0]
    # The SUM of nonzero species column must be the same as the single value
    self.assertEqual(last_row[nonzero_species][0], sum(rref.T[nonzero_species]))

  def testConverMatrixToSOMReactions(self):
    if IGNORE_TEST:
      return
    som_reactions1 = []
    for r in self.games_pp.reactions:
      som_reactions1.append(self.games_pp.convertReactionToSOMReaction(r))
    som_mat = self.games_pp.getStoichiometryMatrix(
        som_reactions1,
        self.games_pp.nodes,
        som=True)
    som_reactions2 = self.games_pp.convertMatrixToSOMReactions(som_mat)
    sr1 = None
    sr2 = None
    for sr in som_reactions1:
      if sr.label == PGA_PROD_VC:
        sr1 = sr
        break
    for sr in som_reactions2:
      if sr.label == PGA_PROD_VC:
        sr2 = sr
        break
    sr1_reactants = {ms.som for ms in sr1.reactants}
    sr2_reactants = {ms.som for ms in sr2.reactants}
    sr1_products = {ms.som for ms in sr1.products}
    sr2_products = {ms.som for ms in sr2.products}
    sr1_reactsom = sum([ms.stoichiometry for ms in sr1.reactants])
    sr2_reactsom = sum([ms.stoichiometry for ms in sr2.reactants])
    sr1_prodsom = sum([ms.stoichiometry for ms in sr1.products])
    sr2_prodsom = sum([ms.stoichiometry for ms in sr2.products])
    self.assertEqual(sr1_reactants, sr2_reactants)
    self.assertEqual(sr1_products, sr2_products)
    self.assertEqual(sr1_reactsom, sr2_reactsom)
    self.assertEqual(sr1_prodsom, sr2_prodsom)

  def testGetNode(self):
    if IGNORE_TEST:
      return
    co2 = self.games_pp.simple.getMolecule(CO2)
    co2_node = self.games_pp.getNode(co2)
    self.assertEqual(type(co2_node), SOM)
    self.assertEqual(co2_node.molecules, {co2})

  def testMergeNodes(self):
    if IGNORE_TEST:
      return
    reaction = self.games_pp.simple.getReaction(PGA_CONS)
    pga = self.games_pp.simple.getMolecule(PGA)
    rubp = self.games_pp.simple.getMolecule(RUBP)
    som_pga = self.games_pp.getNode(pga)
    som_rubp = self.games_pp.getNode(rubp)
    som_pga_rubp = self.games_pp.mergeNodes(som_pga, som_rubp, reaction)
    self.assertTrue(isinstance(som_pga_rubp, SOM))
    self.assertTrue(pga in som_pga_rubp.molecules)
    self.assertTrue(rubp in som_pga_rubp.molecules)
    self.assertTrue(reaction in som_pga_rubp.reactions)

  def testAddReaction(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.games_pp.reactions_lu), 0)
    reaction = self.games_pp.simple.getReaction(PGA_CONS)
    self.games_pp.addReaction(reaction)
    self.assertTrue(reaction in self.games_pp.reactions_lu)
    self.games_pp.addReaction(reaction)
    self.assertTrue(len(self.games_pp.reactions_lu), 1)

  def testProcessUniUniReaction(self):
    if IGNORE_TEST:
      return
    reaction = self.games_pp.simple.getReaction(PGA_CONS)
    som = self.games_pp.processUniUniReaction(reaction)
    self.assertTrue(isinstance(som, SOM))
    self.assertTrue(som in self.games_pp.nodes)
    pga = self.games_pp.simple.getMolecule(PGA)
    rubp = self.games_pp.simple.getMolecule(RUBP)
    self.assertTrue(pga in som.molecules)
    self.assertTrue(rubp in som.molecules)

  def testAddArc(self):
    if IGNORE_TEST:
      return
    reaction = self.games_pp.simple.getReaction(PGA_PROD_VC)
    co2 = self.games_pp.simple.getMolecule(CO2)
    pga = self.games_pp.simple.getMolecule(PGA)
    som_co2 = self.games_pp.getNode(co2)
    som_pga = self.games_pp.getNode(pga)
    self.games_pp.addArc(som_co2, som_pga, reaction)
    self.assertTrue(self.games_pp.has_edge(som_co2, som_pga))
    self.assertEqual(self.games_pp.get_edge_data(som_co2, som_pga)[REACTION], [PGA_PROD_VC])

  def testProcessUniMultiReaction(self):
    if IGNORE_TEST:
      return
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(isinstance(games_pp2, GAMES_PP))
    reaction = games_pp2.simple.getReaction(CH2FH4toHCHO)
    reactant1 = games_pp2.simple.getMolecule(CH2FH4)
    product1 = games_pp2.simple.getMolecule(FH4)
    product2 = games_pp2.simple.getMolecule(HCHO)
    games_pp2.processUniMultiReaction(reaction)
    som_reactant1 = games_pp2.getNode(reactant1)
    som_product1 = games_pp2.getNode(product1)
    som_product2 = games_pp2.getNode(product2)
    self.assertTrue(games_pp2.has_edge(som_product1, som_reactant1))
    self.assertTrue(games_pp2.has_edge(som_product2, som_reactant1))

  def testProcessMultiUniReaction(self):
    if IGNORE_TEST:
      return
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(isinstance(games_pp2, GAMES_PP))
    reaction = games_pp2.simple.getReaction(HCHOtoCH2FH4)
    reactant1 = games_pp2.simple.getMolecule(FH4)
    reactant2 = games_pp2.simple.getMolecule(HCHO)
    product1 = games_pp2.simple.getMolecule(CH2FH4)
    games_pp2.processMultiUniReaction(reaction)
    som_product1 = games_pp2.getNode(product1)
    som_reactant1 = games_pp2.getNode(reactant1)
    som_reactant2 = games_pp2.getNode(reactant2)
    self.assertTrue(games_pp2.has_edge(som_reactant1, som_product1))
    self.assertTrue(games_pp2.has_edge(som_reactant2, som_product1))

  def testProcessEqualSOMReaction(self):
    if IGNORE_TEST:
      return
    print("type three", self.games_pp.type_three_errors)
    self.assertEqual(len(self.games_pp.type_three_errors), ZERO)
    # Add an arc
    reaction = self.games_pp.simple.getReaction(PGA_PROD_VC)
    co2 = self.games_pp.simple.getMolecule(CO2)
    pga = self.games_pp.simple.getMolecule(PGA)
    som_co2 = self.games_pp.getNode(co2)
    som_pga = self.games_pp.getNode(pga)
    self.games_pp.addArc(som_co2, som_pga, reaction)
    # Process a dummy reaction
    soms_co2 = SOMStoichiometry(som_co2, 1.0)
    soms_pga = SOMStoichiometry(som_pga, 1.0)
    som_reaction = SOMReaction([soms_co2], [soms_pga], "dummy")
    self.games_pp.processEqualSOMReaction(som_reaction)
    self.assertTrue(len(self.games_pp.type_three_errors) > ZERO)

  def testProcessUnequalSOMReaction(self):
    if IGNORE_TEST:
      return
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(len(games_pp2.edges)==0)
    reaction = games_pp2.simple.getReaction(CH2FH4toHCHO)
    som_reaction = games_pp2.convertReactionToSOMReaction(reaction)
    games_pp2.processUnequalSOMReaction(som_reaction)
    self.assertTrue(len(games_pp2.edges)==2)
    reactant1 = games_pp2.simple.getMolecule(FH4)
    reactant2 = games_pp2.simple.getMolecule(HCHO)
    product1 = games_pp2.simple.getMolecule(CH2FH4)
    games_pp2.processUniMultiReaction(reaction)
    som_product1 = games_pp2.getNode(product1)
    som_reactant1 = games_pp2.getNode(reactant1)
    som_reactant2 = games_pp2.getNode(reactant2)
    self.assertTrue(games_pp2.has_edge(som_reactant1, som_product1))
    self.assertTrue(games_pp2.has_edge(som_reactant2, som_product1))

  def testAddTypeOneError(self):
    if IGNORE_TEST:
      return
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(len(games_pp2.type_one_errors)==ZERO)
    reaction = games_pp2.simple.getReaction(CH2FH4toHCHO)
    reactant1 = games_pp2.simple.getMolecule(FH4)
    product1 = games_pp2.simple.getMolecule(CH2FH4)
    games_pp2.addTypeOneError(reactant1, product1, reaction)
    self.assertTrue(len(games_pp2.type_one_errors)==ONE)
    error = games_pp2.type_one_errors[ZERO]
    self.assertEqual(error.node1, reactant1.name)
    self.assertEqual(error.node2, product1.name)
    self.assertEqual(error.reactions, [reaction.label])

  def testChekcTypeOneError(self):
    if IGNORE_TEST:
      return
    games_pp1 = GAMES_PP(self.simple1)
    reaction = games_pp1.simple.getReaction(PGA_CONS)
    pga = games_pp1.simple.getMolecule(PGA)
    rubp = games_pp1.simple.getMolecule(RUBP)
    som_pga = games_pp1.getNode(pga)
    som_rubp = games_pp1.getNode(rubp)
    som_pga_rubp = games_pp1.mergeNodes(som_pga, som_rubp, reaction)
    self.assertTrue(len(games_pp1.type_one_errors)==ZERO)
    is_error = games_pp1.checkTypeOneError((pga, rubp), reaction)
    self.assertTrue(is_error)
    error = games_pp1.type_one_errors[ZERO]
    self.assertEqual(error.node1, pga.name)
    self.assertEqual(error.node2, rubp.name)
    self.assertEqual(error.reactions, [reaction.label])

  def testCheckTypeTwoError(self):
    if IGNORE_TEST:
      return    
    # Create a dummy cycle by adding two conflicting arcs
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(isinstance(games_pp2, GAMES_PP))
    unimulti_reaction = games_pp2.simple.getReaction(CH2FH4toHCHO)
    multiuni_reaction = games_pp2.simple.getReaction(HCHOtoCH2FH4)
    fh4 = games_pp2.simple.getMolecule(FH4)
    # hcho = games_pp2.simple.getMolecule(HCHO)
    ch2fh4 = games_pp2.simple.getMolecule(CH2FH4)
    som_fh4 = games_pp2.getNode(fh4)
    som_ch2fh4 = games_pp2.getNode(ch2fh4)
    # do we need the next two methods if we're giving a cycle? 
    games_pp2.addArc(som_fh4, som_ch2fh4, unimulti_reaction)
    games_pp2.addArc(som_ch2fh4, som_fh4, multiuni_reaction) 
    self.assertTrue(len(games_pp2.type_two_errors)==ZERO)
    games_pp2.checkTypeTwoError()
    self.assertTrue(len(games_pp2.type_two_errors)==ONE)
    error = games_pp2.type_two_errors[ZERO]
    self.assertTrue(games_pp2.has_edge(error[ZERO], error[ONE]))
    self.assertTrue(games_pp2.has_edge(error[ONE], error[ZERO]))
    error_reactions = set(games_pp2.get_edge_data(error[ZERO], error[ONE])[REACTION])
    error_reactions = error_reactions.union(set(games_pp2.get_edge_data(error[ONE], error[ZERO])[REACTION]))
    self.assertTrue(CH2FH4toHCHO in error_reactions)
    self.assertTrue(HCHOtoCH2FH4 in error_reactions)

  def testProcessErrorReaction(self):
    if IGNORE_TEST:
      return  
    games_pp1 = GAMES_PP(self.simple1)
    self.assertTrue(len(games_pp1.echelon_errors)==ZERO)
    reaction = self.games_pp.simple.getReaction(PGA_PROD_VO)
    games_pp1.processErrorReaction(reaction)
    self.assertTrue(len(games_pp1.echelon_errors)==ONE)
    self.assertTrue(reaction in games_pp1.echelon_errors)

  def testAnalyze(self):
    if IGNORE_TEST:
      return  
    games_pp1 = GAMES_PP(self.simple1)
    games_pp2 = GAMES_PP(self.simple2)
    self.assertTrue(games_pp1.analyze(rref=False))
    self.assertTrue(games_pp2.analyze())
    self.assertTrue(len(games_pp1.echelon_errors)>ZERO)
    self.assertTrue(len(games_pp1.type_one_errors)==ZERO)
    self.assertTrue(len(games_pp1.type_two_errors)==ZERO)
    self.assertTrue(len(games_pp2.type_one_errors)>ZERO)

if __name__ == '__main__':
  unittest.main()





















































