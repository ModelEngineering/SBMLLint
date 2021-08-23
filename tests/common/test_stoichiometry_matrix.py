"""
Tests for StoichiometryMatrix
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.som import SOM
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common.stoichiometry_matrix import StoichiometryMatrix

from scipy.optimize import linprog
import numpy as np
import os
import pandas as pd
import re
import unittest


IGNORE_TEST = False
REACTION1 = "reaction1"
REACTION4 = "reaction4"
MOLECULE_C = "C"
MOLECULE_M = "M"
REMAINIG_REACTIONS = 4
REMAINING_MOLECULES = 3


#############################
# Tests
#############################
class TestStoichiometryMatrix(unittest.TestCase):

  def setUp(self):
    self.simple_consistent = SimpleSBML()
    self.simple_consistent.initialize(cn.TEST_FILE8)
    self.consistent_matrix = StoichiometryMatrix(self.simple_consistent)
    #
    self.simple_inconsistent = SimpleSBML()
    self.simple_inconsistent.initialize(cn.TEST_FILE9)
    self.inconsistent_matrix = StoichiometryMatrix(self.simple_inconsistent)   
    #
    self.simple_repeated_species = SimpleSBML()
    self.simple_repeated_species.initialize(cn.TEST_FILE_REPEATED_SPECES)
    self.repeated_species_matrix = StoichiometryMatrix(self.simple_repeated_species)   

  def testConstructor(self):
    remaining_reactions = [r.label for r in self.consistent_matrix.reactions]
    print(remaining_reactions)
    remaining_molecules = self.consistent_matrix.molecules
    self.assertFalse(REACTION1 in remaining_reactions)
    self.assertTrue(REACTION4 in remaining_reactions)
    self.assertFalse(MOLECULE_C in remaining_molecules)
    self.assertTrue(MOLECULE_M in remaining_molecules)
    #
    self.assertEqual(len(self.inconsistent_matrix.reactions), REMAINIG_REACTIONS)
    self.assertEqual(len(self.inconsistent_matrix.molecules), REMAINING_MOLECULES)

  def testMakeStoichiometryMatrix(self):
    self.assertEqual(type(self.consistent_matrix), StoichiometryMatrix)
    self.assertEqual(type(self.inconsistent_matrix), StoichiometryMatrix)
    self.assertEqual(self.inconsistent_matrix.stoichiometry_matrix.shape, 
    	(REMAINING_MOLECULES, REMAINIG_REACTIONS))
    # tests repeated species matrix
    self.assertEqual(type(self.repeated_species_matrix), StoichiometryMatrix)
    rs_matrix = self.repeated_species_matrix.stoichiometry_matrix
    self.assertEqual(rs_matrix.shape, (6, 4))
    self.assertEqual(rs_matrix.loc['S0', 'J0'], -1.0)
    self.assertEqual(rs_matrix.loc['S1', 'J0'], 1.0)
    self.assertEqual(rs_matrix.loc['S1', 'J1'],-2.0)
    self.assertEqual(rs_matrix.loc['S2', 'J1'], 1.0)
    self.assertEqual(rs_matrix.loc['S2', 'J2'], -1.0)
    self.assertEqual(rs_matrix.loc['S3', 'J2'], 2.0)
    self.assertEqual(rs_matrix.loc['S4', 'J3'], -2.0)
    self.assertEqual(rs_matrix.loc['S5', 'J3'], 2.0)

  def testIsConsistent(self):
    self.assertTrue(self.consistent_matrix.consistent is None)
    self.assertTrue(self.inconsistent_matrix.consistent is None)
    self.assertTrue(self.consistent_matrix.isConsistent())
    self.assertFalse(self.inconsistent_matrix.isConsistent())

if __name__ == '__main__':
  unittest.main()
    