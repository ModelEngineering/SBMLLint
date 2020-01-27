"""
Tests for Set Of Moletules (SOM)
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.games.som import SOM
from SBMLLint.common import simple_sbml

import numpy as np
import os
import re
import unittest


IGNORE_TEST = False
NUM_MERGED_SOMS = 2
NAMEFILTER = "[0-9a-zA-Z]+"


#############################
# Tests
#############################
class TestSOM(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple.initialize(cn.TEST_FILE3)
    self.molecules = self.simple.molecules
    self.soms = []
    for mole in self.molecules:
      self.soms.append(SOM({mole}))
  
  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.soms), len(self.molecules))
    self.assertEqual(self.soms[0].molecules, {self.molecules[0]})
    
  def testMakeId(self):
    if IGNORE_TEST:
      return
    som = self.soms[0]
    self.assertTrue(som.molecules.intersection(self.molecules))
    molecule = list(re.findall(NAMEFILTER, som.identifier))[0]
    self.assertEqual(list(som.molecules)[0], self.simple.getMolecule(molecule))
  
  def testMerge(self):
    if IGNORE_TEST:
      return
    som1 = self.soms[0]
    som2 = self.soms[1]
    molecule1 = list(som1.molecules)[0]
    molecule2 = list(som2.molecules)[0]
    new_som = som1.merge(som2)
    self.assertEqual(len(new_som.molecules), NUM_MERGED_SOMS)
    self.assertTrue(molecule1 in new_som.molecules)
    self.assertTrue(molecule2 in new_som.molecules)

if __name__ == '__main__':
  unittest.main()
