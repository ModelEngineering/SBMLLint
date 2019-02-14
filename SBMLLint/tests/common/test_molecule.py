"""
Tests for Molecule
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = True
NAME = "name"
NONAME = 'not_a_name'


#############################
# Tests
#############################
class TestMolecule(unittest.TestCase):

  def setUp(self):
    Molecule.molecules = []
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    molecules = []
    molecule = Molecule(NAME, other_molecules=molecules)
    self.assertEqual(molecule.name, NAME)
    self.assertEqual(molecules, [molecule])

  def testGetMolecule(self):
    _ = Molecule(NAME)
    molecule = Molecule.getMolecule(NAME)
    self.assertEqual(molecule, Molecule.molecules[0])
    self.assertIsNone(Molecule.getMolecule(NONAME))

  def testInitialize(self):
    if IGNORE_TEST:
      return
    Molecule.initialize(self.simple)
    self.assertEqual(len(Molecule.molecules), cn.NUM_SPECIES)
    

if __name__ == '__main__':
  unittest.main()
