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


IGNORE_TEST = False
NAME = "name"


#############################
# Tests
#############################
class TestMolecule(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(Molecule(NAME).name, NAME)

  def testInitialize(self):
    if IGNORE_TEST:
      return
    Molecule.initialize(self.simple)
    self.assertEqual(len(Molecule.all), cn.NUM_SPECIES)
    

if __name__ == '__main__':
  unittest.main()
