"""
Tests for Moiety
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml
from SBMLLint.structured_name.moiety import Moiety

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
MOLECULE_NAME = "%s%s%s" % (MOIETY1, cn.MOIETY_SEPARATOR, MOIETY2)


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
    self.assertEqual(len(Molecule.molecules), cn.NUM_SPECIES)
    

if __name__ == '__main__':
  unittest.main()
