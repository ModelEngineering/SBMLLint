"""
Tests for Moiety
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml
from SBMLLint.structured_names.moiety import Moiety

import itertools
import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
MOIETY_NAME3 = "third"
NAMES = [MOIETY_NAME1, MOIETY_NAME2, MOIETY_NAME3]
MOLECULE_NAME = "%s%s%s" % (MOIETY_NAME1, cn.MOIETY_SEPARATOR, 
    MOIETY_NAME2)
iterator = itertools.product([0,1], [0, 1], [0, 1])
MOLECULE_NAME_SET = []  # A set of names from moiety combinations
for item in iterator:
  name = ""
  for idx,ele in enumerate(item):
    if ele == 1:
      if len(name) == 0:
        name = NAMES[idx]
      else:
        name = "%s%s%s" % (name, cn.MOIETY_SEPARATOR, NAMES[idx])
  if len(name) > 0:
    MOLECULE_NAME_SET.append(name)


#############################
# Tests
#############################
class TestMoiety(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(Moiety(MOLECULE_NAME).name, MOLECULE_NAME)

  def testExtract(self):
    if IGNORE_TEST:
      return
    moiety1 = Moiety(MOIETY_NAME1)
    moiety2 = Moiety(MOIETY_NAME2)
    molecule = Molecule(MOLECULE_NAME)
    names = set([m.name for m in Moiety.extract(molecule)])
    self.assertEqual(names, set([MOIETY_NAME1, MOIETY_NAME2]))

  def testAppendToMolecule(self):
    if IGNORE_TEST:
      return
    moiety2 = Moiety(MOIETY_NAME2)
    molecule = Molecule(MOIETY_NAME1)
    new_molecule = moiety2.appendToMolecule(molecule)
    self.assertEqual(new_molecule.name, Molecule(MOLECULE_NAME).name)

  def testCountMoietys(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(MOLECULE_NAME)
    df = Moiety.countMoietys([molecule])
    self.assertEquals(df.columns.tolist(), [cn.VALUE])
    df2 = Moiety.countMoietys([molecule, molecule])
    self.assertTrue(df2.equals(df + df))


if __name__ == '__main__':
  unittest.main()
