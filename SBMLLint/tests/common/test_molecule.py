"""
Tests for Molecule
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = True
NUM1 = 2
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
NAME = "name"
BAD_NAME = 'not_a_name'
MOLECULE_NAME = "%s%s%s" % (MOIETY_NAME1, cn.MOIETY_DOUBLE_SEPARATOR, 
    MOIETY_NAME2)


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
    with self.assertRaises(ValueError):
      Molecule(BAD_NAME)

  def testGetMolecule(self):
    _ = Molecule(NAME)
    molecule = Molecule.getMolecule(NAME)
    self.assertEqual(molecule, Molecule.molecules[0])
    self.assertIsNone(Molecule.getMolecule(BAD_NAME))

  def testAppend(self):
    if IGNORE_TEST:
      return
    moiety2 = Moiety(MOIETY_NAME2)
    molecule = Molecule(MOIETY_NAME1)
    new_molecule = molecule.append(moiety2)
    self.assertEqual(new_molecule.name, Molecule(MOLECULE_NAME).name)

  def testHasMoiety(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(MOLECULE_NAME)
    self.assertTrue(molecule.hasMoiety(Moiety(MOIETY_NAME1)))
    self.assertTrue(molecule.hasMoiety(Moiety(MOIETY_NAME2)))


class TestMoleculeStoichiometry(unittest.TestCase):

  def setUp(self):
    Molecule.molecules = []

  def testConstructor(self):
    molecule = Molecule(MOLECULE_NAME)
    m_s = MoleculeStoichiometry(molecule, NUM1)
    self.assertTrue(m_s.molecule.isEqual(molecule))
    self.assertEqual(m_s.stoichiometry, NUM1)
  
  def testExtractMoietyStoichiometrys(self):
    if IGNORE_TEST:
      return
    for expected, strings in MOLECULE_STOICHIOMETRY_STGS.items():
      moietys = list(set([e[0] for e in expected]))
      moietys.sort()
      for stg in strings:
        result = Moiety.extract(stg)
        trues = [result[n].isEqual(moietys[n]) 
            for n, _ in enumerate(result)]
        self.assertTrue(all(trues))

  # TODO: Update tests
  def testCountMoietys(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(MOLECULE_NAME)
    df = MoleculeStoichiometry.countMoietys([molecule])
    self.assertEquals(df.columns.tolist(), [cn.VALUE])
    df2 = MoietyStoichiometry.countMoietys([molecule, molecule])
    self.assertTrue(df2.equals(df + df))

  def testInitialize(self):
    if IGNORE_TEST:
      return
    Molecule.initialize(self.simple)
    self.assertEqual(len(Molecule.molecules), cn.NUM_SPECIES)
    

if __name__ == '__main__':
  unittest.main()
