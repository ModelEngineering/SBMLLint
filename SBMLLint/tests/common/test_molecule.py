"""
Tests for Molecule
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml

import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
NUM1 = 2
NUM2 = 3
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
NAME = "name"
BAD_NAME = 'not_a_name'
MOLECULE_NAME = "%s%s%s" % (MOIETY_NAME1, cn.MOIETY_DOUBLE_SEPARATOR, 
    MOIETY_NAME2)
MOLECULE_STOICHIOMETRY_STGS = {
    ((Moiety("A"), 1), (Moiety("P"), 1)): 
        [Molecule("A_P"), Molecule("A_1__P"), Molecule("A__P")],
    ((Moiety("AA"), 1), (Moiety("PP"), 1)): 
        [Molecule("AA_PP"), Molecule("AA_1__PP"), Molecule("AA__PP"),
        ],
    }


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
    self.simple = SimpleSBML(cn.TEST_FILE)
    Molecule.molecules = []

  def testConstructor(self):
    if IGNORE_TEST:
      return
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
      m_ss = [MoietyStoichiometry(m, 1) for m in moietys]
      for stg in strings:
        result = Molecule.extractMoietyStoichiometrys(stg)
        trues = [x.isEqual(y) for x, y in zip(result, m_ss)]
        self.assertTrue(all(trues))

  def testExtractMoietyStoichiometrys(self):
    if IGNORE_TEST:
      return
    moiety_stoich = MoietyStoichiometry(Molecule(MOIETY_NAME1),
        NUM1)
    mole_stoich = MoleculeStoichiometry(Molecule(str(moiety_stoich), 
        NUM2))

  def testCountMoietys(self):
    if IGNORE_TEST:
      return
    moiety_stoich = MoietyStoichiometry(Molecule(MOIETY_NAME1),
        NUM1)
    mole_stoich = MoleculeStoichiometry(Molecule(str(moiety_stoich)), 
        NUM2)
    df = mole_stoich.countMoietys()
    self.assertEquals(df.columns.tolist(), [cn.VALUE])
    expected = NUM1 * NUM2
    trues = [expected ==  n for n in df[cn.VALUE]]

  def testInitialize(self):
    if IGNORE_TEST:
      return
    Molecule.initialize(self.simple)
    self.assertEqual(len(Molecule.molecules), cn.NUM_SPECIES)
    

if __name__ == '__main__':
  unittest.main()
