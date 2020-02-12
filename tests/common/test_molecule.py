"""
Tests for Molecule
"""
from SBMLLint.common import constants as cn
from SBMLLint.common import config
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML

import itertools
import numpy as np
import os
import unittest


IGNORE_TEST = False
NUM1 = 2
NUM2 = 3
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
MOIETY_NAME3 = "third"
NAME = "name"
NO_NAME = 'not_a_name'
MOLECULE_NAME = "%s%s%s" % (MOIETY_NAME1, cn.MOIETY_DOUBLE_SEPARATOR, 
    MOIETY_NAME2)
MOLECULE_STOICHIOMETRY_STGS = {
    ((Moiety("A"), 1), (Moiety("P"), 1)): 
        [Molecule("A_P"), Molecule("A_1__P"), Molecule("A__P")],
    ((Moiety("AA"), 1), (Moiety("PP"), 1)): 
        [Molecule("AA_PP"), Molecule("AA_1__PP"), Molecule("AA__PP"),
        ],
    }
NAMES = [MOIETY_NAME1, MOIETY_NAME2, MOIETY_NAME3]
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
TEST_CFG_FILE = os.path.join(cn.TEST_DIR, "test_sbmllint_cfg.yml")


#############################
# Tests
#############################
class TestMolecule(unittest.TestCase):

  def _init(self):
    Molecule.molecules = []
    self.simple = SimpleSBML()
    self.simple = self.simple.initialize(cn.TEST_FILE)

  def setUp(self):
    if IGNORE_TEST:
      return
    self._init()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(NAME)
    self.assertEqual(molecule.name, NAME)

  def testAppend(self):
    if IGNORE_TEST:
      return
    moiety2 = Moiety(MOIETY_NAME2)
    molecule = Molecule(MOIETY_NAME1)
    new_molecule = molecule.append(moiety2)
    self.assertEqual(new_molecule.name, Molecule(MOLECULE_NAME).name)

  def testGetMoietys(self):
    if IGNORE_TEST:
      return
    m_s1 = MoietyStoichiometry(Moiety(MOIETY_NAME1), NUM1)
    m_s2 = MoietyStoichiometry(MOIETY_NAME2, NUM2)
    m_s3 = MoietyStoichiometry(MOIETY_NAME1, NUM2)
    molecule = Molecule(str(m_s1))
    molecule = molecule.append(Moiety(str(m_s2)))
    moietys = molecule.getMoietys()
    expected = [Moiety(MOIETY_NAME1), Moiety(MOIETY_NAME2)]
    for moiety in moietys:
      self.assertTrue(any([moiety.isEqual(e) for e in expected]))

  def testHasMoiety(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(MOLECULE_NAME)
    self.assertTrue(molecule.hasMoiety(Moiety(MOIETY_NAME1)))
    self.assertTrue(molecule.hasMoiety(Moiety(MOIETY_NAME2)))

  def testMoietyStoichiometrys(self):
    if IGNORE_TEST:
      return
    config.setConfiguration(TEST_CFG_FILE)
    #
    molecule = Molecule("Glu6P")
    self.assertEqual(molecule.moiety_stoichiometrys[0].moiety.name,
        "Glu")
    #
    molecule = Molecule("ATP")
    self.assertEqual(len(molecule.moiety_stoichiometrys), 2)


class TestMoleculeStoichiometry(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML()
    self.simple = self.simple.initialize(cn.TEST_FILE)
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
        result = Molecule.getMoietyStoichiometrys(stg)
        trues = [x.isEqual(y) for x, y in zip(result, m_ss)]
        self.assertTrue(all(trues))

  def testExtractMoietyStoichiometrys(self):
    if IGNORE_TEST:
      return
    moiety_stoich = MoietyStoichiometry(Molecule(MOIETY_NAME1),
        NUM1)
    mole_stoich = MoleculeStoichiometry(Molecule(str(moiety_stoich)), 
        NUM2)

  def testCountMoietys(self):
    if IGNORE_TEST:
      return
    moiety_stoich = MoietyStoichiometry(Molecule(MOIETY_NAME1),
        NUM1)
    mole_stoich = MoleculeStoichiometry(Molecule(str(moiety_stoich)), 
        NUM2)
    df = mole_stoich.countMoietys()
    self.assertEqual(df.columns.tolist(), [cn.VALUE])
    expected = NUM1 * NUM2
    trues = [expected ==  n for n in df[cn.VALUE]]

  def testCountMoietys2(self):
    if IGNORE_TEST:
      return
    m_ss = [
        MoleculeStoichiometry(Molecule("A_P_P_P"), 1),
        MoleculeStoichiometry(Molecule("A__P_3"), 1),
        ]
    dfs = []
    for m_s in m_ss:
      dfs.append(m_s.countMoietys())
    self.assertTrue(dfs[0].equals(dfs[1]))

  def testCountMoietysInCollection(self):
    if IGNORE_TEST:
      return
    m_ss = [MoleculeStoichiometry(Molecule(n), NUM1)
        for n in MOLECULE_NAME_SET[:3]]
    df = MoleculeStoichiometry.countMoietysInCollection(m_ss)
    indexes = df.index.tolist()
    self.assertEqual(len(indexes), len(set(indexes)))

  def testGetMolecules(self):
    if IGNORE_TEST:
      return
    names = ["a", "b", "c"]
    full_names = list(names)
    full_names.extend(names)
    m_ss = [MoleculeStoichiometry(Molecule(n), NUM1)
        for n in full_names]
    molecules = MoleculeStoichiometry.getMolecules(m_ss)
    self.assertEqual(len(molecules), len(names))
    

if __name__ == '__main__':
  unittest.main()
