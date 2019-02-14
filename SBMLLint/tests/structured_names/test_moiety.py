"""
Tests for Moiety and MoietyStoichiometry
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml
from SBMLLint.structured_names.moiety  \
    import Moiety, MoietyStoichiometry, \
    _extractFromMoietyStoichiometryString

import itertools
import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = True
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
MOIETY_NAME3 = "third"
NAMES = [MOIETY_NAME1, MOIETY_NAME2, MOIETY_NAME3]
NUM1 = 1.5
NUM2 = 3
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
MOIETY_STOICHIOMETRY_STGS = {
    ("P", 1): ["P", "P_1", "P_1.0"],
    ("PP", 2): ["PP_2", "PP_2.0"],
    }
MOLECULE_STOICHIOMETRY_STGS = {
    ((Moiety("A"), 1), (Moiety("P"), 1)): 
        [Molecule("A_P"), Molecule("A_1__P"), Molecule("A__P")],
    ((Moiety("AA"), 1), (Moiety("PP"), 1)): 
        [Molecule("AA_PP"), Molecule("AA_1__PP"), Molecule("AA__PP"),
        Molecule("AA__PP_1.0")],
    }


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testExtractFromMoietyStoichiometryString(self):
    if IGNORE_TEST:
      return
    for expected, strings in MOIETY_STOICHIOMETRY_STGS.items():
      for stg in strings:
        result = _extractFromMoietyStoichiometryString(stg)
        self.assertEqual(result, expected)


class TestMoiety(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(Moiety(MOIETY_NAME1).name, MOIETY_NAME1)

  def testAppendToMolecule(self):
    if IGNORE_TEST:
      return
    moiety2 = Moiety(MOIETY_NAME2)
    molecule = Molecule(MOIETY_NAME1)
    new_molecule = moiety2.appendToMolecule(molecule)
    self.assertEqual(new_molecule.name, Molecule(MOLECULE_NAME).name)

  def testExtract(self):
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

  def testIsInMolecule(self):
    if IGNORE_TEST:
      return
    moiety = Moiety(MOIETY_NAME1)
    molecule = Molecule(MOLECULE_NAME)
    self.assertTrue(moiety.isInMolecule(molecule))


class TestMoietyStoichiometry(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)

  # TODO: test stoichiometry
  def testExtract(self):
    if IGNORE_TEST:
      return
    NUM1 = 2
    NUM1 = 1.5
    moiety1 = Moiety(MOIETY_NAME1)
    moiety_stoich1 = MoietyStoichiometry(moiety1, NUM1)
    moiety2 = Moiety(MOIETY_NAME2)
    moiety_stoich2 = MoietyStoichiometry(moiety2, NUM2)
    molecule = Molecule(MOIETY_NAME1)
    molecule = moiety2.appendToMolecule(molecule)
    import pdb; pdb.set_trace()
    names = set([m.molecule.name for m in Moiety.extract(molecule)])
    self.assertEqual(names, set([MOIETY_NAME1, MOIETY_NAME2]))

  # TODO: Update tests
  def testCountMoietys(self):
    if IGNORE_TEST:
      return
    molecule = Molecule(MOLECULE_NAME)
    df = MoietyStoichiometry.countMoietys([molecule])
    self.assertEquals(df.columns.tolist(), [cn.VALUE])
    df2 = MoietyStoichiometry.countMoietys([molecule, molecule])
    self.assertTrue(df2.equals(df + df))


if __name__ == '__main__':
  unittest.main()
