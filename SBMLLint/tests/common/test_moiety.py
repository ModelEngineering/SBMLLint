"""
Tests for Moiety and MoietyStoichiometry
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import simple_sbml
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry

import itertools
import numpy as np
import os
import tesbml
import unittest


IGNORE_TEST = False
MOIETY_NAME1 = "first"
MOIETY_NAME2 = "second"
NUM1 = 2
MOIETY_STOICHIOMETRY_STGS = {
    ("P", 1): ["P", "P_1"],
    ("PP", 2): ["PP_2"],
    }

#######################################
class TestMoiety(unittest.TestCase):

  def setUp(self):
    self.simple = SimpleSBML(cn.TEST_FILE)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(Moiety(MOIETY_NAME1).name, MOIETY_NAME1)


class TestMoietyStoichiometry(unittest.TestCase):

  def testStoichiometry(self):
    moiety_stoichiometry = MoietyStoichiometry(Moiety(MOIETY_NAME1),
        NUM1)
    self.assertEqual(moiety_stoichiometry.moiety.name,
        MOIETY_NAME1)
    self.assertEqual(moiety_stoichiometry.stoichiometry, NUM1)
 

  def testMake(self):
    if IGNORE_TEST:
      return
    for expected, strings in MOIETY_STOICHIOMETRY_STGS.items():
      for stg in strings:
        result = MoietyStoichiometry.make(stg)
        self.assertEqual(result.moiety.name, expected[0])
        self.assertEqual(result.stoichiometry, expected[1])


if __name__ == '__main__':
  unittest.main()
