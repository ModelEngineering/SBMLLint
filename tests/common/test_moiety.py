"""
Tests for Moiety and MoietyStoichiometry
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.moiety import Moiety, MoietyStoichiometry
from SBMLLint.common import util

import itertools
import numpy as np
import os
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

  def testGetMoietys(self):
    m_ss = []
    for stgs in MOIETY_STOICHIOMETRY_STGS.values():
      m_ss.extend([MoietyStoichiometry.make(s) for s in stgs])
    moietys = MoietyStoichiometry.getMoietys(m_ss)
    self.assertEqual(len(MOIETY_STOICHIOMETRY_STGS.keys()),
        len(moietys))


if __name__ == '__main__':
  unittest.main()
