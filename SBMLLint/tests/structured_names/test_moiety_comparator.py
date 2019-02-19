"""
Tests for MoietyComparator
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util
from SBMLLint.structured_names.moiety_comparator  \
     import MoietyComparator

import itertools
import numpy as np
import os
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
TEST_FILE4 = "test_file4.xml"
TEST_FILE3 = "test_file3.antimony"
# Antimony as source
PATH= os.path.join(cn.TEST_DIR, TEST_FILE3)
with open(PATH, 'r') as fd:
  lines = fd.readlines()
NUM1 = 2
NUM2 = 3
SBML= util.getXMLFromAntimony(''.join(lines))


######################################
# Auxiliary Functions
######################################
def analyze(simple):
  return MoietyComparator.analyzeReactions(simple)
  

#############################
# Tests
#############################
class TestMoietyComparator(unittest.TestCase):

  def setUp(self):
    self.molecules1 = [MoleculeStoichiometry(Molecule(n), NUM1)
        for n in MOLECULE_NAME_SET[:3]]
    self.molecules2 = [MoleculeStoichiometry(Molecule(n), NUM2)
        for n in MOLECULE_NAME_SET[3:]]
    self.comparator = MoietyComparator(self.molecules1,
        self.molecules2)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(
        len(self.comparator.molecule_stoichiometry_collections), 2)

  def testIsSame(self):
    if IGNORE_TEST:
      return
    self.assertFalse(self.comparator.isSame())
    comparator = MoietyComparator(self.molecules1,
        self.molecules1)
    self.assertTrue(comparator.isSame())

  def testDifference(self):
    if IGNORE_TEST:
      return
    df = self.comparator.difference()
    self.assertLess(df.loc[MOIETY_NAME1].tolist()[0], 0)
  
  def testReportDifference(self):  
    if IGNORE_TEST:
      return
    stg = self.comparator.reportDifference()
    self.assertGreater(len(stg), 0)
    comparator = MoietyComparator(self.molecules1,
        self.molecules1)
    stg = comparator.reportDifference()
    self.assertEqual(len(stg), 0)

  def testAnalyzeReactions(self):
    if IGNORE_TEST:
      return
    simple = SimpleSBML()
    simple.initialize(SBML)
    num, stg = analyze(simple)
    self.assertGreaterEqual(num, 0)
    self.assertTrue('2' in stg)
    self.assertGreater(stg.count('\n'),  5)
    

if __name__ == '__main__':
  unittest.main()
