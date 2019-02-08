"""
Tests for MoietyComparator
"""
from SBMLLint.common import constants as cn
from SBMLLint.common.molecule import Molecule
from SBMLLint.common.reaction import Reaction
from SBMLLint.common.simple_sbml import SimpleSBML
import SBMLLint.common.util_tellurium as util_tellurium
from SBMLLint.structured_names.moiety_comparator  \
     import MoietyComparator

import itertools
import numpy as np
import os
import unittest


IGNORE_TEST = True
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


######################################
# Auxiliary Functions
######################################
def analyze(simple):
  Molecule.initialize(simple)
  Reaction.initialize(simple)
  return MoietyComparator.analyzeReactions()
  

#############################
# Tests
#############################
class TestMoietyComparator(unittest.TestCase):

  def setUp(self):
    self.molecules1 = [Molecule(n) for n in MOLECULE_NAME_SET[:3]]
    self.molecules2 = [Molecule(n) for n in MOLECULE_NAME_SET[4:]]
    self.comparator = MoietyComparator(self.molecules1,
        self.molecules2)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.comparator.list_of_molecules), 2)

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
    stg = self.comparator.reportDifference()
    self.assertGreater(len(stg), 0)
    comparator = MoietyComparator(self.molecules1,
        self.molecules1)
    stg = comparator.reportDifference()
    self.assertEqual(len(stg), 0)


  def testAnalyzeReactions1(self):
    path = os.path.join(cn.TEST_DIR, TEST_FILE3)
    with open(path, 'r') as fd:
      lines = fd.readlines()
    sbml = util_tellurium.getSBMLStringFromAntimony(''.join(lines))
    simple = SimpleSBML(sbml)
    stg = analyze(simple)
    import pdb; pdb.set_trace()

  def testAnalyzeReactions2(self):
    path = os.path.join(cn.TEST_DIR, TEST_FILE4)
    simple = SimpleSBML(path)
    stg = analyze(simple)
    self.assertTrue('10' in stg)
    self.assertGreater(stg.count('\n'),  5)
    

if __name__ == '__main__':
  unittest.main()
