"""
Tests for MoietyComparator
"""
from SBMLLint.common import constants as cn
from SBMLLint.common import exceptions
from SBMLLint.common import config
from SBMLLint.common.molecule import Molecule, MoleculeStoichiometry
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.common import util
from SBMLLint.moiety_analysis.moiety_comparator  \
     import MoietyComparator

import copy
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
MOLECULE_NAMES = []  # A set of names from moiety combinations
for item in iterator:
  name = ""
  for idx,ele in enumerate(item):
    if ele == 1:
      if len(name) == 0:
        name = NAMES[idx]
      else:
        name = "%s%s%s" % (name, cn.MOIETY_SEPARATOR, NAMES[idx])
  if len(name) > 0:
    MOLECULE_NAMES.append(name)
MOLECULE_NAMES1 = MOLECULE_NAMES[:3]
MOLECULE_NAMES2 = MOLECULE_NAMES[3:]
TEST_FILE4 = "test_file4.xml"
TEST_FILE3 = "test_file3.antimony"
# Antimony as source
PATH= os.path.join(cn.TEST_DIR, TEST_FILE3)
with open(PATH, 'r') as fd:
  lines = fd.readlines()
NUM1 = 2
NUM2 = 3


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
        for n in MOLECULE_NAMES1]
    self.molecules2 = [MoleculeStoichiometry(Molecule(n), NUM2)
        for n in MOLECULE_NAMES2]
    self.comparator = MoietyComparator(self.molecules1,
        self.molecules2)
    self.config_dict = copy.deepcopy(config._config_dict)

  def tearDown(self):
    config._config_dict = self.config_dict

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(
        len(self.comparator.molecule_stoichiometry_collections), 2)

  def testIsSame(self):
    if IGNORE_TEST:
      return
    self.assertFalse(self.comparator.isSame())
    comparator = MoietyComparator(self.molecules1, self.molecules1)
    self.assertTrue(comparator.isSame())

  def testDifference1(self):
    if IGNORE_TEST:
      return
    comparator = MoietyComparator(self.molecules1, self.molecules2)
    df = comparator.difference()
    self.assertLess(df.loc[MOIETY_NAME1].tolist()[0], 0)

  def testDifference2(self):
    if IGNORE_TEST:
      return
    config._config_dict[cn.CFG_IGNORED_MOIETIES] = [MOIETY_NAME1]
    comparator = MoietyComparator(self.molecules1,
        self.molecules2)
    df = comparator.difference()
    self.assertFalse(MOIETY_NAME1 in df.index)

  def testDifference3(self):
    if IGNORE_TEST:
      return
    config._config_dict[cn.CFG_IGNORED_MOLECULES] = [MOIETY_NAME1]
    comparator = MoietyComparator(self.molecules1,
        self.molecules2)
    df = comparator.difference()
    expected_count = sum([NUM1 if MOIETY_NAME1 in m else 0 
        for m in MOLECULE_NAMES1 if len(MOIETY_NAME1) != len(m)])
    expected_count -= sum([NUM2 if MOIETY_NAME1 in m else 0 
        for m in MOLECULE_NAMES2 if len(MOIETY_NAME1) != len(m)])
    self.assertEqual(df.loc[MOIETY_NAME1, cn.VALUE], expected_count)

  def testDifference4(self):
    if IGNORE_TEST:
      return
    config._config_dict[cn.CFG_PROCESS_BOUNDARY_REACTIONS] = False
    molecules1 = [MoleculeStoichiometry(Molecule(n), 0)
        for n in MOLECULE_NAMES[:3]]
    comparator = MoietyComparator(molecules1, self.molecules2)
    df = comparator.difference()
    self.assertEqual(df[df.columns[0]].sum(), 0)

  def testDifference5(self):
    if IGNORE_TEST:
      return
    config._config_dict[cn.CFG_PROCESS_BOUNDARY_REACTIONS] = True
    molecules1 = [MoleculeStoichiometry(Molecule(n), 0)
        for n in MOLECULE_NAMES[:3]]
    comparator = MoietyComparator(molecules1, self.molecules2)
    df = comparator.difference()
    self.assertNotEqual(df[df.columns[0]].sum(), 0)
  
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
    try:
      SBML= util.getXMLFromAntimony(''.join(lines))
    except exceptions.MissingTelluriumError:
      return
    simple.initialize(SBML)
    result = analyze(simple)
    self.assertGreaterEqual(result.num_reactions, 0)
    self.assertGreaterEqual(result.num_imbalances, 0)
    self.assertTrue('2' in result.report)
    self.assertGreater(result.report.count('\n'),  5)
    

if __name__ == '__main__':
  unittest.main()
