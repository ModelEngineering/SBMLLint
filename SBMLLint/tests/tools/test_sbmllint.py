from SBMLLint.common import constants as cn
from SBMLLint.common.runner import Runner
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.tools import sbmllint


import numpy as np
import os
import sys
import unittest


IGNORE_TEST = False
TEST_FILE = "test_sbmllint.txt"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_FILE)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def testLint(self):
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(cn.TEST_FILE4, file_out=fd,
          mass_balance_check=sbmllint.STRUCTURED_NAMES)
    self.assertGreaterEqual(
        result.num_reactions, result.num_imbalances)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLint2(self):
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(cn.TEST_FILE2, file_out=fd,
          mass_balance_check=sbmllint.STRUCTURED_NAMES)
    self.assertGreaterEqual(
        result.num_reactions, result.num_imbalances)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLint3(self):
    model = """
    2Glu + 2A_P_P_P -> 2Glu_P + 2A_P_P; 1
    Glu = 0
    A_P_P_P = 0
    Glu_P = 0
    A_P_P = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model, file_out=fd,
          mass_balance_check=sbmllint.STRUCTURED_NAMES)
    self.assertEqual(result.num_reactions, 1)
    self.assertEqual(result.num_imbalances, 0)

  def testLint4(self):
    model = """
    2Glu_DUMMYIMPLICIT + 2A__P_3 -> 2Glu_P + 2A_P_P; 1
    Glu_DUMMYIMPLICIT = 0
    A_P_P = 0
    Glu_P = 0
    A__P_3 = 0
    A = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model, file_out=fd,
          mass_balance_check=sbmllint.STRUCTURED_NAMES)
    self.assertEqual(result.num_reactions, 1)
    self.assertEqual(result.num_imbalances, 0)

  def testLint5(self):
    model = """
    A -> B; 1
    B -> B + DUMMYIMPLICIT;1
    A = 0
    B = 0
    DUMMYIMPLICIT = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model,
                             file_out=fd,
                             mass_balance_check="games",
                             implicit_games=False
                             )
    self.assertTrue(result)
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model,
                             file_out=fd,
                             mass_balance_check="games",
                             implicit_games=True
                             )    
    self.assertFalse(result)

  def testRemoveImplicit(self):
    implicit = "MA"
    path = os.path.join(cn.BIOMODELS_DIR, cn.TEST_FILE13)
    simple = SimpleSBML()
    simple.initialize(path)
    simple = sbmllint.removeImplicit(simple, implicit)
    implicit_reactions = []
    for r in simple.reactions:
      reactants = [reactant.molecule.name for reactant in r.reactants]
      products = [product.molecule.name for product in r.products]
      if (implicit in reactants) or (implicit in products):
        implicit_reactions.append(r.label)
    self.assertTrue(len(implicit_reactions) == 0)

  def testMain(self):
    return
    # FIXME: This test fails in traves
    module_dir = os.path.abspath(os.curdir)
    for ele in ["SBMLLint", "tools"]:
      module_dir = os.path.join(module_dir, ele)
    module_path = os.path.join(module_dir, "sbmllint.py")
    runner = Runner(module_path)
    runner.execute([cn.TEST_FILE4], '')
    self.assertGreater(runner.output.count('\n'), 0)


if __name__ == '__main__':
  unittest.main()
