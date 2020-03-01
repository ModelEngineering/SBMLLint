from SBMLLint.common import constants as cn
from SBMLLint.common import exceptions
from SBMLLint.common.runner import Runner
from SBMLLint.common.simple_sbml import SimpleSBML
from SBMLLint.tools import sbmllint


import numpy as np
import os
import sys
import unittest


IGNORE_TEST = False
IS_REPORT = False
TEST_FILE = "test_sbmllint.txt"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_FILE)
TEST_147_CFG_FILE2 = os.path.join(cn.TEST_DIR,
    "test_BIOMOD147_cfg2.yml")
TEST_147_CFG_FILE = os.path.join(cn.TEST_DIR,
    "test_BIOMOD147_cfg.yml")
TEST_147_SBML_FILE = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000147_url.xml")


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def testLint(self):
    if IGNORE_TEST:
      return
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=cn.TEST_FILE4, 
          file_out=fd,
          is_report=IS_REPORT,
          mass_balance_check=cn.MOIETY_ANALYSIS)
    self.assertGreaterEqual(
        result.num_reactions, result.num_imbalances)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLintWithXMLFileFid(self):
    if IGNORE_TEST:
      return
    fid = open(cn.TEST_FILE4, "r")
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=fid, file_out=fd,
          is_report=IS_REPORT, mass_balance_check=cn.MOIETY_ANALYSIS)
    fid.close()
    self.assertGreater(len(result), 0)

  def testLintWithConfigFid(self):
    if IGNORE_TEST:
      return
    def get(config_fid=None):
      with open(TEST_OUT_PATH, 'w') as fd:
        result = sbmllint.lint(model_reference=cn.TEST_FILE4, file_out=fd,
            config_fid=config_fid,
            is_report=IS_REPORT,
            mass_balance_check=cn.MOIETY_ANALYSIS)
      return result
    #
    result1 = get()
    fid = open(cn.CFG_DEFAULT_PATH)
    result2 = get(config_fid=fid)
    self.assertEqual(result1, result2)
    fid.close()

  def testLint(self):
    if IGNORE_TEST:
      return
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=cn.TEST_FILE4,
          file_out=fd,
          is_report=True,
          mass_balance_check=cn.MOIETY_ANALYSIS)
    self.assertGreaterEqual(
        result.num_reactions, result.num_imbalances)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLint2(self):
    if IGNORE_TEST:
      return
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=cn.TEST_FILE2,
          file_out=fd, is_report=True,
          mass_balance_check=cn.MOIETY_ANALYSIS)
    self.assertGreaterEqual(
        result.num_reactions, result.num_imbalances)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLint3(self):
    if IGNORE_TEST:
      return
    model = """
    2Glu + 2A_P_P_P -> 2Glu_P + 2A_P_P; 1
    Glu = 0
    A_P_P_P = 0
    Glu_P = 0
    A_P_P = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      try:
        result = sbmllint.lint(model_reference=model, file_out=fd,
            is_report=IS_REPORT,
            mass_balance_check=cn.MOIETY_ANALYSIS)
      except exceptions.MissingTelluriumError:
        return
    self.assertEqual(result.num_reactions, 1)
    self.assertEqual(result.num_imbalances, 0)

  def testLint4(self):
    if IGNORE_TEST:
      return
    model = """
    2Glu_DUMMYMOIETY + 2A__P_3 -> 2Glu_P + 2A_P_P; 1
    Glu_DUMMYMOIETY = 0
    A_P_P = 0
    Glu_P = 0
    A__P_3 = 0
    A = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      try:
        result = sbmllint.lint(model_reference=model, file_out=fd,
            is_report=IS_REPORT,
            mass_balance_check=cn.MOIETY_ANALYSIS)
      except exceptions.MissingTelluriumError:
        return
    self.assertEqual(result.num_reactions, 1)
    self.assertEqual(result.num_imbalances, 0)

  def testLint5(self):
    if IGNORE_TEST:
      return
    model = """
    A -> B; 1
    B -> B + DUMMYMOLECULE;1
    A = 0
    B = 0
    DUMMYMOLECULE = 0
    """
    with open(TEST_OUT_PATH, 'w') as fd:
      try:
        result = sbmllint.lint(model_reference=model,
                               file_out=fd,
                               mass_balance_check="games",
                               implicit_games=False,
                               is_report=IS_REPORT
                               )
      except exceptions.MissingTelluriumError:
        return
    self.assertTrue(result)
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=model,
                             file_out=fd,
                             mass_balance_check="games",
                             implicit_games=True,
                             is_report=IS_REPORT
                             )    
    self.assertFalse(result)

  def testRemoveIgnored(self):
    if IGNORE_TEST:
      return
    implicit = "MA"
    path = os.path.join(cn.BIOMODELS_DIR, cn.TEST_FILE13)
    simple = SimpleSBML()
    simple.initialize(path)
    simple = sbmllint.removeIgnored(simple, implicit)
    implicit_reactions = []
    for r in simple.reactions:
      reactants = [reactant.molecule.name for reactant in r.reactants]
      products = [product.molecule.name for product in r.products]
      if (implicit in reactants) or (implicit in products):
        implicit_reactions.append(r.label)
    self.assertTrue(len(implicit_reactions) == 0)

  def testMain(self):
    if IGNORE_TEST:
      return
    return
    # FIXME: This test fails in traves
    module_dir = os.path.abspath(os.curdir)
    for ele in ["SBMLLint", "tools"]:
      module_dir = os.path.join(module_dir, ele)
    module_path = os.path.join(module_dir, "sbmllint.py")
    runner = Runner(module_path)
    runner.execute([cn.TEST_FILE4], '')
    self.assertGreater(runner.output.count('\n'), 0)

  def testBIOMOD147(self):
    if IGNORE_TEST:
      return
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=TEST_147_SBML_FILE,
          config_fid=open(TEST_147_CFG_FILE, "r"),
          file_out=fd,
          mass_balance_check=cn.MOIETY_ANALYSIS,
          implicit_games=False,
          is_report=IS_REPORT,
          )
    # Verify that exclicit declaration is decomposed
    # into moieties
    molecule_name = "IkBeIKKNFkB"
    self.assertFalse("%s:" % molecule_name in result.report)

  def testBIOMOD147_2(self):
    if IGNORE_TEST:
      return
    with open(TEST_OUT_PATH, 'w') as fd:
      result = sbmllint.lint(model_reference=TEST_147_SBML_FILE,
          config_fid=open(TEST_147_CFG_FILE2, "r"),
          file_out=fd,
          mass_balance_check=cn.MOIETY_ANALYSIS,
          implicit_games=False,
          is_report=IS_REPORT,
          )
    # Verify that exclicit declaration is decomposed
    # into moieties
    self.assertFalse("cytoplasm:" in result.report)
    self.assertFalse("nulceus:" in result.report)


if __name__ == '__main__':
  unittest.main()
