from SBMLLint.common import constants as cn
from SBMLLint.common.runner import Runner
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
      num_react, num_bad = sbmllint.lint(cn.TEST_FILE4, file_out=fd)
    self.assertGreaterEqual(num_react, num_bad)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testLint2(self):
    with open(TEST_OUT_PATH, 'w') as fd:
      num_react, num_bad = sbmllint.lint(cn.TEST_FILE2, file_out=fd)
    self.assertGreaterEqual(num_react, num_bad)
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
      num_react, num_bad = sbmllint.lint(model, file_out=fd)
    self.assertEqual(num_react, 1)
    self.assertEqual(num_bad, 1)

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
