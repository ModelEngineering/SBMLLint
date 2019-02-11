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
      sbmllint.lint(cn.TEST_FILE4, file_out=fd)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertGreater(len(lines), 0)

  def testMain(self):
    module_dir = os.path.abspath(os.curdir)
    for ele in ["SBMLLint", "tools"]:
      module_dir = os.path.join(module_dir, ele)
    module_path = os.path.join(module_dir, "sbmllint.py")
    runner = Runner(module_path)
    runner.execute([cn.TEST_FILE4], '')
    self.assertGreater(runner.output.count('\n'), 0)


if __name__ == '__main__':
  unittest.main()
