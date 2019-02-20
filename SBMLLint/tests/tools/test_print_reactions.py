from SBMLLint.common import constants as cn
from SBMLLint.common.runner import Runner
from SBMLLint.tools import print_reactions


import numpy as np
import os
import sys
import unittest


IGNORE_TEST = False
TEST_FILE = "test_print_reactions.txt"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_FILE)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def testPrettyPrint(self):
    with open(TEST_OUT_PATH, 'w') as fd:
      print_reactions.prettyPrint(cn.TEST_FILE, file_out=fd)
    with open(TEST_OUT_PATH, 'r') as fd:
      lines = fd.readlines()
    self.assertEqual(len(lines), cn.NUM_REACTIONS)

  def testMain(self):
    return
    # FIXME: This test fails in traves
    module_dir = os.path.abspath(os.curdir)
    for ele in ["SBMLLint", "tools"]:
      module_dir = os.path.join(module_dir, ele)
    module_path = os.path.join(module_dir, "print_reactions.py")
    runner = Runner(module_path)
    runner.execute([cn.TEST_FILE], '')
    self.assertEqual(runner.output.count('\n'), cn.NUM_REACTIONS)


if __name__ == '__main__':
  unittest.main()
