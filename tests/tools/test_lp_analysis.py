from SBMLLint.common import constants as cn
from SBMLLint.tools import lp_analysis
from SBMLLint.common import simple_sbml


import numpy as np
import os
import sys
import unittest
import yaml


IGNORE_TEST = False
TEST_SBML_INCONSISTENT_PTH = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000147_url.xml")
TEST_SBML_CONSISTENT_PTH = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000145_url.xml")
TEST_SBML_CONSISTENT_PTH = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000010_url.xml")

#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testLPAnalysis(self):
    if IGNORE_TEST:
      return
    def test(path, expected):
      with open(path, "r") as fd:
        self.assertEqual(
            lp_analysis.LPAnalysis(fd, is_report=True), expected)
    #
    test(TEST_SBML_INCONSISTENT_PTH, False)
    test(TEST_SBML_CONSISTENT_PTH, True)
    

if __name__ == '__main__':
  unittest.main()
