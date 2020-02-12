from SBMLLint.common import constants as cn
from SBMLLint.common.runner import Runner
from SBMLLint.tools import analyze_moiety_analysis


import numpy as np
import os
import unittest


IGNORE_TEST = False
TEST_FILE = "test_analyze_moiety_analysis.csv"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_FILE)


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def testIsStructuredName(self):
    self.assertTrue(analyze_moiety_analysis.isStructuredName("a_b"))
    self.assertFalse(analyze_moiety_analysis.isStructuredName("a_1"))
    self.assertFalse(analyze_moiety_analysis.isStructuredName("Species_1"))
    self.assertFalse(analyze_moiety_analysis.isStructuredName("ATP"))

  def testCalcStats(self):
    analyze_moiety_analysis.calcStats(initial=1, final=5,
         out_path=TEST_OUT_PATH, report_progress=False)
    self.assertTrue(os.path.isfile(TEST_OUT_PATH))



if __name__ == '__main__':
  unittest.main()
