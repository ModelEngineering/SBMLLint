from SBMLLint.common import constants as cn
from SBMLLint.tools import lp_analysis
from SBMLLint.common import simple_sbml


import numpy as np
import os
import sys
import unittest
import yaml


IGNORE_TEST = False
TEST_SBML_PTH = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000147_url.xml")

#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testLPAnalysis(self):
    if IGNORE_TEST:
      return
    simple = simple_sbml.SimpleSBML()
    with open(TEST_SBML_PTH, "r") as fd:
      simple.initialize(fd)
    result = lp_analysis.LPAnalysis(simple)
    import pdb; pdb.set_trace()
    

if __name__ == '__main__':
  unittest.main()
