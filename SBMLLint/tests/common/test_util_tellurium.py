from SBMLLint.common import constants as cn
from SBMLLint.common import util_tellurium
import tesbml


import numpy as np
import os
import unittest


IGNORE_TEST = False


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testGetModelFromAntimony(self):
    document = util_tellurium.getSBMLFromAntimony(
        util_tellurium.ANTIMONY_STG)
    self.assertTrue(isinstance(document, str))

  def testMakeSBMLFile(self):
    util_tellurium.makeSBMLFile()
    self.assertTrue(os.path.isfile(cn.TEST_FILE2))



if __name__ == '__main__':
  unittest.main()
