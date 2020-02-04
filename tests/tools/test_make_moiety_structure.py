from SBMLLint.common import constants as cn
from SBMLLint.tools import make_moiety_structure


import numpy as np
import os
import sys
import unittest


IGNORE_TEST = False
TEST_MOIETY_FILE = os.path.join(cn.TEST_DIR,
    "test_moieties.yml")
TEST_CFG_FILE = os.path.join(cn.TEST_DIR,
     "test_cfg_file.yml")
TEST_SBML_FILE = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000147_url.xml")


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def _cleanUp(self):
    if os.path.isfile(TEST_CFG_FILE):
      os.remove(TEST_CFG_FILE)

  def setUp(self):
    self._cleanUp()
    self.xml_fid = open(TEST_SBML_FILE)
    self.moiety_fid = open(TEST_MOIETY_FILE)

  def tearDown(self):
    self._cleanUp()
    self.xml_fid.close()
    self.moiety_fid.close()

  # TODO: 1. Test if substring
  def testGetMoieties(self):
    names = make_moiety_structure.getMoieties(self.moiety_fid)

  def testMain(self):
    pass


if __name__ == '__main__':
  unittest.main()
