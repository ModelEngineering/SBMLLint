from SBMLLint.common import constants as cn
from SBMLLint.tools import make_moiety_structure


import numpy as np
import os
import sys
import unittest
import yaml


IGNORE_TEST = False
TEST_MOIETY_FILE = os.path.join(cn.TEST_DIR,
    "test_moieties.yml")
TEST_MOIETY_REDUN_FILE = os.path.join(cn.TEST_DIR,
    "test_moieties_redun.yml")
# moiety_structure section of config file
TEST_CONFIG_FILE = os.path.join(cn.TEST_DIR,
     "test_cfg_file.yml")
TEST_SBML_FILE = os.path.join(cn.TEST_DIR,
    "test_BIOMD0000000147_url.xml")
MOIETIES = ["IkBa", "IkBb", "IkBe", "IKK", "NFkB",
    "nucleus", "cytoplasm", "mRNA"]
MOIETIES_REDUN = ["A", "AA"]


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def _cleanUp(self):
    for ffile in [TEST_CONFIG_FILE, TEST_MOIETY_FILE,
        TEST_MOIETY_REDUN_FILE]:
      if os.path.isfile(ffile):
        os.remove(ffile)

  def setUp(self):
    self._cleanUp()
    with open(TEST_MOIETY_FILE, "w") as fd:
      yaml.dump(MOIETIES, fd)
    self.xml_fid = open(TEST_SBML_FILE)
    self.moiety_fid = open(TEST_MOIETY_FILE)

  def tearDown(self):
    self._cleanUp()
    self.xml_fid.close()
    self.moiety_fid.close()

  def testGetMoieties(self):
    names = make_moiety_structure.getMoieties(self.moiety_fid)
    diff = set(names).symmetric_difference(MOIETIES)
    self.assertEqual(len(diff), 0)

  def testGetMoietiesSubstrings(self):
    with open(TEST_MOIETY_REDUN_FILE, "w") as fd:
      yaml.dump(MOIETIES_REDUN, fd)
    with open(TEST_MOIETY_REDUN_FILE, "r") as fd:
      try:
        names = make_moiety_structure.getMoieties(fd)
        self.assertTrue(False)
      except ValueError:
        pass
      except:
        self.assertTrue(False)

  def testMain(self):
    def isSubstr(stg, stgs):
      return any([stg in s for s in stgs])
    #
    make_moiety_structure.main(self.xml_fid,
        self.moiety_fid, TEST_CONFIG_FILE)
    self.assertTrue(os.path.isfile(TEST_CONFIG_FILE))
    with open(TEST_CONFIG_FILE, "r") as fd:
      dct = yaml.load(fd)
    moiety_dct = dct[cn.CFG_MOIETY_STRUCTURE]
    missing = [m for m in MOIETIES 
        if not isSubstr(m, moiety_dct.keys())]
    self.assertEqual(len(missing), 0)
    

if __name__ == '__main__':
  unittest.main()
