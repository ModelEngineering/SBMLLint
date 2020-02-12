from SBMLLint.common import config
from SBMLLint.common import constants as cn

import os
import unittest


IGNORE_TEST = True
TEST_BAD_CONFIG_FILE = os.path.join(cn.TEST_DIR,
    "test_cfg_bad.yml")


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def tearDown(self):
    config.setConfiguration()

  def verifyDefaultConfiguration(self):
    result = config.getConfiguration()
    self.assertTrue(isinstance(result, dict))
    for k, v in cn.CFG_DEFAULTS.items():
      self.assertTrue(k in result)
      self.assertEqual(result[k], cn.CFG_DEFAULTS[k])

  def testGetConfiguration(self):
    if IGNORE_TEST:
      return
    self.verifyDefaultConfiguration()

  def testSetConfiguration(self):
    if IGNORE_TEST:
      return
    config._config_dict = {}
    self.assertEqual(len(config.getConfiguration()), 0)
    config.setConfiguration()
    self.verifyDefaultConfiguration()

  def testSetConfiguration2(self):
    # TESTING
    with self.assertRaises(SystemExit):
      config.setConfiguration(path=TEST_BAD_CONFIG_FILE)
      
    

if __name__ == '__main__':
  unittest.main()
