from SBMLLint.common import config
from SBMLLint.common import constants as cn

import unittest


IGNORE_TEST = False


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
    self.verifyDefaultConfiguration()

  def testSetConfiguration(self):
    config._config_dict = {}
    self.assertEqual(len(config.getConfiguration()), 0)
    config.setConfiguration()
    self.verifyDefaultConfiguration()
    

if __name__ == '__main__':
  unittest.main()
