from SBMLLint.common import config

import unittest

IGNORE_TEST = False


#############################
# Tests
#############################
class TestFunctions(unittest.TestCase):

  def testGetConfiguration(self):
    result = config.getConfiguration()
    self.assertTrue(isinstance(result, dict))
    self.assertTrue("implicits" in result.keys())
    

if __name__ == '__main__':
  unittest.main()
