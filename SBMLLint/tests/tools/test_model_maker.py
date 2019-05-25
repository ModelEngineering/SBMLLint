from SBMLLint.common import constants as cn
from SBMLLint.tools.model_maker import ModelMaker

import os
import unittest


IGNORE_TEST = False
TEST_IN_FILE = "test_model_maker.txt"
TEST_IN_PATH = os.path.join(cn.TEST_DIR, TEST_IN_FILE)
TEST_OUT_FILE = "test_model_maker.ant"
TEST_OUT_PATH = os.path.join(cn.TEST_DIR, TEST_OUT_FILE)
SYM1 = "sym1"
SYM2 = "sym2"
SYM3 = "sym3"
SYM4 = "sym4"
REACTIONSTR = "%s + %s -> %s + %s" % (SYM1, SYM2, SYM3, SYM4)


#############################
# Tests
#############################
class TestModelMaker(unittest.TestCase):

  def tearDown(self):
    if os.path.isfile(TEST_OUT_PATH):
      os.remove(TEST_OUT_PATH)

  def setUp(self):
    self.maker = ModelMaker(REACTIONSTR)

  def testConstructor(self):
    self.assertEqual(self.maker._reaction_strs, [REACTIONSTR])
    self.maker = ModelMaker([REACTIONSTR])
    self.assertEqual(self.maker._reaction_strs, [REACTIONSTR])
    self.maker = ModelMaker(TEST_IN_PATH)
    self.assertGreater(len(self.maker._reaction_strs), 0)

  def testExtractSymbols(self):
    symbols = self.maker.extractSymbols()
    b = set(symbols).symmetric_difference([SYM1, SYM2, SYM3, SYM4])
    self.assertEqual(len(b), 0)
    self.maker = ModelMaker([REACTIONSTR, REACTIONSTR])
    b = set(symbols).symmetric_difference([SYM1, SYM2, SYM3, SYM4])
    self.assertEqual(len(b), 0)

  def testMakeModelStr(self):
    model_str = self.maker.makeModelStr()
    self.assertGreater(len(model_str), len(REACTIONSTR))
   
  def testMakeModelFile(self): 
    model_str = self.maker.makeModelFile(TEST_OUT_PATH)


if __name__ == '__main__':
  unittest.main()
